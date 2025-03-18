use crate::utils::{Direction, MethPos, NumberOfNucleotides, PosType, Position};
use anyhow::{Context, Result};

use rust_htslib::bam::Record;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use std::io::prelude::*;

pub fn process_bedgraph_data(bed_graph_path: PathBuf) -> Result<Vec<MethPos>> {
    let bedgraph_file =
        File::open(bed_graph_path).with_context(|| "Unable to open bedGraph file.")?;
    let bedgraph_reader = BufReader::new(bedgraph_file);
    let mut meth_positions: Vec<MethPos> = Vec::new();

    for bed_line in bedgraph_reader.lines() {
        let bed_line = bed_line?;
        if bed_line.starts_with("track") {
            continue;
        }
        let bed_fields: Vec<&str> = bed_line.split('\t').collect();

        let chrom = bed_fields[0]
            .strip_prefix("chr")
            .unwrap_or(bed_fields[0])
            .to_string();

        let pos = bed_fields[1]
            .parse::<usize>()
            .context("Invalid position value")?;
        let methylation: i64 = bed_fields[3]
            .parse::<f64>() // Erst in f64 umwandeln
            .context("Invalid methylation value")?
            .round() as i64;
        let position = MethPos::new(Position::new(chrom.clone(), pos as u32), methylation);

        meth_positions.push(position);
    }
    Ok(meth_positions)
}

// TODO: Test ueber unterschiedliche Chromosome hinweg
pub fn count_bases_in_reads(
    bam_file_path: PathBuf,
    mut candidate_pos: Vec<MethPos>,
) -> Result<HashMap<PosType, NumberOfNucleotides>> {
    let mut indexed_reader =
        bam::IndexedReader::from_path(bam_file_path).context("Unable to read BAM/CRAM file.")?;
    // let mut bam_reader = bam::RecordBuffer::new(indexed_bam, true);

    let mut position_counts = HashMap::new();
    let total_positions = candidate_pos.len();
    for (index, bcf_position) in candidate_pos.iter_mut().enumerate() {
        let chrom = bcf_position.position().chrom();
        let start = *bcf_position.position().pos() as i64;
        let end = start + 1;
        let tid = indexed_reader
            .header()
            .tid(chrom.as_bytes())
            .with_context(|| format!("Chromosome {} not found in BAM file header", chrom))?;
        indexed_reader.fetch((tid, start, end))?;
        for record in indexed_reader.records() {
            let record = record?;
            insert_bases(&record, bcf_position, 0);
            insert_bases(&record, bcf_position, 1);
        }
        // Fill 8 types of methylation position hashmaps
        for ((direction, cpg_position), value) in bcf_position.meth_bases().iter() {
            let pos_type = PosType::new(
                direction.clone(),
                bcf_position.methylation() > &20,
                *cpg_position,
            );
            position_counts
                .entry(pos_type)
                .and_modify(|existing_meth_pos: &mut HashMap<char, usize>| {
                    for (base, count) in value {
                        *existing_meth_pos.entry(*base).or_insert(0) += count;
                    }
                })
                .or_insert_with(HashMap::new);
        }
        if (index + 1) % 1000 == 0 || (index + 1) == total_positions {
            println!("\tProcessed {}/{} positions", index + 1, total_positions);
        }
    }
    Ok(position_counts)
}

pub fn insert_bases(record: &Record, meth_pos: &mut MethPos, offset: u32) {
    if let Some(qpos) = record
        .cigar()
        .read_pos(*meth_pos.position().pos() + offset, false, false)
        .unwrap()
    {
        let base = unsafe { record.seq().decoded_base_unchecked((qpos) as usize) } as char;

        let read_direction = if read_reverse_orientation(record) {
            Direction::Reverse
        } else {
            Direction::Forward
        };
        meth_pos
            .meth_bases_mut()
            .entry((read_direction, offset))
            .or_default()
            .entry(base)
            .and_modify(|count| *count += 1)
            .or_insert(1);
    }
}

pub fn write_assigned_bases(
    output: Option<PathBuf>,
    meth_positions: HashMap<PosType, NumberOfNucleotides>,
) -> Result<()> {
    let mut writer: Box<dyn Write> = match output {
        Some(path) => {
            let file = std::fs::File::create(path.clone())
                .with_context(|| format!("error creating file: {:?}", path))?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    for dir in [Direction::Forward, Direction::Reverse].iter() {
        for meth_type in [true, false].iter() {
            for cpg_pos in [0, 1] {
                let default_bases = HashMap::new();
                let bases = meth_positions
                    .get(&PosType::new(dir.clone(), *meth_type, cpg_pos))
                    .unwrap_or(&default_bases);
                let counts: Vec<usize> = ['A', 'C', 'G', 'T', 'N']
                    .iter()
                    .map(|b| *bases.get(b).unwrap_or(&0))
                    .collect();

                let orig_base = match (dir.as_str(), cpg_pos) {
                    ("forward", 0) => "C",
                    ("forward", 1) => "G",
                    ("reverse", 0) => "C",
                    ("reverse", 1) => "G",
                    _ => panic!("Unbekannte Kombination von dir und cpg_pos"), // Optional, um Fehler zu behandeln
                };
                writeln!(
                    writer,
                    "{},{},{},{}",
                    meth_type,
                    dir.as_str(),
                    orig_base,
                    counts
                        .iter()
                        .map(|c| c.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                )?;
            }
        }
    }

    Ok(())
}

/// Finds out whether the given string is a forward or reverse string.
///
/// # Returns
///
/// True if read given read is a reverse read, false if it is a forward read
pub(crate) fn read_reverse_orientation(read: &Record) -> bool {
    let read_paired = read.is_paired();
    let read_reverse = read.is_reverse();
    let read_first = read.is_first_in_template();
    if read_paired {
        read_reverse && read_first || !read_reverse && !read_first
    } else {
        read_reverse
    }
}

// Extracts all the CpG positions in a bcf to a vector
// Input:
//     bcf_file_path: Path to the bcf filewith methylation candidates
// Output:
//     Result<Vec<(String, u32)>>: Vector mit (Chromosom, Position) der jeweiligen CpG Positionen
// pub fn extract_candidate_pos(bcf_file_path: PathBuf) -> Result<Vec<MethPos>> {
//     let mut bcf = Reader::from_path(bcf_file_path)?;
//     let mut bcf_positions = Vec::new();

//     for record in bcf.records() {
//         let record = record?;
//         let chrom =
//             String::from_utf8_lossy(record.header().rid2name(record.rid().unwrap()).unwrap())
//                 .into_owned();
//         let pos = record.pos() as u32 + 1;
//         bcf_positions.push(MethPos::new(Position::new(chrom, pos, None), None));
//     }
//     println!("BCF Positions: {:?}", bcf_positions);
//     Ok(bcf_positions)
// }

// // TODO vernuenftige csv schreiben
// pub fn write_pos_to_bases(
//     output: Option<PathBuf>,
//     position_counts: HashMap<Position, HashMap<char, u32>>,
// ) -> Result<()> {
//     let output_file =
//         File::create(output.unwrap()).with_context(|| format!("error opening BCF writer"))?;
//     let mut writer = BufWriter::new(output_file);

//     // Schreiben Sie die Ergebnisse in die Datei
//     writeln!(&mut writer, "#CHROM	#POS	#DIR	#A	#C	#G	#T	#N")?;
//     for (position, counts) in &position_counts {
//         write!(
//             &mut writer,
//             "{}	{}	{}	",
//             position.chrom(),
//             position.pos(),
//             position.direction().as_ref().unwrap()
//         )?;
//         write!(
//             &mut writer,
//             "{}	{}	{}	{}	{}",
//             counts.get(&'A').unwrap_or(&0),
//             counts.get(&'C').unwrap_or(&0),
//             counts.get(&'G').unwrap_or(&0),
//             counts.get(&'T').unwrap_or(&0),
//             counts.get(&'N').unwrap_or(&0)
//         )?;
//         writeln!(&mut writer)?;
//     }
//     Ok(())
// }

//  Computes if a given read is valid (Right now we only accept specific flags)

//  # Returns

//  bool: True, if read is valid, else false
// pub fn read_invalid(flag: u16) -> bool {
//     if flag == 0 || flag == 16 || flag == 99 || flag == 83 || flag == 147 || flag == 163 {
//         return false;
//     }
//     true
// }
