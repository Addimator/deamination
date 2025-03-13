use crate::utils::{Direction, MethPos, NumberOfNucleotides, PosType, Position};
use anyhow::{Context, Result};

use rust_htslib::bam::Record;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::rc::Rc;

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

        let methylation = bed_fields[3].parse().context("Invalid methylation value")?;
        let position = MethPos::new(Position::new(chrom.clone(), pos as u32), Some(methylation));

        meth_positions.push(position);
    }
    Ok(meth_positions)
}

pub fn get_relevant_reads(bcf_positions: &[MethPos]) -> (String, u64, u64) {
    let chrom = bcf_positions[0].position().chrom().clone(); // Chromosom vom ersten Eintrag
    let start = (bcf_positions[0].position().pos() - 1).saturating_sub(201) as u64; // Startposition vom ersten Element
    let end = bcf_positions
        .last()
        .ok_or_else(|| anyhow::anyhow!("bcf hat keine Positionen"))
        .map(|x| x.position().pos() + 1)
        .unwrap()
        .saturating_add(200) as u64; // Endposition vom letzten Element

    // let header = bam_reader.tid(); // Kopiere den Header einmal
    // let tid: u32 = header.tid(chrom.as_bytes()).unwrap();
    // first pass, extend reference interval
    // TODO: Chrom oder TID zurueckgeben?
    (chrom, start, end)
}

/// Computes the number of different bases on the CpG positions. Different entries for reads on the forward and reverse strand.
/// Input:
///     sam_file_path: Path to the alignment
///     bcf_positions: Vector with all the CpG positions
/// Output:
///      Result<BTreeMap<(String, u32, char), HashMap<char, u32>>>: Maps a specific position (chromosome, position, read_direction) to a map of bases to their number on this position
// TODO: Test ueber unterschiedliche Chromosome hinweg
// TODO: Hier schon MethPositions verwenden und nicht HashMaps mit Positions und Bases
pub fn count_bases_in_reads(
    bam_file_path: PathBuf,
    mut bcf_positions: Vec<MethPos>,
) -> Result<HashMap<PosType, NumberOfNucleotides>> {
    let mut indexed_reader =
        bam::IndexedReader::from_path(bam_file_path).context("Unable to read BAM/CRAM file.")?;
    // let mut bam_reader = bam::RecordBuffer::new(indexed_bam, true);

    let mut position_counts = HashMap::new();

    for bcf_position in bcf_positions.iter_mut() {
        let chrom = bcf_position.position().chrom();
        let start = *bcf_position.position().pos() as i64;
        let end = start + 1;

        indexed_reader.fetch((
            indexed_reader.header().tid(chrom.as_bytes()).unwrap(),
            start,
            end,
        ))?;

        for record in indexed_reader.records() {
            let record = record?;
            insert_bases(&record, bcf_position, 0);
            insert_bases(&record, bcf_position, 1);
        }
        // Fill 8 types of methylation position hashmaps
        for ((direction, cpg_position), value) in bcf_position.meth_bases().iter() {
            let pos_type = PosType::new(
                direction.clone(),
                bcf_position.methylation().unwrap() > 20.0,
                *cpg_position,
            );
            position_counts
                .entry(pos_type)
                .and_modify(|existing_meth_pos: &mut HashMap<char, usize>| {
                    for (base, count) in value {
                        // println!("Base: {}, Count: {}", base, count);
                        *existing_meth_pos.entry(*base).or_insert(0) += count;
                    }
                })
                .or_insert_with(HashMap::new);
        }
        println!("Position: {:?}", bcf_position);
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
        // dbg!(&read_direction, &offset, &base);
        meth_pos
            .meth_bases_mut()
            .entry((read_direction, offset))
            .or_default() // Falls nicht vorhanden, lege eine neue HashMap an
            .entry(base)
            .and_modify(|count| *count += 1) // Falls vorhanden, erhöhe um 1
            .or_insert(1);
        // dbg!(&meth_pos);
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
    writeln!(writer, "direction,methylation_status,cpg_position,bases")?;
    for dir in [Direction::Forward, Direction::Reverse].iter() {
        for meth_type in [true, false].iter() {
            for cpg_pos in [0, 1] {
                let bases = meth_positions.get(&PosType::new(dir.clone(), *meth_type, cpg_pos));
                writeln!(
                    writer,
                    "{},{},{},{:?}",
                    dir.as_str(),
                    meth_type,
                    cpg_pos,
                    bases
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
// pub fn extract_bcf_positions(bcf_file_path: PathBuf) -> Result<Vec<MethPos>> {
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
