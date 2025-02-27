use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{BufReader, BufRead, Write, BufWriter};
use std::path::PathBuf;
use anyhow::{Context, Result};
use regex::Regex;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::{bam, bcf, bcf::Read};

fn is_only_matches(cigar: &str) -> bool {
    let re = Regex::new(r"^(\d+M)+$").unwrap();
    re.is_match(cigar)
}


/// Extracts all the CpG positions in a vcf to a vector
/// Input: 
///     vcf_file_path: Path to the vcf filewith methylation candidates
/// Output:
///     Result<Vec<(String, u32)>>: Vector mit (Chromosom, Position) der jeweiligen CpG Positionen
pub fn extract_vcf_positions(vcf_file_path: PathBuf) -> Result<Vec<(String, u32)>> {
    let vcf_file = File::open(vcf_file_path)?;
    let vcf_reader = BufReader::new(vcf_file);

    let mut vcf_positions = Vec::new();

    for line in vcf_reader.lines() {
        let line = line?;
        // println!("{:?}", line);
        if !line.starts_with('#') {
            // Skip comment lines

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                if let Ok(chrom) = fields[0].parse::<String>() {
                    if let Ok(position) = fields[1].parse::<u32>() {
                        vcf_positions.push((chrom, position));
                    }
                }
            }
        }
    }

    Ok(vcf_positions)
}

/// Computes the number of different bases on the CpG positions. Different entries for reads on the forward and reverse strand.
/// Input:
///     sam_file_path: Path to the alignment
///     vcf_positions: Vector with all the CpG positions
/// Output:
///      Result<BTreeMap<(String, u32, char), HashMap<char, u32>>>: Maps a specific position (chromosome, position, read_direction) to a map of bases to their number on this position
pub fn count_bases_in_reads(bam_file_path: PathBuf, vcf_positions: &Vec<(String, u32)>) -> Result<BTreeMap<(String, u32, &str), HashMap<char, u32>>> {
    let mut bam_reader = bam::IndexedReader::from_path(bam_file_path)?;

    // let mut bam_reader = Reader::from_path(bam_file_path)?;
    let mut position_counts: BTreeMap<(String, u32, &str), HashMap<char, u32>> = BTreeMap::new();
    // Initialize a HashMap to store the counts
    let mut bed_graph_entries: HashMap<(String, u32, u32, u16, i32), (char, Vec<u8>)> = HashMap::new();
    let mut id = 0; // We need the ID because different aligned reads can look the same

    let chrom = vcf_positions[0].0.clone(); // Chromosom vom ersten Eintrag
    let start = vcf_positions[0].1 - 1; // Startposition vom ersten Element
    let end = vcf_positions
        .last()
        .ok_or_else(|| anyhow::anyhow!("VCF hat keine Positionen")).map(|x| x.1 + 1)?; // Endposition vom letzten Element

    let header = bam_reader.header().to_owned(); // Kopiere den Header einmal
    let tid: u32 = header.tid(chrom.as_bytes()).unwrap();
    // first pass, extend reference interval
    bam_reader.fetch((tid, start, end))?;
    for result in bam_reader.records() {
        let record = result?;
        let chrom = std::str::from_utf8(header.tid2name(record.tid() as u32))
            .unwrap()
            .to_string();
        let position = record.pos() as u32 + 1; // BAM ist 0-basiert, VCF ist 1-basiert
        let cigar = record.cigar().to_string();
        let sequence = record.seq().as_bytes();
        let flag = record.inner.core.flag;
        if read_invalid(flag) || !is_only_matches(&cigar) {
            continue;
        }

        let reverse_read = read_reverse_strand(flag);
        let direction = if reverse_read { 'r' } else { 'f' };
        bed_graph_entries.insert((chrom.clone(), position, position + sequence.len() as u32, record.inner.core.flag, id), (direction, sequence));
        id += 1;
    }
    // Go through all CpG positions and find the bases from the bedGraph entries
    for (chrom, vcf_pos) in vcf_positions.iter() {
        for ((bed_chrom, start_pos, end_pos, _flag, _id), (read_dir, sequence)) in &bed_graph_entries {
            let mut base_pos;
            if *read_dir == 'f' {
                if chrom == bed_chrom && vcf_pos >= start_pos && vcf_pos < end_pos {
                    base_pos = *vcf_pos;
                    let base = char::from(sequence[(base_pos - start_pos) as usize]);
                    let entry = position_counts
                        .entry((chrom.clone(), *vcf_pos, "f_0"))
                        .or_insert(HashMap::new());
                        *entry.entry(base).or_insert(0) += 1;
                }
                if chrom == bed_chrom && (vcf_pos + 1) >= *start_pos && (vcf_pos + 1) < *end_pos {
                    base_pos = *vcf_pos + 1;
                    let base = char::from(sequence[(base_pos - start_pos) as usize]);
                    let entry = position_counts
                        .entry((chrom.clone(), *vcf_pos, "f_1"))
                        .or_insert(HashMap::new());
                        *entry.entry(base).or_insert(0) += 1;
                }
            }
            else {
                if chrom == bed_chrom && vcf_pos  >= start_pos && vcf_pos < end_pos {
                    base_pos = *vcf_pos;
                    let base = char::from(sequence[(base_pos - start_pos) as usize]);
                    let entry = position_counts
                        .entry((chrom.clone(), *vcf_pos, "r_0"))
                        .or_insert(HashMap::new());
                        *entry.entry(base).or_insert(0) += 1;
                }
                if chrom == bed_chrom && (vcf_pos + 1) >= *start_pos && (vcf_pos + 1) < *end_pos {
                    base_pos = *vcf_pos + 1;
                    let base = char::from(sequence[(base_pos - start_pos) as usize]);
                    let entry = position_counts
                        .entry((chrom.clone(), *vcf_pos, "r_1"))
                        .or_insert(HashMap::new());
                        *entry.entry(base).or_insert(0) += 1;
                }

            }
        }
    }
    println!("Reads to dir: {:?}", position_counts);
    Ok(position_counts)
}

pub fn write_pos_to_bases(output: Option<PathBuf>, position_counts: BTreeMap<(String, u32, &str), HashMap<char, u32>>) -> Result<()>{
    let output_file = File::create(output.unwrap()).with_context(|| format!("error opening BCF writer"))?;
    let mut writer = BufWriter::new(output_file);

    // Schreiben Sie die Ergebnisse in die Datei
    writeln!(&mut writer, "#CHROM	#POS	#DIR	#A	#C	#G	#T	#N")?;
    for ((reference_name, position, direction), counts) in &position_counts {
        write!(
            &mut writer,
            "{}	{}	{}	",
            reference_name, position, direction
        )?;
        write!(&mut writer, "{}	{}	{}	{}	{}", counts.get(&'A').unwrap_or(&0), counts.get(&'C').unwrap_or(&0), counts.get(&'G').unwrap_or(&0), counts.get(&'T').unwrap_or(&0), counts.get(&'N').unwrap_or(&0))?;
        writeln!(&mut writer)?;
    }
    Ok(())
}

/// Finds out whether the given string is a forward or reverse string.
///
/// # Returns
///
/// True if read given read is a reverse read, false if it is a forward read
pub(crate) fn read_reverse_strand(flag: u16) -> bool {
    let read_paired = 0b1;
    let read_mapped_porper_pair = 0b01;
    let read_reverse = 0b10000;
    let mate_reverse = 0b100000;
    let first_in_pair = 0b1000000;
    let second_in_pair = 0b10000000;
    if (flag & read_paired != 0
        && flag & read_mapped_porper_pair != 0
        && flag & read_reverse != 0
        && flag & first_in_pair != 0)
        || (flag & read_paired != 0
            && flag & read_mapped_porper_pair != 0
            && flag & mate_reverse != 0
            && flag & second_in_pair != 0)
        || (flag & read_reverse != 0
            && (flag & read_paired == 0 || flag & read_mapped_porper_pair == 0))
    {
        return true;
    }
    false
    // flag & 0x10 != 0
}

/// Computes if a given read is valid (Right now we only accept specific flags)
///
/// # Returns
///
/// bool: True, if read is valid, else false
fn read_invalid(flag: u16) -> bool {
    if flag == 0 || flag == 16 || flag == 99 || flag == 83 || flag == 147 || flag == 163 {
        return false;
    }
    true
}