use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, BufRead, Write, BufWriter};
use std::path::PathBuf;
use anyhow::{Context, Result};


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
pub fn count_bases_in_reads(sam_file_path: PathBuf, vcf_positions: &Vec<(String, u32)>) -> Result<BTreeMap<(String, u32, char), HashMap<char, u32>>> {
    let sam_file = File::open(sam_file_path)?;
    let sam_reader = BufReader::new(sam_file);

    // Initialize a HashMap to store the counts
    let mut position_counts: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();
    let mut bed_graph_entries: HashMap<(String, u32, u32, u16, i32), (char, String)> = HashMap::new();
    let mut id = 0; // We need the ID because different aligned reads can look the same
    
    let mut read_to_dir: HashMap<u16, char> = HashMap::new();
    // Process the SAM and fill bed_graph_entries with the lines
    for line in sam_reader.lines() {
        let line = line.unwrap();
        if line.starts_with('@') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue; // Malformed SAM line, skip
        }
        let flag = fields[1].parse::<u16>().unwrap();
        let chrom = fields[2].to_string();
        let position = fields[3].parse::<u32>().unwrap();
        let sequence = fields[9].to_string();
        if read_invalid(flag) {
            read_to_dir.insert(flag, 'i');
            continue;
        }
        let reverse_read = read_reverse_strand(flag);
        let direction = if reverse_read { 'r' } else { 'f' };
        read_to_dir.insert(flag, direction);
        bed_graph_entries.insert((chrom.clone(), position, position + sequence.len() as u32, flag, id), (direction, sequence));
        id += 1;
    }
    // Go through all CpG positions and find the bases from the bedGraph entries
    for (chrom, vcf_pos) in vcf_positions.iter() {
        for ((bed_chrom, start_pos, end_pos, _flag, _id), (read_dir, sequence)) in &bed_graph_entries {
            let base_pos;
            if *read_dir == 'f' {
                if chrom == bed_chrom && vcf_pos >= start_pos && vcf_pos < end_pos {
                    base_pos = *vcf_pos;
                    
                    let entry = position_counts
                    .entry((chrom.clone(), *vcf_pos, *read_dir))
                    .or_insert(HashMap::new());
                    *entry.entry(sequence.as_bytes()[(base_pos - start_pos) as usize] as char).or_insert(0) += 1;
                }
            }
            else {
                if chrom == bed_chrom && (vcf_pos + 1) >= *start_pos && (vcf_pos + 1) < *end_pos {
                    base_pos = *vcf_pos + 1;
                    
                    let entry = position_counts
                    .entry((chrom.clone(), *vcf_pos, *read_dir))
                    .or_insert(HashMap::new());
                    *entry.entry(sequence.as_bytes()[(base_pos - start_pos) as usize] as char).or_insert(0) += 1;
                }

            }
        }
    }
    println!("Reads to dir: {:?}", read_to_dir);
    Ok(position_counts)
}

pub fn write_pos_to_bases(output: Option<PathBuf>, position_counts: BTreeMap<(String, u32, char), HashMap<char, u32>>) -> Result<()>{
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

pub fn read_reverse_strand(flag:u16) -> bool {
    let read_paired = 0b1;
    let read_mapped_porper_pair = 0b01;
    let read_reverse = 0b10000;
    let mate_reverse = 0b100000;
    let first_in_pair = 0b1000000;
    let second_in_pair = 0b10000000;
    if (flag & read_paired) != 0 && (flag & read_mapped_porper_pair) != 0 {
        if (flag & read_reverse) != 0 && (flag & first_in_pair) != 0 {
            return true
        }
        else if (flag & mate_reverse) != 0 && (flag & second_in_pair) != 0 {
            return true
        }
    }
    else {
        if (flag & read_reverse) != 0 {
            return true
        }
    }
    false
    // read.inner.core.flag == 163 || read.inner.core.flag == 83 || read.inner.core.flag == 16
}

pub fn read_invalid(flag:u16) -> bool {
    // let read_paired = 0b1;
    // let read_mapped_porper_pair = 0b01;
    // let read_unmapped = 0b100;
    // let mate_unmapped = 0b1000;
    // let read_reverse_strand = 0b10000;
    // let mate_reverse_strand = 0b100000;
    
    // let secondary_alignment = 0b100000000;
    // let qc_failed = 0b1000000000;
    // let duplicate = 0b10000000000;
    // let supplemental = 0b100000000000;
    // if (flag & secondary_alignment) != 0 
    // || (flag & qc_failed) != 0 
    // || (flag & duplicate) != 0 
    // || (flag & supplemental) != 0
    // || (flag & read_unmapped) != 0
    // || (flag & mate_unmapped)!= 0
    // {
    //     return true
    // }
    // // invalid if both pairs are reverse or forward
    // if (flag & read_paired) != 0 && (flag & read_mapped_porper_pair) != 0 {
    //     if  (flag & read_reverse_strand) != 0 && (flag & mate_reverse_strand) != 0
    //     || (flag & read_reverse_strand) == 0 && (flag & mate_reverse_strand) == 0
    //     {
    //         return true
    //     }
    // }
    // false    
    if flag == 0 || flag == 16 || flag == 99 || flag == 83 || flag == 147 || flag == 163 {
        return false
    }
    return true

}