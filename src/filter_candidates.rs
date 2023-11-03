
use std::fs::File;
use std::io::{BufReader, Result, BufRead};
use std::path::PathBuf;



pub fn filter_mutations(mut vcf_positions: Vec<(String, u32)>, mutations: PathBuf) -> Result<Vec<(String, u32)>> {
    let mutations_file = File::open(mutations)?;
    let vcf_reader = BufReader::new(mutations_file);
    let mut to_remove = Vec::new();

    for line in vcf_reader.lines() {
        let line = line?;
        if !line.starts_with('#') {
            // Skip comment lines
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                if let [chrom, pos, _, reference, alternate, ..] = &fields[..] {
                    let chrom = chrom.trim_start_matches("chr");

                    if let Ok(pos) = pos.parse::<u32>() {
                        if vcf_positions.contains(&(chrom.to_string(), pos)) {
                            to_remove.push((chrom.to_string(), pos));
                            println!("VCF filtering: Chrom {:?} Pos {:?}", chrom, pos);
                        }
                    }
                }
            }
        }
    }
    // Use the `retain` method to keep only the elements not in `to_remove`
    vcf_positions.retain(|(chrom, pos)| !to_remove.contains(&(chrom.to_owned(), *pos)));

    Ok(vcf_positions)
}


pub fn filter_inconsistent(mut vcf_positions: Vec<(String, u32)>, consistent: PathBuf) -> Result<Vec<(String, u32)>> {
    let bedgraph_file = File::open(consistent).expect("Unable to open bedGraph file");
    let bedgraph_reader = BufReader::new(bedgraph_file);
    let mut vcf_positions_filtered = Vec::new();

    for bed_line in bedgraph_reader.lines() {
        let bed_line = bed_line.expect("Error reading bedGraph line");
        if bed_line.starts_with("track") {
            continue;
        }
        let bed_fields: Vec<&str> = bed_line.split('\t').collect();
    
        if let [chrom_orig, start, end, ..] = &bed_fields[..] {
            let chrom = chrom_orig.trim_start_matches("chr");
            let start = start.parse::<u32>().expect("Invalid position value") + 1;
            let end = end.parse::<u32>().expect("Invalid position value");

            for pos in start..end {
                if vcf_positions.contains(&(chrom.to_string(), pos)) {
                    vcf_positions_filtered.push((chrom.to_string(), pos));
                    println!("Bedgraph filtering: Chrom {:?} Pos {:?}", chrom, pos);
                }
            }
        }
    }

    Ok(vcf_positions_filtered)
        
}