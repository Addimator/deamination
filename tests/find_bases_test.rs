use anyhow::{Result, Context, Ok};
use std::{fs, path::Path, path::PathBuf};
use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{Write, Read, BufReader};

use deamination;


fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}



fn get_true_counts_complete(direction: &str) -> BTreeMap<(String, u32, char), HashMap<char, u32>> {
    let mut true_bases: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();
    
    if direction != "r" {
        let mut value = HashMap::new();
        value.insert('T', 1);
        value.insert('C', 1);
        true_bases.insert(("21".to_string(), 7, 'f'), value);
    }
    if direction != "f" {
        let mut value = HashMap::new();
        value.insert('G', 1);
        value.insert('A', 1);
        true_bases.insert(("21".to_string(), 7, 'r'), value);
    }
    true_bases
}


fn extract_vcf_positions(test: &str) -> Result<()> {
    let true_candidates = vec![("21".to_string(), 7)];    
    let basedir = basedir(test);
    // let output = format!("{}/candidates.bcf", basedir);
    // cleanup_file(&output);
    let vcf_positions = deamination::find_bases::extract_vcf_positions(
        PathBuf::from(format!("{}/candidates.vcf", basedir))
    )
    .with_context(|| format!("error computing candidate positions"))?;
    assert_eq!(
        vcf_positions,
        true_candidates
    );
    Ok(())
}

fn count_bases_in_reads(test: &str, alignment_file: &str, true_position_counts: BTreeMap<(String, u32, char), HashMap<char, u32>>) -> Result<()> {
    let true_candidates = vec![("21".to_string(), 7)];
    let basedir = basedir(test);
    // let output = format!("{}/candidates.bcf", basedir);
    // cleanup_file(&output);
    let position_counts = deamination::find_bases::count_bases_in_reads(
        PathBuf::from(format!("{}/{}.sam", basedir, alignment_file)), &true_candidates
    )
    .with_context(|| format!("error computing the position counts"))?;
    assert_eq!(
        position_counts,
        true_position_counts
    );
    Ok(())
}

fn write_output(test: &str, position_counts: BTreeMap<(String, u32, char), HashMap<char, u32>>) -> Result<()> {
    let basedir = basedir(test);
    let output = format!("{}/output_test.txt", basedir);
    cleanup_file(&output);
    let true_output_content = 

"#CHROM	#POS	#DIR	#A	#C	#G	#T	#N
21	7	f	0	1	0	1	0
21	7	r	1	0	1	0	0
";

    deamination::find_bases::write_output(
        Some(PathBuf::from(output.clone())), position_counts
    )
    .with_context(|| format!("error computing the position counts"))?;
    
    let file = File::open(PathBuf::from(output.clone()))?;
    let mut buf_reader = BufReader::new(file);
    let mut output_content = String::new();
    buf_reader.read_to_string(&mut output_content)?;
    assert_eq!(output_content, true_output_content);

    Ok(())
}



#[test]
fn test_extract_vcf_positions() -> Result<()> {
    extract_vcf_positions("find_bases")?;
    Ok(())
}


#[test]
fn test_count_bases_mid() -> Result<()> {
    let true_bases = get_true_counts_complete("fr");
    count_bases_in_reads("find_bases", "alignment_mid_c", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_starting_c() -> Result<()> {
    let true_bases = get_true_counts_complete("fr");
    count_bases_in_reads("find_bases", "alignment_start_c", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_starting_g() -> Result<()> {
    let true_bases = get_true_counts_complete("r");
    count_bases_in_reads("find_bases", "alignment_start_g", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_ending_c() -> Result<()> {
    let true_bases = get_true_counts_complete("f");
    count_bases_in_reads("find_bases", "alignment_end_c", true_bases)?;
    Ok(())
}

#[test]
fn test_count_bases_ending_g() -> Result<()> {
    let true_bases = get_true_counts_complete("fg");

    count_bases_in_reads("find_bases", "alignment_end_g", true_bases)?;
    Ok(())
}

#[test]
fn test_write_output() -> Result<()> {
    let true_bases = get_true_counts_complete("fg");

    write_output("find_bases", true_bases)?;
    Ok(())
}

