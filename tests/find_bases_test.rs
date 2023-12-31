use anyhow::{Result, Context, Ok};
use std::{fs, path::Path, path::PathBuf};
use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{Read, BufReader};

use deamination;


fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}



fn get_true_counts_complete(direction: &str) -> BTreeMap<(String, u32, &str), HashMap<String, u32>> {
    let mut true_bases: BTreeMap<(String, u32, &str), HashMap<String, u32>> = BTreeMap::new();
    
    if direction.contains("f0") {
        let mut value = HashMap::new();
        value.insert("T".to_string(), 1);
        value.insert("C".to_string(), 1);
        true_bases.insert(("21".to_string(), 7, "f_0"), value);
    }
    if direction.contains("f1") {

        let mut value = HashMap::new();
        value.insert('G'.to_string(), 2);
        true_bases.insert(("21".to_string(), 7, "f_1"), value);
    }
    if direction.contains("r1") {
        let mut value = HashMap::new();
        value.insert('G'.to_string(), 1);
        value.insert('A'.to_string(), 1);
        true_bases.insert(("21".to_string(), 7, "r_1"), value);
    }
    if direction.contains("r0") {
        let mut value = HashMap::new();
        value.insert('C'.to_string(), 2);
        true_bases.insert(("21".to_string(), 7, "r_0"), value);
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

fn count_bases_in_reads(test: &str, alignment_file: &str, true_position_counts: BTreeMap<(String, u32, &str), HashMap<String, u32>>) -> Result<()> {
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


fn read_invalid() -> Result<()> {
    let invalid = deamination::find_bases::read_invalid(81);
    assert_eq!(
        invalid,
        true
    );
    Ok(())
}

fn write_pos_to_bases(test: &str, position_counts: BTreeMap<(String, u32, &str), HashMap<String, u32>>) -> Result<()> {
    let basedir = basedir(test);

    // Read true output file and get its content
    let true_pos_to_bases = format!("{}/pos_to_bases.txt", basedir);
    let file = File::open(PathBuf::from(true_pos_to_bases.clone()))?;
    let mut buf_reader = BufReader::new(file);
    let mut true_output_content = String::new();
    buf_reader.read_to_string(&mut true_output_content)?;

    // Compute output of method
    let output = format!("{}/output_test.txt", basedir);
    cleanup_file(&output);
    deamination::find_bases::write_pos_to_bases(
        Some(PathBuf::from(output.clone())), position_counts
    )
    .with_context(|| format!("error computing the position counts"))?;
    
    // Read method output file and get its content
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
    let true_bases = get_true_counts_complete("f0f1r0r1");
    count_bases_in_reads("find_bases", "alignment_mid_c", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_starting_c() -> Result<()> {
    let true_bases = get_true_counts_complete("f0f1r0r1");
    count_bases_in_reads("find_bases", "alignment_start_c", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_starting_g() -> Result<()> {
    let true_bases = get_true_counts_complete("r1f1");
    count_bases_in_reads("find_bases", "alignment_start_g", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_ending_c() -> Result<()> {
    let true_bases = get_true_counts_complete("f0r0");
    count_bases_in_reads("find_bases", "alignment_end_c", true_bases)?;
    Ok(())
}

#[test]
fn test_count_bases_ending_g() -> Result<()> {
    let true_bases = get_true_counts_complete("f0f1r0r1");

    count_bases_in_reads("find_bases", "alignment_end_g", true_bases)?;
    Ok(())
}

#[test]
fn test_read_invalid() -> Result<()> {
    read_invalid()?;
    Ok(())
}

#[test]
fn test_write_pos_to_bases() -> Result<()> {
    let true_bases = get_true_counts_complete("f0f1r0r1");

    write_pos_to_bases("find_bases", true_bases)?;
    Ok(())
}

