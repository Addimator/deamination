use anyhow::{Result, Context};
use std::{fs, path::Path, path::PathBuf};
use std::collections::{HashMap, BTreeMap};


use deamination;


fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
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
    println!("POSITION COUNTS: {:?}", position_counts);
    assert_eq!(
        position_counts,
        true_position_counts
    );
    Ok(())
}


#[test]
fn test_extract_vcf_positions() -> Result<()> {
    extract_vcf_positions("find_bases")?;
    Ok(())
}


#[test]
fn test_count_bases_mid() -> Result<()> {
    let mut true_bases: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();

    let mut value1 = HashMap::new();
    value1.insert('T', 1);
    value1.insert('C', 1);

    let mut value2 = HashMap::new();
    value2.insert('G', 1);
    value2.insert('A', 1);

    true_bases.insert(("21".to_string(), 7, 'f'), value1);
    true_bases.insert(("21".to_string(), 7, 'r'), value2);
    count_bases_in_reads("find_bases", "alignment_mid_c", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_starting_c() -> Result<()> {
    let mut true_bases: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();

    let mut value1 = HashMap::new();
    value1.insert('T', 1);
    value1.insert('C', 1);

    let mut value2 = HashMap::new();
    value2.insert('G', 1);
    value2.insert('A', 1);

    true_bases.insert(("21".to_string(), 7, 'f'), value1);
    true_bases.insert(("21".to_string(), 7, 'r'), value2);
    count_bases_in_reads("find_bases", "alignment_start_c", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_starting_g() -> Result<()> {
    let mut true_bases: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();


    let mut value2 = HashMap::new();
    value2.insert('G', 1);
    value2.insert('A', 1);

    true_bases.insert(("21".to_string(), 7, 'r'), value2);
    count_bases_in_reads("find_bases", "alignment_start_g", true_bases)?;
    Ok(())
}


#[test]
fn test_count_bases_ending_c() -> Result<()> {
    let mut true_bases: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();

    let mut value1 = HashMap::new();
    value1.insert('T', 1);
    value1.insert('C', 1);

    true_bases.insert(("21".to_string(), 7, 'f'), value1);
    count_bases_in_reads("find_bases", "alignment_end_c", true_bases)?;
    Ok(())
}

#[test]
fn test_count_bases_ending_g() -> Result<()> {
    let mut true_bases: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();

    let mut value1 = HashMap::new();
    value1.insert('T', 1);
    value1.insert('C', 1);

    let mut value2 = HashMap::new();
    value2.insert('G', 1);
    value2.insert('A', 1);

    true_bases.insert(("21".to_string(), 7, 'f'), value1);
    true_bases.insert(("21".to_string(), 7, 'r'), value2);
    count_bases_in_reads("find_bases", "alignment_end_g", true_bases)?;
    Ok(())
}