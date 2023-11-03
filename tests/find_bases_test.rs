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
    let true_candidates = [("chr21".to_string(), 3)];
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

fn count_bases_in_reads(test: &str) -> Result<()> {
    let true_candidates = vec![("chr21".to_string(), 3)];
    let mut true_position_counts: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();

    let basedir = basedir(test);
    // let output = format!("{}/candidates.bcf", basedir);
    // cleanup_file(&output);
    let position_counts = deamination::find_bases::count_bases_in_reads(
        PathBuf::from(format!("{}/alignment.sam", basedir)), &true_candidates
    )
    .with_context(|| format!("error computing the position counts"))?;
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
fn test_count_bases() -> Result<()> {
    count_bases_in_reads("find_bases")?;
    Ok(())
}