use anyhow::{Context, Ok, Result};
use deamination::utils::{Direction, MethPos, NumberOfNucleotides, PosType, Position};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, Read};
use std::{fs, path::Path, path::PathBuf};

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn extract_bed_positions(
    test: &str,
    candidate_file: &str,
    true_candidates: Vec<MethPos>,
) -> Result<()> {
    let basedir = basedir(test);
    // let output = format!("{}/candidates.bcf", basedir);
    // cleanup_file(&output);
    let bcf_positions = deamination::find_bases::process_bedgraph_data(PathBuf::from(format!(
        "{}/{candidate_file}.bed",
        basedir
    )))
    .with_context(|| format!("error computing candidate positions"))?;
    assert_eq!(bcf_positions, true_candidates);
    Ok(())
}

fn get_true_counts_complete(
    direction: &str,
) -> (Vec<MethPos>, HashMap<PosType, NumberOfNucleotides>) {
    let mut candidates: Vec<MethPos> = Vec::new();
    let mut base_map = HashMap::new();

    if direction.contains("f0") {
        let mut meth_bases = HashMap::new();
        meth_bases.insert('T', 1);
        meth_bases.insert('C', 1);
        let position = Position::new("21".to_string(), 7);
        let meth_pos = MethPos::new(position, 0);
        candidates.push(meth_pos);
        base_map.insert(PosType::new(Direction::Forward, false, 0), meth_bases);
    }
    if direction.contains("f1") {
        let mut meth_bases = HashMap::new();
        meth_bases.insert('G', 2);
        let position = Position::new("21".to_string(), 8);
        let meth_pos = MethPos::new(position, 70);
        candidates.push(meth_pos);
        base_map.insert(PosType::new(Direction::Forward, true, 1), meth_bases);
    }
    if direction.contains("r1") {
        let mut meth_bases = HashMap::new();
        meth_bases.insert('G', 1);
        meth_bases.insert('A', 1);
        let position = Position::new("21".to_string(), 8);
        let meth_pos = MethPos::new(position, 0);
        candidates.push(meth_pos);
        base_map.insert(PosType::new(Direction::Reverse, false, 1), meth_bases);
    }
    if direction.contains("r0") {
        let mut meth_bases = HashMap::new();
        meth_bases.insert('C', 2);
        let position = Position::new("21".to_string(), 7);
        let meth_pos = MethPos::new(position, 0);
        candidates.push(meth_pos);
        base_map.insert(PosType::new(Direction::Reverse, false, 0), meth_bases);
    }
    (candidates, base_map)
}

fn count_bases_in_reads(
    test: &str,
    alignment_file: &str,
    candidates: Vec<MethPos>,
    true_position_counts: HashMap<PosType, NumberOfNucleotides>,
) -> Result<()> {
    let basedir = basedir(test);
    // Todo: Better name?
    let position_counts = deamination::find_bases::count_bases_in_reads(
        PathBuf::from(format!("{}/{}.bam", basedir, alignment_file)),
        candidates,
    )
    .with_context(|| format!("error computing the position counts"))?;
    assert_eq!(position_counts, true_position_counts);
    Ok(())
}

#[test]
fn test_extract_candidates() -> Result<()> {
    let true_candidates = vec![
        MethPos::new(Position::new("J02458".to_string(), 2), 0),
        MethPos::new(Position::new("J02458".to_string(), 8), 100),
        MethPos::new(Position::new("J02459".to_string(), 3), 14),
        MethPos::new(Position::new("J02459".to_string(), 6), 35),
        MethPos::new(Position::new("J02459".to_string(), 12), 55),
    ];

    extract_bed_positions("find_bases", "candidates", true_candidates)?;
    Ok(())
}

#[test]
fn test_count_bases_mid() -> Result<()> {
    let (candidates, meth_pos) = get_true_counts_complete("f0f1r0r1");
    count_bases_in_reads("find_bases", "alignment_mid_c", candidates, meth_pos)?;
    Ok(())
}

#[test]
fn test_count_bases_starting_c() -> Result<()> {
    let (candidates, meth_pos) = get_true_counts_complete("f0f1r0r1");
    count_bases_in_reads("find_bases", "alignment_start_c", candidates, meth_pos)?;
    Ok(())
}

#[test]
fn test_count_bases_starting_g() -> Result<()> {
    let (candidates, meth_pos) = get_true_counts_complete("r1f1");
    count_bases_in_reads("find_bases", "alignment_start_g", candidates, meth_pos)?;
    Ok(())
}

#[test]
fn test_count_bases_ending_c() -> Result<()> {
    let (candidates, meth_pos) = get_true_counts_complete("f0r0");
    count_bases_in_reads("find_bases", "alignment_end_c", candidates, meth_pos)?;
    Ok(())
}

#[test]
fn test_count_bases_ending_g() -> Result<()> {
    let (candidates, meth_pos) = get_true_counts_complete("f0f1r0r1");

    count_bases_in_reads("find_bases", "alignment_end_g", candidates, meth_pos)?;
    Ok(())
}
