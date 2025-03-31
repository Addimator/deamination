use anyhow::{Context, Ok, Result};
use deamination::utils::{Direction, MethPos, PosType, Position};
use std::collections::HashMap;

use std::{fs, path::Path, path::PathBuf};

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn get_true_counts(direction: &str) -> HashMap<PosType, HashMap<char, usize>> {
    let mut true_bases: HashMap<PosType, HashMap<char, usize>> = HashMap::new();

    if direction.contains("f0") {
        let mut value = HashMap::new();
        value.insert('T', 1);
        value.insert('C', 1);
        true_bases.insert(PosType::new(Direction::Forward, true, 0), value);
    }
    if direction.contains("f1") {
        let mut value = HashMap::new();
        value.insert('G', 2);
        true_bases.insert(PosType::new(Direction::Forward, true, 1), value);
    }
    if direction.contains("r1") {
        let mut value = HashMap::new();
        value.insert('G', 1);
        value.insert('A', 1);
        true_bases.insert(PosType::new(Direction::Reverse, true, 1), value);
    }
    if direction.contains("r0") {
        let mut value = HashMap::new();
        value.insert('C', 2);
        true_bases.insert(PosType::new(Direction::Reverse, true, 0), value);
    }
    if direction == "complete" {
        true_bases.insert(
            PosType::new(Direction::Reverse, true, 0),
            HashMap::from([('C', 1)]),
        );
        true_bases.insert(
            PosType::new(Direction::Reverse, true, 1),
            HashMap::from([('G', 1)]),
        );
        true_bases.insert(
            PosType::new(Direction::Reverse, false, 0),
            HashMap::from([('C', 3)]),
        );
        true_bases.insert(
            PosType::new(Direction::Reverse, false, 1),
            HashMap::from([('A', 1), ('G', 1), ('T', 1)]),
        );
        true_bases.insert(
            PosType::new(Direction::Forward, true, 0),
            HashMap::from([('T', 1)]),
        );
        true_bases.insert(
            PosType::new(Direction::Forward, true, 1),
            HashMap::from([('G', 1)]),
        );
        true_bases.insert(
            PosType::new(Direction::Forward, false, 0),
            HashMap::from([('C', 2), ('T', 1)]),
        );
        true_bases.insert(
            PosType::new(Direction::Forward, false, 1),
            HashMap::from([('G', 3)]),
        );
    }
    true_bases
}

fn process_bedgraph_data(
    test: &str,
    candidate_file: &str,
    true_candidates: Vec<MethPos>,
) -> Result<()> {
    let basedir = basedir(test);
    let bcf_positions = deamination::find_bases::process_bedgraph_data(PathBuf::from(format!(
        "{}/{candidate_file}.bed",
        basedir
    )))
    .with_context(|| "error computing candidate positions")?;
    assert_eq!(bcf_positions, true_candidates);
    Ok(())
}

fn count_bases_in_reads(
    test: &str,
    alignment_file: &str,
    true_candidates: Vec<MethPos>,
    true_position_counts: HashMap<PosType, HashMap<char, usize>>,
) -> Result<()> {
    let basedir = basedir(test);
    let base_counter = deamination::find_bases::count_bases_in_reads(
        PathBuf::from(format!("{}/{}.bam", basedir, alignment_file)),
        true_candidates,
    )
    .with_context(|| "error computing the position counts")?;
    assert_eq!(base_counter, true_position_counts);
    Ok(())
}

#[test]
fn test_process_bedgraph_data_simple() -> Result<()> {
    let true_candidates = vec![MethPos::new(Position::new("21".to_string(), 7), 35)];
    process_bedgraph_data("short_tests", "candidates", true_candidates)?;
    Ok(())
}

#[test]
fn test_process_bedgraph_data_complex() -> Result<()> {
    let true_candidates = vec![
        MethPos::new(Position::new("1".to_string(), 7123), 0),
        MethPos::new(Position::new("9".to_string(), 532), 100),
        MethPos::new(Position::new("21".to_string(), 1), 14),
        MethPos::new(Position::new("21".to_string(), 325), 35),
    ];
    process_bedgraph_data("short_tests", "candidates_complex", true_candidates)?;
    Ok(())
}

#[test]
fn test_count_bases_mid() -> Result<()> {
    let true_candidates = vec![MethPos::new(Position::new("21".to_string(), 7), 35)];
    let true_bases = get_true_counts("f0f1r0r1");

    count_bases_in_reads(
        "short_tests",
        "alignment_mid_c",
        true_candidates,
        true_bases,
    )?;
    Ok(())
}

#[test]
fn test_count_bases_starting_c() -> Result<()> {
    let true_candidates = vec![MethPos::new(Position::new("21".to_string(), 7), 35)];
    let true_bases = get_true_counts("f0f1r0r1");

    count_bases_in_reads(
        "short_tests",
        "alignment_start_c",
        true_candidates,
        true_bases,
    )?;
    Ok(())
}

#[test]
fn test_count_bases_starting_g() -> Result<()> {
    let true_candidates = vec![MethPos::new(Position::new("21".to_string(), 7), 35)];
    let true_bases = get_true_counts("r1f1");

    count_bases_in_reads(
        "short_tests",
        "alignment_start_g",
        true_candidates,
        true_bases,
    )?;
    Ok(())
}

#[test]
fn test_count_bases_ending_c() -> Result<()> {
    let true_candidates = vec![MethPos::new(Position::new("21".to_string(), 7), 35)];
    let true_bases = get_true_counts("f0r0");

    count_bases_in_reads(
        "short_tests",
        "alignment_end_c",
        true_candidates,
        true_bases,
    )?;
    Ok(())
}

#[test]
fn test_count_bases_ending_g() -> Result<()> {
    let true_bases = get_true_counts("f0f1r0r1");
    let true_candidates = vec![MethPos::new(Position::new("21".to_string(), 7), 35)];

    count_bases_in_reads(
        "short_tests",
        "alignment_end_g",
        true_candidates,
        true_bases,
    )?;
    Ok(())
}

#[test]
fn test_complete_workflow() -> Result<()> {
    let true_candidates = vec![
        MethPos::new(Position::new("20".to_string(), 4), 90),
        MethPos::new(Position::new("20".to_string(), 8), 18),
        MethPos::new(Position::new("21".to_string(), 3), 47),
        MethPos::new(Position::new("21".to_string(), 9), 3),
    ];
    process_bedgraph_data("test_complete", "bed_avg", true_candidates.clone())?;

    let true_bases = get_true_counts("complete");

    count_bases_in_reads("test_complete", "alignment", true_candidates, true_bases)?;
    Ok(())
}
