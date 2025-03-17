// use anyhow::{Context, Ok, Result};
// use deamination::utils::Position;
// use std::collections::{BTreeMap, HashMap};
// use std::fs::File;
// use std::hash::Hash;
// use std::io::{BufReader, Read};
// use std::{fs, path::Path, path::PathBuf};

// fn basedir(test: &str) -> String {
//     format!("tests/resources/{}", test)
// }

// fn cleanup_file(f: &str) {
//     if Path::new(f).exists() {
//         fs::remove_file(f).unwrap();
//     }
// }

// fn get_true_counts_complete(direction: &str) -> HashMap<(Position), HashMap<char, u32>> {
//     let mut true_bases: HashMap<(Position), HashMap<char, u32>> = HashMap::new();

//     if direction.contains("f0") {
//         let mut value = HashMap::new();
//         value.insert('T', 1);
//         value.insert('C', 1);
//         true_bases.insert(
//             Position::new("21".to_string(), 7, Some("f_0".to_string())),
//             value,
//         );
//     }
//     if direction.contains("f1") {
//         let mut value = HashMap::new();
//         value.insert('G', 2);
//         true_bases.insert(
//             Position::new("21".to_string(), 8, Some("f_1".to_string())),
//             value,
//         );
//     }
//     if direction.contains("r1") {
//         let mut value = HashMap::new();
//         value.insert('G', 1);
//         value.insert('A', 1);
//         true_bases.insert(
//             Position::new("21".to_string(), 8, Some("r_1".to_string())),
//             value,
//         );
//     }
//     if direction.contains("r0") {
//         let mut value = HashMap::new();
//         value.insert('C', 2);
//         true_bases.insert(
//             Position::new("21".to_string(), 7, Some("r_0".to_string())),
//             value,
//         );
//     }
//     true_bases
// }

// fn extract_bcf_positions(
//     test: &str,
//     candidate_file: &str,
//     true_candidates: Vec<(Position)>,
// ) -> Result<()> {
//     let basedir = basedir(test);
//     // let output = format!("{}/candidates.bcf", basedir);
//     // cleanup_file(&output);
//     let bcf_positions = deamination::find_bases::extract_bcf_positions(PathBuf::from(format!(
//         "{}/{candidate_file}.bcf",
//         basedir
//     )))
//     .with_context(|| format!("error computing candidate positions"))?;
//     assert_eq!(bcf_positions, true_candidates);
//     Ok(())
// }

// fn count_bases_in_reads(
//     test: &str,
//     alignment_file: &str,
//     true_position_counts: HashMap<(Position), HashMap<char, u32>>,
// ) -> Result<()> {
//     let true_candidates = vec![Position::new("21".to_string(), 7, None)];
//     let basedir = basedir(test);
//     // Todo: Better name?
//     let position_counts = deamination::find_bases::count_bases_in_reads(
//         PathBuf::from(format!("{}/{}.bam", basedir, alignment_file)),
//         &true_candidates,
//     )
//     .with_context(|| format!("error computing the position counts"))?;
//     assert_eq!(position_counts, true_position_counts);
//     Ok(())
// }

// // Todo: Can be removed?
// fn read_invalid() -> Result<()> {
//     let invalid = deamination::find_bases::read_invalid(81);
//     assert_eq!(invalid, true);
//     Ok(())
// }

// // fn write_pos_to_bases(
// //     test: &str,
// //     position_counts: HashMap<(Position), HashMap<char, u32>>,
// // ) -> Result<()> {
// //     let basedir = basedir(test);

// //     // Read true output file and get its content
// //     let true_pos_to_bases = format!("{}/pos_to_bases.txt", basedir);
// //     let file = File::open(PathBuf::from(true_pos_to_bases.clone()))?;
// //     let mut buf_reader = BufReader::new(file);
// //     let mut true_output_content = String::new();
// //     buf_reader.read_to_string(&mut true_output_content)?;

// //     // Compute output of method
// //     let output = format!("{}/output_test.txt", basedir);
// //     cleanup_file(&output);
// //     deamination::find_bases::write_pos_to_bases(
// //         Some(PathBuf::from(output.clone())),
// //         position_counts,
// //     )
// //     .with_context(|| format!("error computing the position counts"))?;

// //     // Read method output file and get its content
// //     let file = File::open(PathBuf::from(output.clone()))?;
// //     let mut buf_reader = BufReader::new(file);
// //     let mut output_content = String::new();
// //     buf_reader.read_to_string(&mut output_content)?;
// //     assert_eq!(output_content, true_output_content);

// //     Ok(())
// // }

// #[test]
// fn test_extract_bcf_positions() -> Result<()> {
//     let true_candidates = vec![Position::new("21".to_string(), 7, None)];
//     extract_bcf_positions("find_bases", "candidates", true_candidates)?;
//     let true_candidates_complex = vec![
//         Position::new("1".to_string(), 7123, None),
//         Position::new("9".to_string(), 532, None),
//         Position::new("21".to_string(), 1, None),
//         Position::new("21".to_string(), 325, None),
//     ];
//     extract_bcf_positions("find_bases", "candidates_complex", true_candidates_complex)?;
//     Ok(())
// }

// #[test]
// fn test_count_bases_mid() -> Result<()> {
//     let true_bases = get_true_counts_complete("f0f1r0r1");
//     count_bases_in_reads("find_bases", "alignment_mid_c", true_bases)?;
//     Ok(())
// }

// #[test]
// fn test_count_bases_starting_c() -> Result<()> {
//     let true_bases = get_true_counts_complete("f0f1r0r1");
//     count_bases_in_reads("find_bases", "alignment_start_c", true_bases)?;
//     Ok(())
// }

// #[test]
// fn test_count_bases_starting_g() -> Result<()> {
//     let true_bases = get_true_counts_complete("r1f1");
//     count_bases_in_reads("find_bases", "alignment_start_g", true_bases)?;
//     Ok(())
// }

// #[test]
// fn test_count_bases_ending_c() -> Result<()> {
//     let true_bases = get_true_counts_complete("f0r0");
//     count_bases_in_reads("find_bases", "alignment_end_c", true_bases)?;
//     Ok(())
// }

// #[test]
// fn test_count_bases_ending_g() -> Result<()> {
//     let true_bases = get_true_counts_complete("f0f1r0r1");

//     count_bases_in_reads("find_bases", "alignment_end_g", true_bases)?;
//     Ok(())
// }

// // #[test]
// // fn test_write_pos_to_bases() -> Result<()> {
// //     let true_bases = get_true_counts_complete("f0f1r0r1");

// //     write_pos_to_bases("find_bases", true_bases)?;
// //     Ok(())
// // }
