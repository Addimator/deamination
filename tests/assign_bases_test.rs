// use anyhow::{Context, Ok, Result};
// use std::collections::{BTreeMap, HashMap};
// use std::fs::File;
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

// fn create_true_baseline_maps() -> Result<(
//     HashMap<(String, usize), Vec<String>>,
//     HashMap<(String, usize), Vec<String>>,
//     HashMap<(String, usize), Vec<String>>,
//     HashMap<(String, usize), Vec<String>>,
// )> {
//     let mut forward_baseline_map_0 = HashMap::new();
//     let mut forward_baseline_map_1 = HashMap::new();
//     let mut reverse_baseline_map_0 = HashMap::new();
//     let mut reverse_baseline_map_1 = HashMap::new();
//     forward_baseline_map_0.insert(
//         ("21".to_string(), 7),
//         vec![
//             "0".to_string(),
//             "1".to_string(),
//             "0".to_string(),
//             "1".to_string(),
//             "0".to_string(),
//         ],
//     );
//     forward_baseline_map_1.insert(
//         ("21".to_string(), 7),
//         vec![
//             "0".to_string(),
//             "0".to_string(),
//             "2".to_string(),
//             "0".to_string(),
//             "0".to_string(),
//         ],
//     );
//     reverse_baseline_map_0.insert(
//         ("21".to_string(), 7),
//         vec![
//             "0".to_string(),
//             "2".to_string(),
//             "0".to_string(),
//             "0".to_string(),
//             "0".to_string(),
//         ],
//     );
//     reverse_baseline_map_1.insert(
//         ("21".to_string(), 7),
//         vec![
//             "1".to_string(),
//             "0".to_string(),
//             "1".to_string(),
//             "0".to_string(),
//             "0".to_string(),
//         ],
//     );
//     Ok((
//         forward_baseline_map_0,
//         forward_baseline_map_1,
//         reverse_baseline_map_0,
//         reverse_baseline_map_1,
//     ))
// }

// fn create_true_position_maps() -> Result<(
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
//     HashMap<String, BTreeMap<char, usize>>,
// )> {
//     let meth_pos_forward_0 = HashMap::new();
//     let meth_pos_reverse_0 = HashMap::new();
//     let mut unmeth_pos_forward_0 = HashMap::new();
//     let mut unmeth_pos_reverse_0 = HashMap::new();
//     let meth_pos_forward_1 = HashMap::new();
//     let meth_pos_reverse_1 = HashMap::new();
//     let mut unmeth_pos_forward_1 = HashMap::new();
//     let mut unmeth_pos_reverse_1 = HashMap::new();

//     let inner_unmeth_pos_forward = [('A', 0), ('C', 1), ('G', 0), ('T', 1), ('N', 0)]
//         .iter()
//         .cloned()
//         .collect();
//     unmeth_pos_forward_0.insert("21".to_string(), inner_unmeth_pos_forward);
//     let inner_unmeth_pos_forward = [('A', 0), ('C', 0), ('G', 2), ('T', 0), ('N', 0)]
//         .iter()
//         .cloned()
//         .collect();
//     unmeth_pos_forward_1.insert("21".to_string(), inner_unmeth_pos_forward);
//     let inner_unmeth_pos_reverse = [('A', 0), ('C', 2), ('G', 0), ('T', 0), ('N', 0)]
//         .iter()
//         .cloned()
//         .collect();
//     unmeth_pos_reverse_0.insert("21".to_string(), inner_unmeth_pos_reverse);
//     let inner_unmeth_pos_reverse = [('A', 1), ('C', 0), ('G', 1), ('T', 0), ('N', 0)]
//         .iter()
//         .cloned()
//         .collect();
//     unmeth_pos_reverse_1.insert("21".to_string(), inner_unmeth_pos_reverse);
//     Ok((
//         meth_pos_forward_0,
//         meth_pos_reverse_0,
//         unmeth_pos_forward_0,
//         unmeth_pos_reverse_0,
//         meth_pos_forward_1,
//         meth_pos_reverse_1,
//         unmeth_pos_forward_1,
//         unmeth_pos_reverse_1,
//     ))
// }

// fn read_bases_file(test: &str) -> Result<()> {
//     let (
//         forward_baseline_map_true_0,
//         forward_baseline_map_true_1,
//         reverse_baseline_map_true_0,
//         reverse_baseline_map_true_1,
//     ) = create_true_baseline_maps().with_context(|| "error computing true baseline maps")?;
//     let basedir = basedir(test);
//     let (
//         forward_baseline_map_0,
//         forward_baseline_map_1,
//         reverse_baseline_map_0,
//         reverse_baseline_map_1,
//     ) = deamination::assign_bases::read_bases_file(PathBuf::from(format!(
//         "{}/pos_to_bases.txt",
//         basedir
//     )))
//     .with_context(|| "error computing baseline maps")?;
//     assert_eq!(forward_baseline_map_0, forward_baseline_map_true_0);
//     assert_eq!(forward_baseline_map_1, forward_baseline_map_true_1);
//     assert_eq!(reverse_baseline_map_0, reverse_baseline_map_true_0);
//     assert_eq!(reverse_baseline_map_1, reverse_baseline_map_true_1);
//     Ok(())
// }

// fn process_bedgraph_data(test: &str) -> Result<()> {
//     let (
//         meth_pos_forward_0_true,
//         meth_pos_reverse_0_true,
//         unmeth_pos_forward_0_true,
//         unmeth_pos_reverse_0_true,
//         meth_pos_forward_1_true,
//         meth_pos_reverse_1_true,
//         unmeth_pos_forward_1_true,
//         unmeth_pos_reverse_1_true,
//     ) = create_true_position_maps().with_context(|| "error computing true baseline maps")?;

//     let (
//         forward_baseline_map_true_0,
//         forward_baseline_map_true_1,
//         reverse_baseline_map_true_0,
//         reverse_baseline_map_true_1,
//     ) = create_true_baseline_maps().with_context(|| "error computing true baseline maps")?;
//     let basedir = basedir(test);
//     let (
//         meth_pos_forward_0,
//         meth_pos_reverse_0,
//         unmeth_pos_forward_0,
//         unmeth_pos_reverse_0,
//         meth_pos_forward_1,
//         meth_pos_reverse_1,
//         unmeth_pos_forward_1,
//         unmeth_pos_reverse_1,
//     ) = deamination::assign_bases::process_bedgraph_data(
//         PathBuf::from(format!("{}/avg_bed_graph.bedGraph", basedir)),
//         &forward_baseline_map_true_0,
//         &forward_baseline_map_true_1,
//         &reverse_baseline_map_true_0,
//         &reverse_baseline_map_true_1,
//     )?;

//     assert_eq!(meth_pos_forward_0_true, meth_pos_forward_0);
//     assert_eq!(meth_pos_forward_1_true, meth_pos_forward_1);
//     assert_eq!(unmeth_pos_forward_0_true, unmeth_pos_forward_0);
//     assert_eq!(unmeth_pos_forward_1_true, unmeth_pos_forward_1);
//     assert_eq!(meth_pos_reverse_0_true, meth_pos_reverse_0);
//     assert_eq!(meth_pos_reverse_1_true, meth_pos_reverse_1);
//     assert_eq!(unmeth_pos_reverse_0_true, unmeth_pos_reverse_0);
//     assert_eq!(unmeth_pos_reverse_1_true, unmeth_pos_reverse_1);
//     Ok(())
// }

// pub fn write_pos_to_bases(
//     test: &str,
//     meth_pos_forward_0: HashMap<String, BTreeMap<char, usize>>,
//     meth_pos_reverse_0: HashMap<String, BTreeMap<char, usize>>,
//     unmeth_pos_forward_0: HashMap<String, BTreeMap<char, usize>>,
//     unmeth_pos_reverse_0: HashMap<String, BTreeMap<char, usize>>,
//     meth_pos_forward_1: HashMap<String, BTreeMap<char, usize>>,
//     meth_pos_reverse_1: HashMap<String, BTreeMap<char, usize>>,
//     unmeth_pos_forward_1: HashMap<String, BTreeMap<char, usize>>,
//     unmeth_pos_reverse_1: HashMap<String, BTreeMap<char, usize>>,
// ) -> Result<()> {
//     let basedir = basedir(test);

//     // Read true output file and get its content
//     let true_assigned_bases = format!("{}/assigned_bases.txt", basedir);
//     let file = File::open(PathBuf::from(true_assigned_bases.clone()))?;
//     let mut buf_reader = BufReader::new(file);
//     let mut true_output_content = String::new();
//     buf_reader.read_to_string(&mut true_output_content)?;

//     // Compute output of method
//     let output = format!("{}/output_test.txt", basedir);
//     cleanup_file(&output);
//     deamination::assign_bases::write_assigned_bases(
//         Some(PathBuf::from(output.clone())),
//         meth_pos_forward_0,
//         meth_pos_reverse_0,
//         unmeth_pos_forward_0,
//         unmeth_pos_reverse_0,
//         meth_pos_forward_1,
//         meth_pos_reverse_1,
//         unmeth_pos_forward_1,
//         unmeth_pos_reverse_1,
//     )
//     .with_context(|| "error computing the position counts")?;

//     // Read method output file and get its content
//     let file = File::open(PathBuf::from(output.clone()))?;
//     let mut buf_reader = BufReader::new(file);
//     let mut output_content = String::new();
//     buf_reader.read_to_string(&mut output_content)?;
//     assert_eq!(output_content, true_output_content);

//     Ok(())
// }

// #[test]
// fn test_read_bases_file() -> Result<()> {
//     read_bases_file("assign_bases")?;
//     Ok(())
// }

// #[test]
// fn test_process_bedgraph_data() -> Result<()> {
//     process_bedgraph_data("assign_bases")?;
//     Ok(())
// }

// #[test]
// fn test_write_pos_to_bases() -> Result<()> {
//     let (
//         meth_pos_forward_0,
//         meth_pos_reverse_0,
//         unmeth_pos_forward_0,
//         unmeth_pos_reverse_0,
//         meth_pos_forward_1,
//         meth_pos_reverse_1,
//         unmeth_pos_forward_1,
//         unmeth_pos_reverse_1,
//     ) = create_true_position_maps().with_context(|| "error computing true baseline maps")?;

//     write_pos_to_bases(
//         "assign_bases",
//         meth_pos_forward_0,
//         meth_pos_reverse_0,
//         unmeth_pos_forward_0,
//         unmeth_pos_reverse_0,
//         meth_pos_forward_1,
//         meth_pos_reverse_1,
//         unmeth_pos_forward_1,
//         unmeth_pos_reverse_1,
//     )?;
//     Ok(())
// }
