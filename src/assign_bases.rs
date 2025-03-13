// use std::collections::{HashMap, BTreeMap};
// use std::fs::File;
// use std::io::{Write, BufReader, BufWriter};
// use std::path::PathBuf;
// use anyhow::{Result, Context};
// use std::io::prelude::*;
// use std::collections::HashSet;

// use crate::utils::{Position, MethPos, PosType} ;

// /// Computes for the different positions f_0, r_0, f_1, r_1 the respective number of bases at these positions
// pub fn read_bases_file(bases_file_path: PathBuf) -> Result<HashMap<Position, Vec<MethPos>>> {
//     let mut meth_positions: HashMap<Position, Vec<MethPos>> = HashMap::new();

//     let bases_file = File::open(bases_file_path)
//         .with_context(|| format!("Unable to open bases file"))?;
//     let bases_reader = BufReader::new(bases_file);

//     for base_line in bases_reader.lines() {
//         let base_line = base_line?;
//         if base_line.starts_with("#") {
//             continue;
//         }
//         let base_fields: Vec<String> = base_line.split('\t').map(|s| s.to_string()).collect();

//         let chrom = base_fields[0].to_string();
//         let pos = base_fields[1].parse::<u32>().context("Invalid position value")?;
//         let direction = base_fields[2].to_string();
//         let position = Position::new(chrom.clone(), pos, Some(direction.clone()));
//         let position_no_direction = Position::new(chrom.clone(), pos, None);
//         // TODO: Methylation value is not used in the current implementation, should be hashmap read from file
//         println!("{:?}",  base_fields[3..].to_vec());
//         let meth_pos = MethPos::new(position.clone(), None, HashMap::new());
//         meth_positions.entry(position_no_direction).or_default().push(meth_pos);

//     }

//     Ok(meth_positions)
// }

// pub fn process_bedgraph_data(bed_graph_path: PathBuf, mut meth_positions: HashMap<Position, Vec<MethPos>>)
//  ->  Result<HashMap<PosType, Vec<MethPos>>> {
//     let bedgraph_file = File::open(bed_graph_path)
//         .with_context(|| format!("Unable to open bedGraph file."))?;
//     let bedgraph_reader = BufReader::new(bedgraph_file);
//     let mut type_to_pos: HashMap<PosType, Vec<MethPos>> = HashMap::new();

//     for bed_line in bedgraph_reader.lines() {
//         let bed_line = bed_line?;
//         if bed_line.starts_with("track") {
//             continue;
//         }
//         let bed_fields: Vec<&str> = bed_line.split('\t').collect();

//         // Extrahiere relevante Informationen aus dem Bedgraph
//         let chrom_orig = bed_fields[0].to_string();
//         let chrom = if let Some(s) = chrom_orig.strip_prefix("chr") {
//             s.to_string()
//         } else {
//             chrom_orig
//         };
//         // TODO: Was genau ist Startposition
//         let pos = [bed_fields[1].parse::<usize>().context("Invalid position value")?, bed_fields[2].parse::<usize>().context("Invalid position value")?];
//         let methylation = bed_fields[3].parse::<f64>().context("Invalid methylation value")?;
//         let position = Position::new(chrom.clone(), pos[0] as u32, None);
//         if meth_positions.contains_key(&position) {
//             for meth_position in meth_positions[&position].iter_mut() {
//                 meth_position.set_methylation(Some(methylation));
//                 type_to_pos.entry(PosType::new(meth_position.position().direction().as_ref().unwrap(),  methylation > 20.0)).or_default().push(meth_position.to_owned());

//             }
//         }
//     }
//     Ok(type_to_pos)
// }

// fn bases_percentage(target_map: &HashMap<String, BTreeMap<char, usize>>, chrom: &String, base: char) -> f64 {
//     if let Some(map_chrom) = target_map.get(chrom) {
//         let sum_bases: usize = map_chrom.values().sum();
//         if let Some(&number_base) = map_chrom.get(&base) {
//             if sum_bases > 0 {
//                 println!("{:?}, {:?}", number_base, sum_bases);
//                 return (number_base as f64 / sum_bases as f64) * 100.0;
//             }
//         }
//     }
//     0.0 // Return 0.0 if something goes wrong or if sum_bases is 0.
// }

// pub fn write_assigned_bases(output: Option<PathBuf>, meth_positions: HashMap<PosType, Vec<MethPos>>) -> Result<()>{

//     let mut writer: Box<dyn Write> = match output {
//         Some(path) => {
//             let file = std::fs::File::create(path.clone())
//                 .with_context(|| format!("error creating file: {:?}", path))?;
//             Box::new(BufWriter::new(file))
//         }
//         None => Box::new(BufWriter::new(std::io::stdout())),
//     };

//     for dir in ["f_0", "r_0", "f_1", "r_1"].iter() {
//         for meth_typ in ["meth", "unmeth"].iter() {
//             // TODO: metype is &&str not bool
//             for meth_pos in meth_positions[&PosType::new(dir.to_string(), meth_typ)].iter() {
//                 // Write methylation info per position (Debug?)
//                 writeln!(writer, "meth_pos{:?}", meth_pos)?;

//             }
//             // Write whole methylation info for each direction
//         }
//     }
//     writeln!(writer)?;
//         // Erstelle ein HashSet, um Schlüssel ohne Wiederholungen zu speichern
//     let mut unique_keys: HashSet<String> = HashSet::new();

//     // Füge die Schlüssel aus allen HashMaps zum HashSet hinzu
//     for map in &[&meth_pos_forward_0, &meth_pos_reverse_0, &unmeth_pos_forward_0, &unmeth_pos_reverse_0, &meth_pos_forward_1, &meth_pos_reverse_1, &unmeth_pos_forward_1, &unmeth_pos_reverse_1] {
//         for key in map.keys() {
//             unique_keys.insert(key.clone());
//         }
//     }

//     for chrom in unique_keys {
//         writeln!(writer,
//             "0: Methylated positions
//             As in forward string: {} \tCs in forward string: {} \tGs in forward string: {} \tTs IN FORWARD STRING: {}
//             As in reverse String: {} \tCs in reverse String: {} \tGs in reverse String: {} \tTs in reverse String: {}",
//             bases_percentage(&meth_pos_forward_0, &chrom, 'A'), bases_percentage(&meth_pos_forward_0, &chrom, 'C'), bases_percentage(&meth_pos_forward_0, &chrom, 'G'), bases_percentage(&meth_pos_forward_0, &chrom, 'T'),
//             bases_percentage(&meth_pos_reverse_0, &chrom, 'A'), bases_percentage(&meth_pos_reverse_0, &chrom, 'C'), bases_percentage(&meth_pos_reverse_0, &chrom, 'G'), bases_percentage(&meth_pos_reverse_0, &chrom, 'T'))?;
//         writeln!(writer,
//             "1: Methylated positions
//             As in forward string: {} \tCs in forward string: {} \tGs in forward string: {} \tTs in forward string: {}
//             As IN REVERSE STRING: {} \tCs in reverse String: {} \tGs in reverse String: {} \tTs in reverse String: {}",
//             bases_percentage(&meth_pos_forward_1, &chrom, 'A'), bases_percentage(&meth_pos_forward_1, &chrom, 'C'), bases_percentage(&meth_pos_forward_1, &chrom, 'G'), bases_percentage(&meth_pos_forward_1, &chrom, 'T'),
//             bases_percentage(&meth_pos_reverse_1, &chrom, 'A'), bases_percentage(&meth_pos_reverse_1, &chrom, 'C'), bases_percentage(&meth_pos_reverse_1, &chrom, 'G'), bases_percentage(&meth_pos_reverse_1, &chrom, 'T'))?;
//         writeln!(writer,
//             "0: Unmethylated positions
//             As in forward string: {} \tCs in forward string: {} \tGs in forward string: {} \tTs IN FORWARD STRING: {}
//             As in reverse String: {} \tCs in reverse String: {} \tGs in reverse String: {} \tTs in reverse String: {}",
//             bases_percentage(&unmeth_pos_forward_0, &chrom, 'A'), bases_percentage(&unmeth_pos_forward_0, &chrom, 'C'), bases_percentage(&unmeth_pos_forward_0, &chrom, 'G'), bases_percentage(&unmeth_pos_forward_0, &chrom, 'T'),
//             bases_percentage(&unmeth_pos_reverse_0, &chrom, 'A'), bases_percentage(&unmeth_pos_reverse_0, &chrom, 'C'), bases_percentage(&unmeth_pos_reverse_0, &chrom, 'G'), bases_percentage(&unmeth_pos_reverse_0, &chrom, 'T'))?;
//         writeln!(writer,
//             "1: Unmethylated positions
//             As in forward string: {} \tCs in forward string: {} \tGs in forward string: {} \tTs in forward string: {}
//             As IN REVERSE STRING: {} \tCs in reverse String: {} \tGs in reverse String: {} \tTs in reverse String: {}",
//             bases_percentage(&unmeth_pos_forward_1, &chrom, 'A'), bases_percentage(&unmeth_pos_forward_1, &chrom, 'C'), bases_percentage(&unmeth_pos_forward_1, &chrom, 'G'), bases_percentage(&unmeth_pos_forward_1, &chrom, 'T'),
//             bases_percentage(&unmeth_pos_reverse_1, &chrom, 'A'), bases_percentage(&unmeth_pos_reverse_1, &chrom, 'C'), bases_percentage(&unmeth_pos_reverse_1, &chrom, 'G'), bases_percentage(&unmeth_pos_reverse_1, &chrom, 'T'))?;
//     }
//     Ok(())
// }
