
use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{Write, BufReader, BufWriter};
use std::path::PathBuf;
use anyhow::{Result, Context};
use std::io::prelude::*;
use std::collections::HashSet;




pub fn 

read_bases_file(bases_file_path: PathBuf) -> Result<(HashMap<(String, usize), Vec<String>>, HashMap<(String, usize), Vec<String>>, HashMap<(String, usize), Vec<String>>, HashMap<(String, usize), Vec<String>>)> {
    let mut forward_baseline_map_0: HashMap<(String, usize), Vec<String>> = HashMap::new();
    let mut forward_baseline_map_1: HashMap<(String, usize), Vec<String>> = HashMap::new();
    let mut reverse_baseline_map_0: HashMap<(String, usize), Vec<String>> = HashMap::new();
    let mut reverse_baseline_map_1: HashMap<(String, usize), Vec<String>> = HashMap::new();
    
    let bases_file = File::open(bases_file_path)
        .with_context(|| format!("Unable to open bases file"))?;
    let bases_reader = BufReader::new(bases_file);

    for base_line in bases_reader.lines() {
        let base_line = base_line?;
        if base_line.starts_with("#") {
            continue;
        }
        let base_fields: Vec<String> = base_line.split('\t').map(|s| s.to_string()).collect();

        let chrom = base_fields[0].to_string();
        let position = base_fields[1].parse::<usize>().context("Invalid position value")?;
        let direction = base_fields[2].to_string();

        if direction == "f_0" {
            forward_baseline_map_0.insert((chrom.clone(), position), base_fields[3..].to_vec());
        } 
        if direction == "f_1" {
            forward_baseline_map_1.insert((chrom.clone(), position), base_fields[3..].to_vec());
        } 
        if direction == "r_0" {
            reverse_baseline_map_0.insert((chrom.clone(), position), base_fields[3..].to_vec());
        } 
        if direction == "r_1" {
            reverse_baseline_map_1.insert((chrom.clone(), position), base_fields[3..].to_vec());
        }

    }

    Ok((forward_baseline_map_0, forward_baseline_map_1, reverse_baseline_map_0 ,reverse_baseline_map_1))
}

pub fn process_bedgraph_data(bed_graph_path: PathBuf, forward_baseline_map_0: &HashMap<(String, usize), Vec<String>>,forward_baseline_map_1: &HashMap<(String, usize), Vec<String>>, reverse_baseline_map_0: &HashMap<(String, usize), Vec<String>>, reverse_baseline_map_1: &HashMap<(String, usize), Vec<String>>)
 -> Result<(HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>, HashMap<String, BTreeMap<char, usize>>)> {
    let bedgraph_file = File::open(bed_graph_path)
        .with_context(|| format!("Unable to open bedGraph file."))?;
    let bedgraph_reader = BufReader::new(bedgraph_file);

    let mut meth_pos_forward_0: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut meth_pos_reverse_0: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut unmeth_pos_forward_0: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut unmeth_pos_reverse_0: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut meth_pos_forward_1: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut meth_pos_reverse_1: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut unmeth_pos_forward_1: HashMap<String, BTreeMap<char, usize>> = HashMap::new();
    let mut unmeth_pos_reverse_1: HashMap<String, BTreeMap<char, usize>> = HashMap::new();

    for bed_line in bedgraph_reader.lines() {
        let bed_line = bed_line?;
        if bed_line.starts_with("track") {
            continue;
        }
        let bed_fields: Vec<&str> = bed_line.split('\t').collect();

        // Extrahiere relevante Informationen aus dem Bedgraph
        let chrom_orig = bed_fields[0].to_string();
        let chrom = if let Some(s) = chrom_orig.strip_prefix("chr") {
            s.to_string()
        } else {
            chrom_orig
        };
        let position = [bed_fields[1].parse::<usize>().context("Invalid position value")?, bed_fields[2].parse::<usize>().context("Invalid position value")?];
        let methylation = bed_fields[3].parse::<f64>().context("Invalid methylation value")?;

        if forward_baseline_map_0.contains_key(&(chrom.clone(), position[1])) {
            let bases_fields_forward = &forward_baseline_map_0[&(chrom.clone(), position[1])];
            if methylation > 20.0 {
                update_base_counts(&mut meth_pos_forward_0, chrom.clone(), bases_fields_forward.to_vec());
            } else {
                update_base_counts(&mut unmeth_pos_forward_0, chrom.clone(), bases_fields_forward.to_vec());
            }
        }
        if forward_baseline_map_1.contains_key(&(chrom.clone(), position[1])) {
            let bases_fields_forward = &forward_baseline_map_1[&(chrom.clone(), position[1])];
            if methylation > 20.0 {
                update_base_counts(&mut meth_pos_forward_1, chrom.clone(), bases_fields_forward.to_vec());
            } else {
                update_base_counts(&mut unmeth_pos_forward_1, chrom.clone(), bases_fields_forward.to_vec());
            }
        }
        if reverse_baseline_map_0.contains_key(&(chrom.clone(), position[0])) {
            let bases_fields_reverse = &reverse_baseline_map_0[&(chrom.clone(), position[0])];
            if methylation > 20.0 {
                update_base_counts(&mut meth_pos_reverse_0, chrom.clone(), bases_fields_reverse.to_vec());
            } else {
                update_base_counts(&mut unmeth_pos_reverse_0, chrom.clone(), bases_fields_reverse.to_vec());
            }
        }
        if reverse_baseline_map_1.contains_key(&(chrom.clone(), position[0])) {
            let bases_fields_reverse = &reverse_baseline_map_1[&(chrom.clone(), position[0])];
            if methylation > 20.0 {
                update_base_counts(&mut meth_pos_reverse_1, chrom.clone(), bases_fields_reverse.to_vec());
            } else {
                update_base_counts(&mut unmeth_pos_reverse_1, chrom.clone(), bases_fields_reverse.to_vec());
            }
        }
    }

    Ok((meth_pos_forward_0, meth_pos_reverse_0, unmeth_pos_forward_0, unmeth_pos_reverse_0, 
        meth_pos_forward_1, meth_pos_reverse_1, unmeth_pos_forward_1, unmeth_pos_reverse_1))
}


fn update_base_counts(
    target_map: &mut HashMap<String, BTreeMap<char, usize>>,
    chrom: String,
    bases_fields: Vec<String>,
) {
    let bases_chars = vec!['A', 'C', 'G', 'T', 'N'];
    let bases_entries = target_map
        .entry(chrom.clone())
        .or_insert_with(|| BTreeMap::new());
    for (base, count_str) in bases_chars.iter().zip(bases_fields.iter()) {
        let count = count_str.parse::<usize>().expect("Invalid base count");
        bases_entries.entry(base.clone())
        .and_modify(|c| *c += count)
        .or_insert(count);
    }
}


fn bases_percentage(target_map: &HashMap<String, BTreeMap<char, usize>>, chrom: &String, base: char) -> f64 {
    if let Some(map_chrom) = target_map.get(chrom) {
        let sum_bases: usize = map_chrom.values().sum();
        if let Some(&number_base) = map_chrom.get(&base) {
            if sum_bases > 0 {
                println!("{:?}, {:?}", number_base, sum_bases);
                return (number_base as f64 / sum_bases as f64) * 100.0;
            }
        }
    }
    0.0 // Return 0.0 if something goes wrong or if sum_bases is 0.
}

pub fn write_assigned_bases(output: Option<PathBuf>, meth_pos_forward_0:HashMap<String, BTreeMap<char, usize>>, meth_pos_reverse_0:HashMap<String, BTreeMap<char, usize>>, unmeth_pos_forward_0:HashMap<String, BTreeMap<char, usize>>, unmeth_pos_reverse_0: HashMap<String, BTreeMap<char, usize>>, meth_pos_forward_1:HashMap<String, BTreeMap<char, usize>>, meth_pos_reverse_1:HashMap<String, BTreeMap<char, usize>>, unmeth_pos_forward_1:HashMap<String, BTreeMap<char, usize>>, unmeth_pos_reverse_1: HashMap<String, BTreeMap<char, usize>>) -> Result<()>{
            
    let mut writer: Box<dyn Write> = match output {
        Some(path) => {
            let file = std::fs::File::create(path.clone())
                .with_context(|| format!("error creating file: {:?}", path))?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(std::io::stdout())),
    };


    writeln!(writer, "meth_pos_forward Pos 0: {:?}", meth_pos_forward_0)?;
    writeln!(writer, "meth_pos_reverse Pos 0: {:?}", meth_pos_reverse_0)?;
    writeln!(writer, "unmeth_pos_forward Pos 0: {:?}", unmeth_pos_forward_0)?;
    writeln!(writer, "unmeth_pos_reverse Pos 0: {:?}", unmeth_pos_reverse_0)?;
    writeln!(writer, "meth_pos_forward Pos 1: {:?}", meth_pos_forward_1)?;
    writeln!(writer, "meth_pos_reverse Pos 1: {:?}", meth_pos_reverse_1)?;
    writeln!(writer, "unmeth_pos_forward Pos 1: {:?}", unmeth_pos_forward_1)?;
    writeln!(writer, "unmeth_pos_reverse Pos 1: {:?}", unmeth_pos_reverse_1)?;
    writeln!(writer)?;


    // writeln!(writer, "meth_pos_forward Pos 0: {:?}, Sum: {:?}", meth_pos_forward_0, meth_pos_forward_0.values().sum())?;
    // writeln!(writer, "meth_pos_reverse Pos 0: {:?}, Sum: {:?}", meth_pos_reverse_0, meth_pos_reverse_0.values().sum())?;
    // writeln!(writer, "unmeth_pos_forward Pos 0: {:?}, Sum: {:?}", unmeth_pos_forward_0, unmeth_pos_forward_0.values().sum())?;
    // writeln!(writer, "unmeth_pos_reverse Pos 0: {:?}, Sum: {:?}", unmeth_pos_reverse_0, unmeth_pos_reverse_0.values().sum())?;
    // writeln!(writer, "meth_pos_forward Pos 1: {:?}, Sum: {:?}", meth_pos_forward_1, meth_pos_forward_1.values().sum())?;
    // writeln!(writer, "meth_pos_reverse Pos 1: {:?}, Sum: {:?}", meth_pos_reverse_1, meth_pos_reverse_1.values().sum())?;
    // writeln!(writer, "unmeth_pos_forward Pos 1: {:?}, Sum: {:?}", unmeth_pos_forward_1, unmeth_pos_forward_1.values().sum())?;
    // writeln!(writer, "unmeth_pos_reverse Pos 1: {:?}, Sum: {:?}", unmeth_pos_reverse_1, unmeth_pos_reverse_1.values().sum())?;
    // writeln!(writer)?;
    // Erstelle ein HashSet, um Schlüssel ohne Wiederholungen zu speichern
    let mut unique_keys: HashSet<String> = HashSet::new();

    // Füge die Schlüssel aus allen HashMaps zum HashSet hinzu
    for map in &[&meth_pos_forward_0, &meth_pos_reverse_0, &unmeth_pos_forward_0, &unmeth_pos_reverse_0, &meth_pos_forward_1, &meth_pos_reverse_1, &unmeth_pos_forward_1, &unmeth_pos_reverse_1] {
        for key in map.keys() {
            unique_keys.insert(key.clone());
        }
    }

    for chrom in unique_keys {
        writeln!(writer, "0: Methylated positions \n\tTs in forward string: {} \n\tAs in reverse String: {}", bases_percentage(&meth_pos_forward_0, &chrom, 'T'), bases_percentage(&meth_pos_reverse_0, &chrom, 'A'))?;
        writeln!(writer, "1: Methylated positions \n\tTs in forward string: {} \n\tAs in reverse String: {}", bases_percentage(&meth_pos_forward_1, &chrom, 'T'), bases_percentage(&meth_pos_reverse_1, &chrom, 'A'))?;
        writeln!(writer, "0: Unmethylated positions \n\tTs in forward string: {} \n\tAs in reverse String: {}", bases_percentage(&unmeth_pos_forward_0, &chrom, 'T'), bases_percentage(&unmeth_pos_reverse_0, &chrom, 'A'))?;
        writeln!(writer, "1: Unmethylated positions \n\tTs in forward string: {} \n\tAs in reverse String: {}", bases_percentage(&unmeth_pos_forward_1, &chrom, 'T'), bases_percentage(&unmeth_pos_reverse_1, &chrom, 'A'))?;
    }
    Ok(())
}

