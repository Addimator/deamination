
use std::collections::HashMap;


pub fn update_base_counts(
    target_map: &mut HashMap<String, HashMap<char, usize>>,
    chrom: String,
    bases_fields: Vec<String>,
) {
    let bases_chars = vec!['A', 'C', 'G', 'T', 'N'];
    let bases_entries = target_map
        .entry(chrom.clone())
        .or_insert_with(|| HashMap::new());
    for (base, count_str) in bases_chars.iter().zip(bases_fields.iter().skip(3)) {
        let count = count_str.parse::<usize>().expect("Invalid base count");
        bases_entries.entry(base.clone())
        .and_modify(|c| *c += count)
        .or_insert(count);
    }
}