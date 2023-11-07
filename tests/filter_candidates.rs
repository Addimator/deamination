use anyhow::{Result, Context, Ok};
use std::path::PathBuf;

use deamination;


fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn filter_mutations(test: &str) -> Result<()> {
    let basedir = basedir(test);
    let candidates = deamination::find_bases::extract_vcf_positions(
        PathBuf::from(format!("{}/candidates.vcf", basedir))
    )?;
    let true_filtered_candidates: Vec<(String, u32)> = vec![
        ("21".to_string(), 7),
        ("21".to_string(), 13),
        ("21".to_string(), 15),
        ("21".to_string(), 19),
        ("21".to_string(), 25),
    ];
    // let output = format!("{}/candidates.bcf", basedir);
    // cleanup_file(&output);

    let filtered_candidates = deamination::filter_candidates::filter_mutations(
       candidates,  PathBuf::from(format!("{}/mutations.vcf", basedir))
    )
    .with_context(|| format!("error filtering mutations"))?;
    assert_eq!(
        true_filtered_candidates,
        filtered_candidates
    );
    Ok(())
}


fn filter_inconsistent(test: &str) -> Result<()> {
    let basedir = basedir(test);
    let candidates = deamination::find_bases::extract_vcf_positions(
        PathBuf::from(format!("{}/candidates.vcf", basedir))
    )?;
    let true_filtered_candidates: Vec<(String, u32)> = vec![
        ("21".to_string(), 4), ("21".to_string(), 7), ("21".to_string(), 15)
        ];
    // let output = format!("{}/candidates.bcf", basedir);
    // cleanup_file(&output);

    let filtered_candidates = deamination::filter_candidates::filter_inconsistent(
       candidates,  PathBuf::from(format!("{}/inconsistent.bed", basedir))
    )
    .with_context(|| format!("error filtering mutations"))?;
    assert_eq!(
        true_filtered_candidates,
        filtered_candidates
    );
    Ok(())
}

// fn write_filtered_candidates(test: &str) -> Result<()> {
//     let basedir = basedir(test);
//     let filtered_candidates = deamination::find_bases::extract_vcf_positions(
//         PathBuf::from(format!("{}/candidates.vcf", basedir))
//     )?;
//     // Read true output file and get its content
//     let true_output = format!("{}/candidates_out.vcf", basedir);
//     let file = File::open(PathBuf::from(true_output.clone()))?;
//     let mut buf_reader = BufReader::new(file);
//     let mut true_output_content = String::new();
//     buf_reader.read_to_string(&mut true_output_content)?;

//     // Compute output of method
//     let output = format!("{}/output_test.txt", basedir);
//     cleanup_file(&output);
//     deamination::filter_candidates::write_filtered_candidates(
//         Some(PathBuf::from(output.clone())), filtered_candidates
//     )
//     .with_context(|| format!("error computing the position counts"))?;
    
//     // Read method output file and get its content
//     let file = File::open(PathBuf::from(output.clone()))?;
//     let mut buf_reader = BufReader::new(file);
//     let mut output_content = String::new();
//     buf_reader.read_to_string(&mut output_content)?;
//     assert_eq!(output_content, true_output_content);

//     Ok(())
// }



#[test]
fn test_filter_mutations() -> Result<()> {
    filter_mutations("filter_candidates")?;
    Ok(())
}


#[test]
fn test_filter_inconsistent() -> Result<()> {
    filter_inconsistent("filter_candidates")?;
    Ok(())
}


// #[test]
// fn test_write_filtered_candidates() -> Result<()> {
//     write_filtered_candidates("filter_candidates")?;
//     Ok(())
// }
