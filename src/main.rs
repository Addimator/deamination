use anyhow::{Ok, Result};
use find_bases::{count_bases_in_reads, process_bedgraph_data, write_assigned_bases};
use std::path::PathBuf;
use structopt::StructOpt;

mod find_bases;
mod utils;

// How does this work?
// 1. Read the bedGraph file (bedgraph file because we use the average bedgraph in order  to find methylation info of each position) and extract the positions of interest.
//     Create for each position a MethPos object with the position and methylation level.
// 2. Read the BAM file and extract the bases at the positions of interest. Insert the number of bases in the meth_bases field of the MethPos object.
// 3. Write the bases to a file or print them to stdout.
// 4. The output file will contain the bases at the positions of interest in the format:
//    <chromosome> <position> <base1> <base2> ... <baseN>

#[derive(Debug, StructOpt)]
struct FindBasesArgs {
    #[structopt(
        name = "aligned-reads",
        parse(from_os_str),
        required = true,
        help = "Aligned reads in SAM format"
    )]
    sam_file_path: PathBuf,

    #[structopt(
        name = "candidates",
        parse(from_os_str),
        required = true,
        help = "Bedgraph with methylation level of CpG positions"
    )]
    bed_graph_path: PathBuf,

    #[structopt(
        name = "output",
        parse(from_os_str),
        help = "Path of output TXT file, if not given output is printed to stdout"
    )]
    output: Option<PathBuf>,
}

fn find_bases(
    sam_file_path: PathBuf,
    bed_graph_path: PathBuf,
    output: Option<PathBuf>,
) -> Result<()> {
    println!("Step(1/3): Finding positions in bedGraph file");
    let meth_positions = process_bedgraph_data(bed_graph_path)?;

    println!("Step(2/3): Finding nucleotides at positions under interest");
    let position_counts = count_bases_in_reads(sam_file_path, meth_positions)?;

    println!("Step(3/3): Writing output to file");
    write_assigned_bases(output, position_counts)?;

    Ok(())
}

fn main() -> Result<()> {
    let args = FindBasesArgs::from_args();
    find_bases(args.sam_file_path, args.bed_graph_path, args.output)
}
