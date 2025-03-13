use anyhow::{Ok, Result};
use find_bases::{count_bases_in_reads, process_bedgraph_data, write_assigned_bases};
use std::path::PathBuf;
use structopt::StructOpt;

mod assign_bases;
mod filter_candidates;
mod find_bases;
mod utils;

// TODO: Tests mit nicht nur matching Cigar
// TODO: Filter vcf verstehen
#[derive(Debug, StructOpt, Clone)]
pub enum Deamination {
    #[structopt(
        name = "filter-candidates",
        about = "Filter CpG positions according to mutations and consistent part from GIAB",
        usage = "cargo run -- filter-candidates candidates.bcf mutations.bcf consistent.bedGraph candidates_filtered.bcf"
    )]
    CandidateFilter {
        #[structopt(
            name = "candidates-all",
            parse(from_os_str),
            required = true,
            help = "bcf file with all CpG positions in the genome"
        )]
        candidates: PathBuf,
        #[structopt(
            name = "mutations",
            parse(from_os_str),
            required = true,
            help = "bcf file with all the positions of mutations in the genome"
        )]
        mutations: PathBuf,
        #[structopt(
            name = "consistent-bedGraph",
            parse(from_os_str),
            required = true,
            help = "Bedgraph with all the consistent parts"
        )]
        consistent: PathBuf,
        #[structopt(
            name = "output",
            parse(from_os_str),
            help = "Path of output bcf file, if not given output is printed to stdout"
        )]
        output: Option<PathBuf>,
    },
    #[structopt(
        name = "find-bases",
        about = "Find bases at CpG sites",
        usage = "cargo run -- filter-bases aligned_reads.sam candidates.bcf pos_to_bases.txt"
    )]
    BaseFinder {
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
    },
    // #[structopt(
    //     name = "assign-bases",
    //     about = "Assign bases to methylated and unmethylated cases",
    //     usage = "cargo run -- assign-bases ref.bedGraph pos_to_bases.txt output.txt"
    // )]
    // BaseAssigner {
    //     #[structopt(
    //         name = "ref-bedgraph",
    //         parse(from_os_str),
    //         required = true,
    //         help = "Bedgraph for reference which positions are (un)methylated"
    //     )]
    //     bed_graph_path: PathBuf,
    //     #[structopt(
    //         name = "pos-to-bases",
    //         parse(from_os_str),
    //         required = true,
    //         help = "Bases on each CpG position"
    //     )]
    //     bases_file_path: PathBuf,
    //     #[structopt(
    //         name = "output",
    //         parse(from_os_str),
    //         help = "Path of output TXT file, if not given output is printed to stdout"
    //     )]
    //     output: Option<PathBuf>,
    // },
}

pub fn main() -> Result<()> {
    let opt = Deamination::from_args();
    match opt {
        // Filters candidates on GIAB annotations of mutations and consistent regions
        Deamination::CandidateFilter {
            candidates,
            mutations,
            consistent,
            output,
        } => {
            // let bcf_positions = extract_bcf_positions(candidates)?;
            // let filtered_bcf_positions = filter_mutations(bcf_positions, mutations)?;
            // let filtered_bcf_positions = filter_inconsistent(filtered_bcf_positions, consistent)?;
            // write_filtered_candidates(output, filtered_bcf_positions)?;
        }
        Deamination::BaseFinder {
            sam_file_path,
            bed_graph_path,
            // bcf_file_path,
            output,
        } => {
            // let bcf_positions = extract_bcf_positions(bcf_file_path)?;
            // println!("Debug 1");
            let meth_positions = process_bedgraph_data(bed_graph_path)?;
            // println!("Debug 2");
            let position_counts = count_bases_in_reads(sam_file_path, meth_positions)?;
            // dbg!(&position_counts);
            write_assigned_bases(output, position_counts)?;
            // write_pos_to_bases(output, position_counts)?;
        } // Deamination::BaseAssigner {
          //     bases_file_path,
          //     output,
          // } => {
          //     let pos_to_info = read_bases_file(bases_file_path)?;
          // }
    }
    Ok(())
}
