use structopt::StructOpt;
use std::path::PathBuf;
use find_bases::{extract_vcf_positions, count_bases_in_reads, write_pos_to_bases};
use assign_bases::{read_bases_file, process_bedgraph_data, write_assigned_bases};
use filter_candidates::{filter_mutations, filter_inconsistent, write_filtered_candidates};
use anyhow::{Ok, Result};



mod find_bases;
mod assign_bases;
mod filter_candidates;

#[derive(Debug, StructOpt, Clone)]
pub enum Deamination {
    #[structopt(
        name = "filter-candidates",
        about = "Filter CpG positions according to mutations and consistent part from GIAB",
        usage = "cargo run -- filter-candidates candidates.vcf mutations.vcf consistent.bedGraph candidates_filtered.vcf"
    )]
    CandidateFilter {
        #[structopt(
            name = "candidates-all",
            parse(from_os_str),
            required = true,
            help = "VCF file with all CpG positions in the genome"
        )]
        candidates: PathBuf,
        #[structopt(
            name = "mutations",
            parse(from_os_str),
            required = true,
            help = "VCF file with all the positions of mutations in the genome"
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
            help = "Path of output VCF file, if not given output is printed to stdout")]
        output: Option<PathBuf>,
    },
    #[structopt(
        name = "find-bases",
        about = "Find bases at CpG sites",
        usage = "cargo run -- filter-bases aligned_reads.sam candidates.vcf pos_to_bases.txt"
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
            help = "Candidates of CpG positions in vcf format"
        )]
        vcf_file_path: PathBuf,
        #[structopt(
            name = "output", 
            parse(from_os_str), 
            help = "Path of output TXT file, if not given output is printed to stdout")]
        output: Option<PathBuf>,
    },
    #[structopt(
        name = "assign-bases",
        about = "Assign bases to methylated and unmethylated cases",
        usage = "cargo run -- assign-bases ref.bedGraph pos_to_bases.txt output.txt"
    )]
    BaseAssigner {
        #[structopt(
            name = "ref-bedgraph",
            parse(from_os_str),
            required = true,
            help = "Bedgraph for reference which positions are (un)methylated"
        )]
        bed_graph_path: PathBuf,
        #[structopt(
            name = "pos-to-bases",
            parse(from_os_str),
            required = true,
            help = "Bases on each CpG position"
        )]
        bases_file_path: PathBuf,
        #[structopt(
            name = "output", 
            parse(from_os_str), 
            help = "Path of output TXT file, if not given output is printed to stdout")]
        output: Option<PathBuf>,
    },
}




pub fn main() -> Result<()> {
    let opt = Deamination::from_args();
    match opt {
        Deamination::CandidateFilter { candidates, mutations, consistent, output } => {
            let vcf_positions = extract_vcf_positions(candidates)?;
            let filtered_vcf_positions = filter_mutations(vcf_positions, mutations)?;
            let filtered_vcf_positions = filter_inconsistent(filtered_vcf_positions, consistent)?;       
            write_filtered_candidates(output, filtered_vcf_positions)?;
        }
        Deamination::BaseFinder { sam_file_path, vcf_file_path, output } => {
            let vcf_positions = extract_vcf_positions(vcf_file_path)?;
            let position_counts = count_bases_in_reads(sam_file_path, &vcf_positions)?;
            write_pos_to_bases(output, position_counts)?;
        }
        Deamination::BaseAssigner { bed_graph_path, bases_file_path, output } => {
            let (forward_baseline_map_0, forward_baseline_map_1, reverse_baseline_map_0, reverse_baseline_map_1) = read_bases_file(bases_file_path)?;
            let (meth_pos_forward_0, meth_pos_reverse_0, unmeth_pos_forward_0, unmeth_pos_reverse_0, 
                meth_pos_forward_1, meth_pos_reverse_1, unmeth_pos_forward_1, unmeth_pos_reverse_1) = process_bedgraph_data(bed_graph_path, &forward_baseline_map_0, &forward_baseline_map_1, &reverse_baseline_map_0, &reverse_baseline_map_1)?;
            write_assigned_bases(output, 
                meth_pos_forward_0, meth_pos_reverse_0, unmeth_pos_forward_0, unmeth_pos_reverse_0, 
                meth_pos_forward_1, meth_pos_reverse_1, unmeth_pos_forward_1, unmeth_pos_reverse_1)?;
        }
    }
    Ok(())
}

