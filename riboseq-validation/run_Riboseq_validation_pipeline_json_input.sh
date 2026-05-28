#----- This script maps Ribo-seq data to the genome using the supplied annotation ----- #
# ----- then an intersection with unique regions in genome coords from the split-ORF pipeline ----- #
# ----- is performed as well as background regions of 3' UTRs, an empirical  ----- #
# ----- background distribution is used to determine which unique regions are  ----- #
# ----- and this is summarized in an Rmd report ----- #

#!/bin/bash -l

CONFIG=$1
echo "CONFIG is: $CONFIG"
ls -l "$CONFIG"

output_star=$(jq -r '.output_star' "$CONFIG")
unique_region_dir=$(jq -r '.unique_region_dir' "$CONFIG")
ensembl_gtf=$(jq -r '.ensembl_gtf' "$CONFIG")
genome_fasta=$(jq -r '.genome_fasta' "$CONFIG")
input_fastq_path=$(jq -r '.input_fastq_path // empty' "$CONFIG")
three_primes=$(jq -r '.three_primes' "$CONFIG")
cds_coordinates=$(jq -r '.cds_coordinates' "$CONFIG")
output_dir=$(jq -r '.output_dir' "$CONFIG")
# script_path=$(jq -r '.script_path' "$CONFIG")
ribo_pipe_path=$(jq -r '.ribo_pipe_path' "$CONFIG")
input_name=$(jq -r '.input_name' "$CONFIG")
region_type=$(jq -r '.region_type' "$CONFIG")
bam_ending=$(jq -r '.bam_ending' "$CONFIG")
tmp_dir=$(jq -r '.tmp_dir' "$CONFIG")
report=$(jq -r '.report // empty' "$CONFIG")
duplicated=$(jq -r '.duplicated' "$CONFIG")

script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"/ribo_cov_scripts

# requirement: BAM files need end with: bam_ending!!!!!

# output_star="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter"
# unique_region_dir="/home/ckalk/tools/SplitORF_pipeline/Output/run_13.01.2026-11.44.16_NMD_CDS_subtraction_minAlignLength_15"
# ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
# filtered_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# input_fastq_path=/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples
# three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged_numbered.bed"
# upf10_dedup_dir="/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated"
# cds_coordinates="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"
# # script_path="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random/expressed_genes_multiple_testing_approach_November_2025"
# ribo_pipe_path="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline"

# input_name="NMD"
# region_type="NMD"

bash "${script_path}"/run_Riboseq_validation_pipeline.sh \
    -b "$bam_ending" \
    -c "$cds_coordinates" \
    -d "$duplicated" \
    -e "$ensembl_gtf" \
    -g "$genome_fasta" \
    -i "$input_fastq_path" \
    -n "$input_name" \
    -o "$output_star" \
    -p "$tmp_dir" \
    -r "$region_type" \
    -s "$script_path" \
    -t "$three_primes" \
    -u "$unique_region_dir"
   


 

if [ -n "$report" ]; then
    export script_path
    export output_star
    export region_type
    export unique_region_dir
    export report
    R -e 'library(rmarkdown); rmarkdown::render(input = file.path(Sys.getenv("script_path"), "RiboSeqReportGenomic_iteration_update_expression_filter_multiple_test_correction.Rmd"), output_file = Sys.getenv("report"), params=list(args = c(file.path(Sys.getenv("output_star"), paste(Sys.getenv("region_type"), "genome", sep = "_")), Sys.getenv("unique_region_dir"), Sys.getenv("region_type"), Sys.getenv("script_path"))))'
fi


