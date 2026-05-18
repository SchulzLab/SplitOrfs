#!/bin/bash


# =============================================================================
# Script Name: generate_CDS_contamination_subtraction_coordinates.sh
# Description: generates a BED file of all contamination RNAs as well as
#               CDS coordinates which will be subtracted from unqiue regions
# Usage:       bash generate_CDS_contamination_subtraction_coordinates.sh
# Author:      Christina Kalk
# Date:        2025-09-29
# =============================================================================

input_dir="/home/ckalk/tools/SplitORF_pipeline/Input2023"
cds_coordinates="${input_dir}/TSL_eq_filtered_29_09_25/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"
contamination_dir="${input_dir}/intronic_RNAs"


# ----- initialize conda ----- #
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pygtftk

# If line number within current file (FNR) is one, but it is not the first file
# then skip: skipping of unnecessary Ensmebl headers 
awk 'FNR==1 && NR!=1 {next} {print}' "${contamination_dir}"/*.txt > "${contamination_dir}"/all_contaminants.txt


python make_contaminant_bed_from_Ensembl.py \
 --ensembl_genomic_coords_txt "${contamination_dir}"/all_contaminants.txt \
 --output_bed_file "${contamination_dir}"/all_contaminants.bed

cat "${contamination_dir}"/all_contaminants.bed "${cds_coordinates}" > "${input_dir}"/CDS_110_filtered_with_contaminants.bed



