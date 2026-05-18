#!/bin/bash


# =============================================================================
# Script Name: 
# Description: Filter the protein coding prot and cDNA sequences
# Usage:       bash 
# Author:      Christina Kalk
# Date:        2025-09-29
# =============================================================================


# ----- initialize conda ----- #
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pygtftk


cd /home/ckalk/tools/SplitORF_pipeline/Input_scripts
python filter_fasta_by_gtf.py \
/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf \
/home/ckalk/tools/SplitORF_pipeline/Input2023/protein_coding_peptide_sequences.fa \
/home/ckalk/tools/SplitORF_pipeline/Input2023/TSL_eq_filtered_29_09_25/protein_coding_peptide_sequences_tsl_eq_filtered_29_09_25.fa


python filter_fasta_by_gtf.py \
/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf \
/home/ckalk/tools/SplitORF_pipeline/Input2023/protein_coding_transcript_and_gene_cDNA.fa \
/home/ckalk/tools/SplitORF_pipeline/Input2023/TSL_eq_filtered_29_09_25/protein_coding_transcript_and_gene_cDNA_tsl_eq_filtered_29_09_25.fa
