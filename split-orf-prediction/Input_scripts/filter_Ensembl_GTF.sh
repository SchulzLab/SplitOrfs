#!/bin/bash


# =============================================================================
# Script Name: filter_Ensembl_GTF.sh
# Description: This script filters the Ensembl GTF annotation 110 for 
#               TSL support level 1 or 2 or equality of the intron chain with RefSeq:
#              - Step 1: Filter GTF for TSL 1,2 (correct, also if TSL 1 (...) or 2 (...))
#              - Step 2: GFFCompare with RefSeq and filter GTF
#              - Step 3: final filtering of Ensembl GTf by transcripts in either 
#                       TSL 1,2 or RefSeq equality filtered GTF
# Usage:       bash filter_Ensembl_GTF.sh
# Author:      Christina Kalk
# Date:        2025-09-29
# =============================================================================

ENSEMBL_GTF="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
REFSEQ_GTF="/projects/splitorfs/work/reference_files/GRCh38_latest_genomic.gtf"
REF_DIR=$(dirname $REFSEQ_GTF)

OUT_DIR="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25"

SCRIPT_DIR="/home/ckalk/scripts/Reference_scripts"


if [[ ! -d "$OUT_DIR" ]]; then
    mkdir $OUT_DIR
fi


# ----- initialize conda ----- #
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pygtftk



python ${SCRIPT_DIR}/clean_ensembl_gtf_tsl1_and_2_correct_29_09_25.py \
    $ENSEMBL_GTF \
    ${OUT_DIR}/Ensembl_110_filtered_transcripts_tsl1_and_2_correct_29_09_25.gtf \
    ${OUT_DIR}/Ensembl_110_filtered_genes_tsl1_and_2_correct_29_09_25.gtf


# combine genes with transcripts
cat ${OUT_DIR}/Ensembl_110_filtered_transcripts_tsl1_and_2_correct_29_09_25.gtf \
    ${OUT_DIR}/Ensembl_110_filtered_genes_tsl1_and_2_correct_29_09_25.gtf \
    > ${OUT_DIR}/Ensembl_110_filtered_tsl1_and_2_correct_29_09_25.gtf



python ${SCRIPT_DIR}/filter_out_gene_entries.py \
 ${REFSEQ_GTF} \
 ${REF_DIR}/$(basename $REFSEQ_GTF .gtf)_no_genes.gtf 

python ${SCRIPT_DIR}/get_chr_names_for_refseq.py \
 ${REF_DIR}/$(basename $REFSEQ_GTF .gtf)_no_genes.gtf \
 ${REF_DIR}/GRCh38_latest_assembly_report.txt \
 ${REF_DIR}/GRCh38_latest_genomic_no_genes_chrs.gtf


conda activate pacbio
gffcompare -o gffcompare_refseq -r ${REF_DIR}/GRCh38_latest_genomic_no_genes_chrs.gtf ${ENSEMBL_GTF}

mv ${REF_DIR}/gffcompare_refseq.Homo_sapiens.GRCh38.110.chr.gtf.refmap \
 ${OUT_DIR}/gffcompare_refseq.Homo_sapiens.GRCh38.110.chr.gtf.refmap

mv ${REF_DIR}/gffcompare_refseq.Homo_sapiens.GRCh38.110.chr.gtf.tmap \
 ${OUT_DIR}/gffcompare_refseq.Homo_sapiens.GRCh38.110.chr.gtf.tmap

conda activate pygtftk

python ${SCRIPT_DIR}/filter_by_equality.py ${ENSEMBL_GTF} ${OUT_DIR}/gffcompare_refseq.Homo_sapiens.GRCh38.110.chr.gtf.tmap \
    ${OUT_DIR}/Ensembl_equality_filtered.gtf



# decided to use the two GTFs to get transcript IDs
# take all exons, CDS and transcripts plus corresponding genes in separate GTF
# concat afterwards
python ${SCRIPT_DIR}/filter_by_equality_and_TSL_GTF_approach_29_09_25.py \
    ${ENSEMBL_GTF} \
    ${OUT_DIR}/Ensembl_110_filtered_tsl1_and_2_correct_29_09_25.gtf \
    ${OUT_DIR}/Ensembl_equality_filtered.gtf \
    ${OUT_DIR}/Ensembl_110_filtered_equality_and_tsl1_2_transcripts_correct_29_09_25.gtf \
    ${OUT_DIR}/Ensembl_110_filtered_equality_and_tsl1_2_genes_correct_29_09_25.gtf


cat ${OUT_DIR}/Ensembl_110_filtered_equality_and_tsl1_2_transcripts_correct_29_09_25.gtf \
    ${OUT_DIR}/Ensembl_110_filtered_equality_and_tsl1_2_genes_correct_29_09_25.gtf \
    > ${OUT_DIR}/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf