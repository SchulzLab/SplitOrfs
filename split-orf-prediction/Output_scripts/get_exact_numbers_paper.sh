#!/bin/bash -l



eval "$(conda shell.bash hook)"
conda activate Riboseq

so_categorization_path='/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv'
region_type='NMD'
unique_prot_orf_pairs_anno_path='/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs_annotated.txt'

python exact_numbers_paper.py \
 --so_categorization_path $so_categorization_path \
 --region_type $region_type\
 --unique_prot_orf_pairs_anno_path $unique_prot_orf_pairs_anno_path


so_categorization_path='/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/so_categorization_df.csv'
region_type='RI'
unique_prot_orf_pairs_anno_path='/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/UniqueProteinORFPairs_annotated.txt'

python exact_numbers_paper.py \
 --so_categorization_path $so_categorization_path \
 --region_type $region_type\
 --unique_prot_orf_pairs_anno_path $unique_prot_orf_pairs_anno_path