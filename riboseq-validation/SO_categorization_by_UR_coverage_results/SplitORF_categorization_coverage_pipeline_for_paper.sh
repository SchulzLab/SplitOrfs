#!/bin/bash -l
# this script runs the downstream analysis steps for the validated regions per sample
# teh NMd inhibited samples are treated separately, for these validation only if 2 ORFs are Ribo-covered
# 1. get all validated regions union and statistics for NMD and RI
# 2. check for RBP validation from RBPDB: the database needs to be downloaded and paths 
# indicated for this to run


eval "$(conda shell.bash hook)"
conda activate Riboseq

rbpbase=true
pfam=false

nmd_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome"

python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${nmd_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${nmd_dir}" \
    --region_type "NMD" \
    --sample_type "NMD_inhibition"

python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${nmd_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${nmd_dir}" \
    --region_type "NMD" \
    --sample_type "HCT_control"


python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${nmd_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${nmd_dir}" \
    --region_type "NMD" \
    --sample_type "cancer"


ri_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/RI_genome"

# python so_categorization_coverage_pipeline.py \
#     --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/UniqueProteinORFPairs.txt" \
#     --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/so_categorization_df.csv" \
#     --ribo_coverage_path "${ri_dir}" \
#     --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/Unique_DNA_Regions_genomic_final.bed" \
#     --result_dir "${ri_dir}" \
#     --region_type "RI" \
#     --sample_type "control"

# python so_categorization_coverage_pipeline.py \
#     --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/UniqueProteinORFPairs.txt" \
#     --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/so_categorization_df.csv" \
#     --ribo_coverage_path "${ri_dir}" \
#     --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/Unique_DNA_Regions_genomic_final.bed" \
#     --result_dir "${ri_dir}" \
#     --region_type "RI" \
#     --sample_type "NMD_inhibition"


# combine genes of interesting candidates per NMD and RI transcripts
#find "${nmd_dir}/SO_coverage_categorization" -type f -name "*genes_interesting_candidates.txt" -print0 | xargs -0 cat | sort | uniq > "${nmd_dir}/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt"

# find "${ri_dir}/SO_coverage_categorization" -type f -name "*genes_interesting_candidates.txt" -print0 | xargs -0 cat | sort | uniq > "${ri_dir}/SO_coverage_categorization/combined_RI_interesting_candidate_genes.txt"

# cat "${nmd_dir}/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt" \
# "${ri_dir}/SO_coverage_categorization/combined_RI_interesting_candidate_genes.txt" \
# | sort | uniq > "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/union_RI_NMD_interesting_candidate_genes.txt"

# echo "Total number of interesting candidate genes NMD and RI"

# wc -l "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/union_RI_NMD_interesting_candidate_genes.txt"


# # get the total number of SO transcripts per RI and NMD transcripts that have ribo-cov
# python get_nr_ribocov_so_trans_across_all_samples.py \
#  --control_cat_csv "${nmd_dir}/SO_coverage_categorization/NMD_control/control_NMD_interesting_candidates.csv"\
#  --nmd_inh_cat_csv "${nmd_dir}/SO_coverage_categorization/NMD_NMD_inhibition/NMD_inhibition_NMD_interesting_candidates.csv" \
#  --region_type NMD


# # get the total number of SO transcripts per RI and NMD transcripts that have ribo-cov
# python get_nr_ribocov_so_trans_across_all_samples.py \
#  --control_cat_csv "${ri_dir}/SO_coverage_categorization/RI_control/control_RI_interesting_candidates.csv"\
#  --nmd_inh_cat_csv "${ri_dir}/SO_coverage_categorization/RI_NMD_inhibition/NMD_inhibition_RI_interesting_candidates.csv" \
#  --region_type RI


# if [[ "${rbpbase}" == true ]]; then 
#     rbpbase_script_dir="/home/ckalk/scripts/SplitORF_pipeline/riboseq-validation/downstream_analysis_validated_URs"
#     bash ../downstream_analysis_validated_URs/run_rbpbase_analysis_gene_level.sh "${rbpbase_script_dir}"> ../downstream_analysis_validated_URs/outreports_of_runs/run_rbpbase_analysis.out 2>&1
# fi

# if [[ "${pfam}" == true ]]; then 
#     bash ../downstream_analysis_validated_URs/check_validated_pfam_domains.sh > ../downstream_analysis_validated_URs/outreports_of_runs/check_validated_pfam_domains.out 2>&1
# fi