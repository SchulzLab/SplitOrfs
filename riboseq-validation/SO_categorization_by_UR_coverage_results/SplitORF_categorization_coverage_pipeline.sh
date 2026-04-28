#!/bin/bash -l
# this script runs the downstream analysis steps for the validated regions per sample
# teh NMd inhibited samples are treated separately, for these validation only if 2 ORFs are Ribo-covered
# 1. get all validated regions union and statistics for NMD and RI
# 2. check for RBP validation from RBPDB: the database needs to be downloaded and paths 
# indicated for this to run


eval "$(conda shell.bash hook)"
conda activate Riboseq

rbpdb=false
pfam=false

nmd_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome"

python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${nmd_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${nmd_dir}" \
    --region_type "NMD" \
    --sample_type "control"

python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${nmd_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${nmd_dir}" \
    --region_type "NMD" \
    --sample_type "NMD_inhibition"

ri_dir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/RI_genome"

python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${ri_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${ri_dir}" \
    --region_type "RI" \
    --sample_type "control"

python so_categorization_coverage_pipeline.py \
    --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/UniqueProteinORFPairs.txt" \
    --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/so_categorization_df.csv" \
    --ribo_coverage_path "${ri_dir}" \
    --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/Unique_DNA_Regions_genomic_final.bed" \
    --result_dir "${ri_dir}" \
    --region_type "RI" \
    --sample_type "NMD_inhibition"


# combine genes of interesting candidates per NMD and RI transcripts
find "${nmd_dir}/SO_coverage_categorization" -type f -name "*genes_interesting_candidates.txt" -print0 | xargs -0 cat > "${nmd_dir}/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt"

find "${ri_dir}/SO_coverage_categorization" -type f -name "*genes_interesting_candidates.txt" -print0 | xargs -0 cat > "${ri_dir}/SO_coverage_categorization/combined_RI_interesting_candidate_genes.txt"


if [[ "${rbpdb}" == true ]]; then 
    # bash ../downstream_analysis_validated_URs/run_rbpdb_analysis.sh > ../downstream_analysis_validated_URs/outreports_of_runs/run_rbpdb_analysis.out 2>&1
    bash ../downstream_analysis_validated_URs/run_rbp2go_analysis.sh > ../downstream_analysis_validated_URs/outreports_of_runs/run_rbp2go_analysis.out 2>&1
fi

if [[ "${pfam}" == true ]]; then 
    bash ../downstream_analysis_validated_URs/check_validated_pfam_domains.sh > ../downstream_analysis_validated_URs/outreports_of_runs/check_validated_pfam_domains.out 2>&1
fi