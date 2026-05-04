#!/bin/bash -l

eval "$(conda shell.bash hook)"
conda activate Riboseq

script_dir=$1
# New run for the SO categorization approach:

################################ NMD ###########################################
python "${script_dir}"/rbp_anlaysis_rbp_base_gene_level.py \
    --rbp_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/RBPBase/RBPbase_Hs_DescriptiveID.xlsx' \
    --interesting_candidate_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome/SO_coverage_categorization/combined_NMD_interesting_candidate_genes.txt' \
    --region_type 'NMD'


################################ RI ###########################################
python "${script_dir}"/rbp_anlaysis_rbp_base_gene_level.py \
    --rbp_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/RBPBase/RBPbase_Hs_DescriptiveID.xlsx' \
    --interesting_candidate_file '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/RI_genome/SO_coverage_categorization/combined_RI_interesting_candidate_genes.txt' \
    --region_type 'RI'


python "${script_dir}"/ri_nmd_union_intersection_rbps.py \
    --nmd_validated_rbp_file "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome/SO_coverage_categorization/rbpbase/interesting_candidates_ribocov_NMD_rbp_genes_rbpbase.txt" \
    --ri_validated_rbp_file "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/RI_genome/SO_coverage_categorization/rbpbase/interesting_candidates_ribocov_RI_rbp_genes_rbpbase.txt"