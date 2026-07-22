#!/bin/bash -l
# this script runs the downstream analysis steps for the validated regions per sample
# it creates a sunburst for all significantly ribo-cov regions of the samples
# provided to the ribo-cov module
# needs to be run from the split-orf conda environment


ribocov_dir=$1
unique_region_dir=$2
region_type=$3
script_path=$4

so_results="$unique_region_dir"/UniqueProteinORFPairs.txt
so_categorization_df="$unique_region_dir"/so_categorization_df.csv
ur_path="$unique_region_dir"/Unique_DNA_Regions_genomic_final.bed


cd "$script_path"

python so_categorization_coverage_pipeline_generalized_21_07_26.py \
    --so_results "$so_results" \
    --so_categorization_df "$so_categorization_df" \
    --ribo_coverage_path "$ribocov_dir" \
    --ur_path "$ur_path" \
    --result_dir "$ribocov_dir" \
    --region_type "$region_type" \
    --sample_type "all"

cd -

