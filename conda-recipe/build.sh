#!/bin/bash
# Exit on error
set -euo pipefail
: "${PREFIX:?not running under conda build}"
: "${SRC_DIR:?not running under conda build}"

echo "SRC_DIR is: $SRC_DIR"
ls -la $SRC_DIR

# Copy your scripts to $PREFIX/bin so they are available in PATH
mkdir -p $PREFIX/bin
cp $SRC_DIR/split-orf-prediction/run_splitorfs_pipeline.sh $PREFIX/bin/split-orf-prediction
cp -r $SRC_DIR/split-orf-prediction/Genomic_scripts_18_10_24 $PREFIX/bin/Genomic_scripts_18_10_24
cp -r $SRC_DIR/split-orf-prediction/SplitOrfs-master $PREFIX/bin/SplitOrfs-master
cp -r $SRC_DIR/split-orf-prediction/Uniqueness_scripts $PREFIX/bin/Uniqueness_scripts
cp $SRC_DIR/split-orf-prediction/*.Rmd $PREFIX/bin/
chmod +x $PREFIX/bin/split-orf-prediction

mkdir -p $PREFIX/bin/ribo_cov_scripts
mkdir -p $PREFIX/bin/ribo_cov_scripts/SO_categorization_by_UR_coverage_results
cp $SRC_DIR/riboseq-validation/run_Riboseq_validation_pipeline_json_input.sh $PREFIX/bin/ribo-cov
cp $SRC_DIR/riboseq-validation/*.py $PREFIX/bin/ribo_cov_scripts
cp $SRC_DIR/riboseq-validation/*.R $PREFIX/bin/ribo_cov_scripts
cp $SRC_DIR/riboseq-validation/*.Rmd $PREFIX/bin/ribo_cov_scripts
cp $SRC_DIR/riboseq-validation/*.sh $PREFIX/bin/ribo_cov_scripts
cp $SRC_DIR/riboseq-validation/SO_categorization_by_UR_coverage_results/SplitORF_categorization_coverage_pipeline_generalized_21_07_26.sh $PREFIX/bin/ribo_cov_scripts/SO_categorization_by_UR_coverage_results
cp $SRC_DIR/riboseq-validation/SO_categorization_by_UR_coverage_results/*.py $PREFIX/bin/ribo_cov_scripts/SO_categorization_by_UR_coverage_results
chmod +x $PREFIX/bin/ribo-cov