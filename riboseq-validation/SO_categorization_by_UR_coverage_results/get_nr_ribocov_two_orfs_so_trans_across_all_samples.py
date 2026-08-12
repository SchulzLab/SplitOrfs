# This script takes ribo-cov transcripts of the NMD inhibition and the other samples
# and forms a set
# ------------------ IMPORTS ------------------ #
import os
import os.path
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description='Summarize the trans.'
    )

    # Required positional arguments

    parser.add_argument('--control_cat_csv',
                        help='interesting candidate CSV of the HCT control samples with 2 URs covered')
    parser.add_argument('--nmd_inh_cat_csv',
                        help='interesting candidate CSV of the NMD inhibited samples with 2 URs covered')
    parser.add_argument('--cancer_cat_csv',
                        help='interesting candidate CSV of the cancer samples with 2 URs covered')
    parser.add_argument('--region_type',
                        help='NMD or RI')

    return parser.parse_args()


def main(control_cat_csv, nmd_inh_cat_csv, cancer_cat_csv, region_type):
    outdir = os.path.dirname(os.path.dirname(nmd_inh_cat_csv))

    control_df = pd.read_csv(control_cat_csv, header=0)
    nmd_inh_df = pd.read_csv(nmd_inh_cat_csv, header=0)
    cancer_cat_csv = pd.read_csv(cancer_cat_csv, header=0)

    print(f'number of genes with 2 ORFs ribo-cov {region_type} cancer:', len(
        cancer_cat_csv['geneID'].unique()))
    print(f'number of genes with 2 ORFs ribo-cov {region_type} HCT control:', len(
        control_df['geneID'].unique()))
    print(f'number of genes with 2 ORFs ribo-cov {region_type} HCT NMD inhibition:', len(
        nmd_inh_df['geneID'].unique()))

    union_df = pd.concat([control_df, nmd_inh_df, cancer_cat_csv])

    assert (union_df['ribocov_analyzable'] == 'ribocov analyzable').all()
    assert (union_df['covDistinctUr'] > 0).all()

    ribo_cov_trans = union_df['OrfTransID'].unique()

    with open(os.path.join(
            outdir, f'ribocov_transcripts_two_orfs_union_{region_type}.txt'), 'w') as f:
        f.write("\n".join(ribo_cov_trans))

    ribo_cov_genes = union_df['geneID'].unique()

    with open(os.path.join(
            outdir, f'ribocov_genes_two_orfs_union_{region_type}.txt'), 'w') as f:
        f.write("\n".join(ribo_cov_genes))

    union_df.to_csv(os.path.join(
        outdir, f'ribocov_union_two_orfs_{region_type}.csv'))

    print(
        f'Number of transcripts with 2 ORFs ribo-cov {region_type}', len(ribo_cov_trans))
    print(
        f'Number of genes with 2 ORFs ribo-cov {region_type}', len(ribo_cov_genes))


if __name__ == '__main__':
    args = parse_args()

    control_cat_csv = args.control_cat_csv
    nmd_inh_cat_csv = args.nmd_inh_cat_csv
    cancer_cat_csv = args.cancer_cat_csv
    region_type = args.region_type

    # control_cat_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_HCT_control/so_categorization_two_orfs_cov_df_HCT_control_NMD.csv"
    # nmd_inh_cat_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_NMD_inhibition/so_categorization_two_orfs_cov_df_NMD_inhibition_NMD.csv"
    # cancer_cat_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_cancer/so_categorization_two_orfs_cov_df_cancer_NMD.csv"
    # region_type = 'NMD'

    main(control_cat_csv, nmd_inh_cat_csv, cancer_cat_csv, region_type)
