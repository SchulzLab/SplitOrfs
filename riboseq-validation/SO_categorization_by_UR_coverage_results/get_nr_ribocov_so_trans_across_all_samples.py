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
                        help='interesting candidate CSV of the control samples')
    parser.add_argument('--nmd_inh_cat_csv',
                        help='interesting candidate CSV of the NMd inhibited samples')
    parser.add_argument('--region_type',
                        help='NMD or RI')

    return parser.parse_args()


def main(control_cat_csv, nmd_inh_cat_csv, region_type):
    outdir = os.path.dirname(nmd_inh_cat_csv)

    control_df = pd.read_csv(control_cat_csv, header=0)
    nmd_inh_df = pd.read_csv(nmd_inh_cat_csv, header=0)

    union_df = pd.concat([control_df, nmd_inh_df])

    assert (union_df['ribocov_analyzable'] == 'ribocov analyzable').all()
    assert (union_df['covDistinctUr'] > 0).all()

    ribo_cov_trans = union_df['OrfTransID'].unique()

    with open(os.path.join(
            outdir, f'ribocov_transcripts_union_{region_type}.txt'), 'w') as f:
        f.write("\n".join(ribo_cov_trans))

    union_df.to_csv(os.path.join(
        outdir, f'ribocov_union_{region_type}_NMD_control.txt'))

    print(f'Number of transcripts ribocov {region_type}', len(ribo_cov_trans))

    nr_so_ribocov = union_df.groupby('OrfTransID').agg(
        {'covDistinctUr': 'max'}).sum()

    print(f'Number of Split-ORFs ribocov {region_type}', nr_so_ribocov)


if __name__ == '__main__':
    args = parse_args()

    control_cat_csv = args.control_cat_csv
    nmd_inh_cat_csv = args.nmd_inh_cat_csv
    region_type = args.region_type

    # control_cat_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome/SO_coverage_categorization/NMD_control/control_NMD_interesting_candidates.csv"
    # nmd_inh_cat_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome/SO_coverage_categorization/NMD_NMD_inhibition/NMD_inhibition_NMD_interesting_candidates.csv"
    # region_type = 'NMD'

    main(control_cat_csv, nmd_inh_cat_csv, region_type)
