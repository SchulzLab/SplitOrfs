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

    parser.add_argument('--control_so_csv',
                        help='CSV of the HCT control samples Split-ORFs (including overlapping URs)')
    parser.add_argument('--nmd_inh_so_csv',
                        help='CSV of the NMD inhibited samples Split-ORFs (including overlapping URs)')
    parser.add_argument('--cancer_so_csv',
                        help='CSV of the cancer samples Split-ORFs (including overlapping URs)')
    parser.add_argument('--region_type',
                        help='NMD or RI')

    return parser.parse_args()


def main(control_so_csv, nmd_inh_so_csv, cancer_so_csv, region_type):
    outdir = os.path.dirname(nmd_inh_so_csv)

    control_df = pd.read_csv(control_so_csv, header=0)
    nmd_inh_df = pd.read_csv(nmd_inh_so_csv, header=0)
    cancer_df = pd.read_csv(cancer_so_csv, header=0)

    union_df = pd.concat([control_df, nmd_inh_df, cancer_df])

    print(f'Total number of Split-ORFs ribocov across all samples including overlapping and identical URs {region_type}', len(
        union_df['OrfID'].unique()))


if __name__ == '__main__':
    args = parse_args()

    control_so_csv = args.control_so_csv
    nmd_inh_so_csv = args.nmd_inh_so_csv
    cancer_so_csv = args.cancer_so_csv
    region_type = args.region_type

    # control_so_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/RI_genome/SO_coverage_categorization/RI_HCT_control/validated_so_df.csv"
    # nmd_inh_so_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/RI_genome/SO_coverage_categorization/RI_NMD_inhibition/validated_so_df.csv"
    # cancer_so_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/RI_genome/SO_coverage_categorization/RI_cancer/validated_so_df.csv"
    # region_type = 'RI'

    main(control_so_csv, nmd_inh_so_csv, cancer_so_csv, region_type)
