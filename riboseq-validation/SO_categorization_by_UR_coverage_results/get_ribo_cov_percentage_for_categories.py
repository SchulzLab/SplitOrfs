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

    parser.add_argument('--ribocov_cat_csv')

    return parser.parse_args()


def main(ribocov_cat_csv):
    outdir = os.path.dirname(ribocov_cat_csv)
    filename = os.path.basename(ribocov_cat_csv)
    filename = filename.removesuffix(".csv")
    filename = filename.removeprefix("so_categorization_df_")

    ribocov_cat_df = pd.read_csv(ribocov_cat_csv, header=0)

    ribocov_cat_df = ribocov_cat_df[ribocov_cat_df["ribocov_analyzable"]
                                    == "ribocov analyzable"]

    ribocov_cat_df["UR present"].value_counts().to_csv(os.path.join(
        outdir, f"{filename}_ribocov_analyzable_value_counts.csv"))


if __name__ == '__main__':
    args = parse_args()

    ribocov_cat_csv = args.ribocov_cat_csv

    # ribocov_cat_csv = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/conda_package_ribocov_test/NMD_genome/SO_coverage_categorization/NMD_HCT_control/so_categorization_df_HCT_control_NMD.csv"

    main(ribocov_cat_csv)
