# For the validated RBPs with NMD and RI get the union and intersection
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

    parser.add_argument('--ri_validated_rbp_file',
                        help='Path to validated RBPs in RI')
    parser.add_argument('--nmd_validated_rbp_file',
                        help='Path to validated RBPs in NMD')

    return parser.parse_args()


def main(ri_validated_rbp_file, nmd_validated_rbp_file):
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(ri_validated_rbp_file))), 'rbpbase_nmd_ri')
    os.makedirs(outdir, exist_ok=True)

    ri_validated_rbp_df = pd.read_csv(
        ri_validated_rbp_file, index_col=0, header=0)
    nmd_validated_rbp_df = pd.read_csv(
        nmd_validated_rbp_file, index_col=0, header=0)

    ri_intersection_df = ri_validated_rbp_df[ri_validated_rbp_df['geneID'].isin(
        nmd_validated_rbp_df['geneID'])].copy()
    nmd_intersection_df = nmd_validated_rbp_df[nmd_validated_rbp_df['geneID'].isin(
        ri_validated_rbp_df['geneID'])].copy()
    intersection_df = pd.concat(
        [ri_intersection_df, nmd_intersection_df]).sort_values(by='geneID')

    intersection_df.to_csv(os.path.join(
        outdir, 'intersection_nmd_ri_rbp_df.csv'))
    intersection_df['geneID'].unique().tofile(os.path.join(
        outdir, 'intersection_nmd_ri_rbp_genes.txt'), sep='\n')
    print('Nr of RBPs validated in NMD and RI intersection',
          len(intersection_df['geneID'].unique()))

    union_df = pd.concat([ri_validated_rbp_df, nmd_validated_rbp_df])
    union_df.to_csv(os.path.join(outdir, 'union_nmd_ri_rbp_df.csv'))
    union_df['geneID'].unique().tofile(os.path.join(
        outdir, 'union_nmd_ri_rbp_genes.txt'), sep='\n')
    print('Nr of RBPs validated in NMD and RI union',
          len(union_df['geneID'].unique()))


if __name__ == '__main__':
    args = parse_args()

    ri_validated_rbp_file = args.ri_validated_rbp_file
    nmd_validated_rbp_file = args.nmd_validated_rbp_file

    # ri_validated_rbp_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/NMD/rbp_val_df.csv"
    # nmd_validated_rbp_file = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/RI/rbp_val_df.csv"

    main(ri_validated_rbp_file, nmd_validated_rbp_file)
