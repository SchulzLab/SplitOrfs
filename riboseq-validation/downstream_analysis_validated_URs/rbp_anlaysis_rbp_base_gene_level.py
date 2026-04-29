# This script takes the validated unique regions and checks whether their genes are
# RBPs, it write a TXT file of validated RBP genes as well as the validated_so_df
# filtered for RBPs and with the RBP information usign the RBPbase database
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

    parser.add_argument('--rbp_file',
                        help='Path to RNA-binding protein CSV from RBPbase')
    parser.add_argument('--interesting_candidate_file',
                        help='TXT fiel of interesting candidate genes')
    parser.add_argument('--region_type',
                        help='NMD or RI')

    return parser.parse_args()


def main(rbp_file, interesting_candidate_file, region_type):
    outdir = os.path.join(os.path.dirname(
        interesting_candidate_file), f'rbpbase')
    os.makedirs(outdir, exist_ok=True)

    rbp_df = pd.read_excel(rbp_file, engine="openpyxl")
    # subset for genes that are RBP in human in at least one experiment!
    rbp_df = rbp_df[rbp_df['any_Hs'] == 'YES']
    rbp_df = rbp_df.rename(columns=lambda x: x.split('\n')[0])
    columns_keep = ['UNIQUE', 'ID', 'Description',
                    'hasRBD-Pfam-Hs', 'RBPDB-Hs',
                    'RNAbinding-GO-Hs', 'Pfam-Id-Hs', 'Pfam-Name-Hs',
                    'Pfam-Description-Hs', 'Gene-type-Hs', 'any_Hs',
                    'hits_Hs', 'hits_allOrganisms']
    rbp_df = rbp_df[columns_keep].copy()

    # load Ribo-seq validated Split-ORFs
    interesting_candidate_df = pd.read_csv(
        interesting_candidate_file, header=None, names=['geneID'])
    # gene IDs are unique in interesting candidates
    assert len(interesting_candidate_df['geneID'].unique()) == len(
        interesting_candidate_df['geneID'])
    # gene IDs in RBPBase are unqiue
    assert len(rbp_df['ID']) == len(rbp_df['ID'].unique())
    rbp_df_filtered = rbp_df[rbp_df['ID'].isin(
        interesting_candidate_df['geneID'])].copy()
    print(f'Number of interesting candidates {region_type}', len(
        interesting_candidate_df['geneID'].unique()))
    print(f'Number of different RBPs interesting candidates {region_type}', len(
        rbp_df_filtered.index))

    rbp_df_filtered = rbp_df_filtered.rename(
        columns={'ID': 'geneID'})

    rbp_val_df = interesting_candidate_df.merge(
        rbp_df_filtered, on='geneID', how='right')
    rbp_val_df.to_csv(os.path.join(
        outdir, region_type, 'rbp_val_df_rbpbase.csv'))

    rbp_val_df['geneID'].unique().tofile(os.path.join(
        outdir, region_type, f'interesting_candidates_ribocov_{region_type}_rbp_genes_rbpbase.txt'), sep='\n')


if __name__ == '__main__':
    args = parse_args()

    rbp_file = args.rbp_file
    interesting_candidate_file = args.interesting_candidate_file
    region_type = args.region_type

    # rbp_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/RBP_analysis/RBPBase/RBPbase_Hs_DescriptiveID.xlsx'
    # interesting_candidate_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/RI_genome/SO_coverage_categorization/combined_RI_interesting_candidate_genes.txt'
    # region_type = 'RI'

    main(rbp_file, interesting_candidate_file, region_type)
