# ----- This script calculates the number of reads mapping to mRNA ----- #
# ----- if it is below a certain threshold then the sample will be filtered out ----- #

import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    parser.add_argument('--so_categorization_path',
                        help='Path to so_categorization_path')

    parser.add_argument('--unique_prot_orf_pairs_anno_path',
                        help='Path to UniqueProteinORFPairs_annotated.txt')

    parser.add_argument('--region_type',
                        help='NMd or RI')

    return parser.parse_args()


def main(so_categorization_path, region_type, unique_prot_orf_pairs_anno_path):
    so_categorization_df = pd.read_csv(
        so_categorization_path, index_col=0)
    nr_trans_no_ur = sum(so_categorization_df['NrDistinctURs'] == 0)
    nr_trans_2_urs = sum(so_categorization_df['NrDistinctURs'] == 1)
    nr_trans_get_2_urs = sum(so_categorization_df['NrDistinctURs'] > 1)
    print(f'Nr of trans no UR for {region_type}', nr_trans_no_ur)
    print(f'Nr of trans one distinct UR for {region_type}', nr_trans_2_urs)
    print(f'Nr of trans 2 distinct URs for {region_type}', nr_trans_get_2_urs)

    nr_trans_no_ur = sum(so_categorization_df['nrOrfsWithUR'] == 0)
    nr_trans_2_urs = sum(so_categorization_df['nrOrfsWithUR'] == 1)
    nr_trans_get_2_urs = sum(so_categorization_df['nrOrfsWithUR'] > 1)
    print(f'Nr of trans no UR for {region_type}', nr_trans_no_ur)
    print(f'Nr of trans one overlap UR for {region_type}', nr_trans_2_urs)
    print(f'Nr of trans 2 overlap URs for {region_type}', nr_trans_get_2_urs)

    unique_prot_orf_pairs_df = pd.read_csv(
        unique_prot_orf_pairs_anno_path, sep='\t')
    nr_no_anno = sum(unique_prot_orf_pairs_df['NumORFAnnot'] == 0)
    nr_one_anno = sum(unique_prot_orf_pairs_df['NumORFAnnot'] == 1)
    nr_get_2_anno = sum(unique_prot_orf_pairs_df['NumORFAnnot'] > 1)
    nr_any_anno = sum(unique_prot_orf_pairs_df['NumORFAnnot'] > 0)
    print(f'Nr of trans no PFAM for {region_type}', nr_no_anno)
    print(f'Nr of trans one PFAM for {region_type}', nr_one_anno)
    print(f'Nr of trans 2 PFAMs for {region_type}', nr_get_2_anno)
    print(f'Nr of trans with any PFAM for {region_type}', nr_any_anno)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    so_categorization_path = args.so_categorization_path
    region_type = args.region_type
    unique_prot_orf_pairs_anno_path = args.unique_prot_orf_pairs_anno_path

    # so_categorization_path = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv'
    # region_type = 'NMD'
    # unique_prot_orf_pairs_anno_path= '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs_annotated.txt'

    main(so_categorization_path, region_type, unique_prot_orf_pairs_anno_path)
