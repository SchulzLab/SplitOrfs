"""
Script: categorize_so_transcripts_by_urs.py

Description:
    Categorize Split-ORF transcripts by their unique regions in distinct ORFs
    as well as the position of the ORFs (first, middle, last). Create a dataframe
    of distinct unique regions as well based on Split-ORf predictions.
    Write results to CSV files for later analyses and plotting.

Main Steps:
    1. Load raw data from input files
    2. Group URs by ORF
    3. Get position of ORFs
    4. Generate df of distinct URs
    5. Generate df with information in which (first,middle,last) ORF URs are present
    6. Write results to CSV

Usage:
    python categorize_so_transcripts_by_urs.py --so_results UniqueProteinORFPairs.txt --ur_path Unique_DNA_Regions_genomic_final.bed

Arguments:
    --so_results     Path to the UniqueProteinORFPairs.txt file from the Split-ORf pipeline
    --ur_path   Path toUnique_DNA_Regions_genomic_final.bed file from the Split-ORf pipeline


Dependencies:
    - pandas

Author:
    Christina Kalk

Date:
    2026-04-09
"""


import os
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Form validated set of URs from Ribo-seq data taking into account positioning of the URs."
    )

    # Required positional arguments
    parser.add_argument("--so_results", help="Path to SO results file")
    parser.add_argument("--ur_path",
                        help="Path to Unique region genomic BED file")

    return parser.parse_args()


def load_so_results(so_results):
    predicted_so_orfs = pd.read_csv(so_results, header=0, sep='\t')
    so_transcripts = predicted_so_orfs['OrfTransID'].to_list()
    predicted_so_orfs['OrfPos'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: x.split(','))
    predicted_so_orfs['OrfStarts'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: [y.split('-')[0] for y in x])
    predicted_so_orfs['nr_SO_starts'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: len(x))
    return predicted_so_orfs, so_transcripts


def explode_so_df(predicted_so_orfs):
    predicted_so_orfs = predicted_so_orfs[[
        'OrfTransID', 'OrfIDs', 'OrfStarts', 'geneID']].copy()
    predicted_so_orfs['OrfID'] = predicted_so_orfs.apply(
        lambda x: x['OrfIDs'].split(','), axis=1)
    predicted_so_orfs['OrfStart'] = predicted_so_orfs['OrfStarts']
    all_predicted_so_orfs = predicted_so_orfs.explode(
        ['OrfID', 'OrfStart'], ignore_index=True).copy()
    return all_predicted_so_orfs, predicted_so_orfs


def load_dna_ur_df(UR_path):
    dna_ur_df = pd.read_csv(UR_path, sep='\t', header=None, names=[
                            'chr', 'start', 'stop', 'ID', 'score', 'strand'])
    dna_ur_df['OrfID'] = dna_ur_df['ID'].str.split(
        ':').apply(lambda x: x[1])
    dna_ur_df['OrfTransID'] = dna_ur_df['ID'].str.split(
        ':').apply(lambda x: x[0])

    return dna_ur_df


def classify_ur_per_orf_position(dna_ur_df, all_predicted_so_orfs):
    # group several exonic URs per ORF together, several URs per ORF are also grouped together
    dna_ur_df = dna_ur_df.groupby('OrfID').agg({'start': 'min',
                                                'stop': 'max',
                                                        'chr': 'first',
                                                        'ID': lambda x: ','.join(x),
                                                        'OrfTransID': 'first'}).reset_index().copy()

    # map ORF positions to IDs
    orf_id_position_map = all_predicted_so_orfs.set_index('OrfID')[
        'OrfPosition']
    dna_ur_df['OrfPosition'] = dna_ur_df['OrfID'].map(orf_id_position_map
                                                      )

    # concatenate genomic regions
    dna_ur_df['genomic_UR'] = dna_ur_df['chr'].astype(
        str) + '_' + dna_ur_df['start'].astype(str) + '_' + dna_ur_df['stop'].astype(str)

    dna_ur_df['OverlapPercentage'] = 0.0
    dna_ur_df['OrfPositionsOverlapping'] = dna_ur_df['OrfPosition'].apply(
        lambda x: set([x]))
    dna_ur_df['OrfIDsOverlapping'] = dna_ur_df['OrfID'].apply(
        lambda x: set([x]))
    return dna_ur_df


def calculate_overlapping_region_percentage(start1, end1, start2, end2):
    if end1 <= start2 or end2 <= start1:
        return 0
    elif start1 < end2 and start2 < end1:
        overlap_start = max(start2, start1)
        overlap_end = min(end1, end2)
        nr_bp_overlap = overlap_end - overlap_start
        shorter_region = min(end2-start2, end1-start1)
        return nr_bp_overlap/shorter_region


def get_max_overlap_of_regions_in_df(chr_df, threshold=0.8):
    for index1 in chr_df.index:
        for index2 in chr_df.index:
            if index1 != index2:
                start1 = float(chr_df.iloc[index1]['start'])
                end1 = float(chr_df.iloc[index1]['stop'])
                start2 = float(chr_df.iloc[index2]['start'])
                end2 = float(chr_df.iloc[index2]['stop'])
                overlap = calculate_overlapping_region_percentage(
                    start1, end1, start2, end2)
                if overlap >= threshold:
                    chr_df.iloc[index2]['OrfPositionsOverlapping'].add(
                        chr_df.iloc[index1]['OrfPosition'])
                    chr_df.iloc[index2]['OrfIDsOverlapping'].add(
                        chr_df.iloc[index1]['OrfID'])
                    chr_df.iloc[index1]['OrfPositionsOverlapping'].add(
                        chr_df.iloc[index2]['OrfPosition'])
                    chr_df.iloc[index1]['OrfIDsOverlapping'].add(
                        chr_df.iloc[index2]['OrfID'])
                if overlap >= float(chr_df.iloc[index1]['OverlapPercentage']):
                    chr_df.loc[index1, 'OverlapPercentage'] = overlap
                if overlap >= float(chr_df.iloc[index2]['OverlapPercentage']):
                    chr_df.loc[index2, 'OverlapPercentage'] = overlap
    return chr_df


def summarize_overlapping_urs(gene_df):
    '''
    check per gene for overlapping unique regions and only keep the first instance
    this is done because not always the same ORFs overlap for overlapping unique regions
    '''
    if len(gene_df.index) > 1:
        gene_df_return = gene_df.copy()
        # search for pairwise overlaps of the OrfIDsOverlapping
        for index1 in gene_df.index:
            # compare 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
            index2 = index1
            while index2 < len(gene_df.index) - 1:
                index2 = index2 + 1
                orf_id_overlap_1 = gene_df.iloc[index1, 6]
                orf_id_overlap_2 = gene_df.iloc[index2, 6]
                # if ORF IDs do overlap
                if len(orf_id_overlap_1 & orf_id_overlap_2) > 0:
                    # check if index still exists or is already removed
                    if index2 in gene_df_return.index:
                        # always keep index1: ensure that one region of the overlapping
                        # ones is kept in the end!
                        gene_df_return = gene_df_return.drop(index=index2)
        return gene_df_return
    else:
        # return gene df if not several URs per gene
        return gene_df


def get_so_position_in_transcript(so_df):
    # sort the ORF starts by position
    so_df['OrfStarts'] = so_df.apply(
        lambda x: sorted([int(start) for start in x['OrfStarts']]), axis=1)
    # map the ORF start to the respective position in the sorted list
    # indicate whether it is the first or a later (first, middle, last)
    so_df['OrfIndex'] = so_df.apply(
        lambda x: x['OrfStarts'].index(int(x['OrfStart'])), axis=1)
    so_df['OrfPosition'] = so_df.apply(lambda x: 'first' if x['OrfIndex'] == 0 else (
        'last' if x['OrfIndex'] == len(x['OrfStarts'])-1 else 'middle'), axis=1)
    return so_df


def identify_middle_unique_regions(row):
    if len(row['OrfPosition']) > 3:
        middle_indices = [index for index, position in enumerate(
            row['OrfPosition']) if position == 'middle']
        return sum([row['hasUR'][index] for index in middle_indices])
    elif len(row['OrfPosition']) == 3:
        return int(row['hasUR'][row['OrfPosition'].index('middle')])
    else:
        return 0


def format_categorization_df(so_categorization_df):
    def get_list_cols(so_categorization_df):
        list_cols = []
        for col in so_categorization_df.columns:
            if so_categorization_df[col].apply(lambda x: isinstance(x, list)).any():
                list_cols.append(col)
        return list_cols

    list_cols = get_list_cols(so_categorization_df)
    for col in list_cols:
        so_categorization_df[col] = so_categorization_df[col].apply(
            lambda x: ','.join(map(str, x)))
    return so_categorization_df


def so_transcript_categorization(dna_overlapping_ur_df, all_predicted_so_orfs):
    '''
    categorize Split-ORF transcripts by number of unique regions and whether these
    are in the first, middle or last ORF. Also give information about distinct unique 
    regions per transcript (multiple transcript isoforms are counted multiple times).
    Write CSV file of the results.
    '''

    all_predicted_so_orfs['hasUR'] = all_predicted_so_orfs['OrfID'].isin(
        dna_overlapping_ur_df['OrfID'])

    genomic_ur_dict = dict(
        zip(dna_overlapping_ur_df['OrfID'], dna_overlapping_ur_df['genomic_UR']))

    all_predicted_so_orfs['genomic_UR'] = all_predicted_so_orfs['OrfID'].map(
        genomic_ur_dict)

    # aggregating together conserves teh order!
    so_categorization_df = all_predicted_so_orfs.groupby('OrfTransID').agg({
        'OrfID': list,
        'OrfStart': list,
        'geneID': 'first',
        'OrfPosition': list,
        'genomic_UR': list,
        'hasUR': list}).reset_index().copy()

    so_categorization_df['nrOrfs'] = so_categorization_df['OrfID'].apply(
        lambda x: len(x))
    so_categorization_df['nrOrfsWithUR'] = so_categorization_df['hasUR'].apply(
        lambda x: sum(x))

    so_categorization_df['URInFirstORF'] = so_categorization_df.apply(
        lambda x: int(x.loc['hasUR'][x['OrfPosition'].index('first')]), axis=1)
    so_categorization_df['URInLastORF'] = so_categorization_df.apply(
        lambda x: int(x.loc['hasUR'][x['OrfPosition'].index('last')]), axis=1)

    so_categorization_df['URInMiddleORF'] = so_categorization_df.apply(
        lambda x: identify_middle_unique_regions(x), axis=1)

    # ur == ur filters out nn values
    so_categorization_df['NrDistinctURs'] = so_categorization_df['genomic_UR'].apply(
        lambda x: len(set(ur for ur in x if ur == ur)))

    return so_categorization_df, all_predicted_so_orfs


def identify_overlapping_unique_regions(all_predicted_so_orfs, dna_ur_df, outdir):
    '''
    Categorize Split-ORF transcripts by their unique regions in distinct ORFs
    as well as the position of the ORFs (first, middle, last). Create a dataframe
    of distinct unique regions as well based on Split-ORf predictions.
    Write results to CSV files for later analyses and plotting.
    '''
    dna_ur_df = classify_ur_per_orf_position(dna_ur_df, all_predicted_so_orfs)

    # get completely overlapping URs
    chr_dfs = {chr: chr_df.reset_index(drop=True).copy(
    ) for chr, chr_df in dna_ur_df.groupby('chr')}
    chr_dfs = {chr: get_max_overlap_of_regions_in_df(
        chr_df, 0.8) for chr, chr_df in chr_dfs.items()}
    dna_overlapping_ur_df = pd.concat(
        chr_dfs.values()).reset_index(drop=True).copy()

    dna_overlapping_ur_df['ORFs_sharing_region'] = dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: len(x))
    dna_overlapping_ur_df['shared_region_type'] = dna_overlapping_ur_df['OrfPositionsOverlapping'].apply(
        lambda x: len(x))

    # frozenset: order within the set does not matter!
    dna_overlapping_ur_df.loc[:, 'OrfIDsOverlapping'] = dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: frozenset(x))

    # aggregate ORFs that have overlapping URs with the exact same set of ORFs
    dna_distinct_ur_df = dna_overlapping_ur_df.groupby('OrfIDsOverlapping').agg(
        {'genomic_UR': 'first',
            'ORFs_sharing_region': 'first',
            'OrfPosition': 'first',
            'ID': lambda x: ','.join(x),
            'OrfTransID': lambda x: ','.join(x),
            'OrfPositionsOverlapping': 'first',
            'OrfIDsOverlapping': 'first',
            'OverlapPercentage': 'max',
         }).reset_index(drop=True)

    dna_distinct_ur_df['geneID'] = dna_distinct_ur_df['ID'].apply(
        lambda x: x.split('|')[0])
    gene_dfs = {gene: gene_df.reset_index(drop=True).copy(
    ) for gene, gene_df in dna_distinct_ur_df.groupby('geneID')}
    gene_dfs = {gene: summarize_overlapping_urs(
        gene_df) for gene, gene_df in gene_dfs.items()}
    dna_distinct_ur_df = pd.concat(
        gene_dfs.values()).reset_index(drop=True).copy()

    dna_distinct_ur_df['OrfPosition'].value_counts(
    ).reset_index().to_csv(os.path.join(outdir, 'distinct_URs_per_position.csv'))

    dna_distinct_ur_df.to_csv(os.path.join(
        outdir, 'dna_distinct_ur_df.csv'))

    so_categorization_df, all_predicted_so_orfs = so_transcript_categorization(
        dna_overlapping_ur_df, all_predicted_so_orfs)

    so_categorization_df = format_categorization_df(so_categorization_df)

    so_categorization_df.to_csv(os.path.join(
        outdir, 'so_categorization_df.csv'))

    all_predicted_so_orfs.to_csv(os.path.join(
        outdir, 'all_predicted_so_orfs_position.csv'))

    return dna_distinct_ur_df


def main(so_results, ur_path):
    outdir = os.path.dirname(ur_path)
    # ------------------ LOAD DNA UNIQUE REGIONS ------------------ #
    dna_ur_df = load_dna_ur_df(
        ur_path)
    predicted_so_orfs, so_transcripts = load_so_results(so_results)
    all_predicted_so_orfs, predicted_so_orfs = explode_so_df(
        predicted_so_orfs)
    all_predicted_so_orfs = get_so_position_in_transcript(
        all_predicted_so_orfs)
    dna_distinct_ur_df = identify_overlapping_unique_regions(
        all_predicted_so_orfs, dna_ur_df, outdir)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    so_results = args.so_results
    ur_path = args.ur_path

    # ur_path = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed'
    # so_results = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt'

    main(so_results, ur_path)
