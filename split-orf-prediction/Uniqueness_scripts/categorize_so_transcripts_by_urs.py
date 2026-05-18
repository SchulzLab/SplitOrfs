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
import copy


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
    predicted_so_orfs = copy.deepcopy(predicted_so_orfs[[
        'OrfTransID', 'OrfIDs', 'OrfStarts', 'geneID']])
    predicted_so_orfs['OrfID'] = predicted_so_orfs.apply(
        lambda x: x['OrfIDs'].split(','), axis=1)
    predicted_so_orfs['OrfStart'] = predicted_so_orfs['OrfStarts']
    all_predicted_so_orfs = copy.deepcopy(predicted_so_orfs.explode(
        ['OrfID', 'OrfStart'], ignore_index=True))
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
    dna_ur_df = copy.deepcopy(dna_ur_df.groupby('OrfID').agg({'start': 'min',
                                                              'stop': 'max',
                                                              'chr': 'first',
                                                              'ID': lambda x: ','.join(x),
                                                              'OrfTransID': 'first'}).reset_index())

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


# def get_max_overlap_of_regions_in_df(chr_df, threshold=0.2):
#     starts = chr_df['start'].to_numpy()
#     ends = chr_df['stop'].to_numpy()

#     for i in range(len(starts)):
#         for j in range(i + 1, len(starts)):  # avoid duplicate + self-comparison
#             overlap = calculate_overlapping_region_percentage(
#                 starts[i], ends[i], starts[j], ends[j]
#             )
#             if overlap >= threshold:
#                 chr_df.iloc[j]['OrfPositionsOverlapping'].add(
#                     chr_df.iloc[i]['OrfPosition'])
#                 chr_df.iloc[j]['OrfIDsOverlapping'].add(
#                     chr_df.iloc[i]['OrfID'])
#                 chr_df.iloc[i]['OrfPositionsOverlapping'].add(
#                     chr_df.iloc[j]['OrfPosition'])
#                 chr_df.iloc[i]['OrfIDsOverlapping'].add(
#                     chr_df.iloc[j]['OrfID'])
#                 if overlap >= float(chr_df.iloc[i]['OverlapPercentage']):
#                     chr_df.loc[i, 'OverlapPercentage'] = overlap
#                 if overlap >= float(chr_df.iloc[j]['OverlapPercentage']):
#                     chr_df.loc[j, 'OverlapPercentage'] = overlap
#     return chr_df


def get_max_overlap_of_regions_in_df(chr_df, threshold=0.2):
    starts = chr_df['start'].to_numpy()
    ends = chr_df['stop'].to_numpy()

    for i in range(len(starts)):
        for j in range(i + 1, len(starts)):  # avoid duplicate + self-comparison
            overlap = calculate_overlapping_region_percentage(
                starts[i], ends[i], starts[j], ends[j]
            )
            if overlap >= threshold:
                chr_df.at[j, 'OrfPositionsOverlapping'] = set(chr_df.loc[j,
                                                                         'OrfPositionsOverlapping']) | {chr_df.loc[i, 'OrfPosition']}
                chr_df.at[j, 'OrfIDsOverlapping'] = set(
                    chr_df.loc[j, 'OrfIDsOverlapping']) | {chr_df.loc[i, 'OrfID']}
                chr_df.at[i, 'OrfPositionsOverlapping'] = set(chr_df.loc[i, 'OrfPositionsOverlapping']) | {
                    chr_df.loc[j, 'OrfPosition']}
                chr_df.at[i, 'OrfIDsOverlapping'] = set(
                    chr_df.loc[i, 'OrfIDsOverlapping']) | {chr_df.loc[j, 'OrfID']}
                if overlap >= float(chr_df.loc[i, 'OverlapPercentage']):
                    chr_df.loc[i, 'OverlapPercentage'] = overlap
                if overlap >= float(chr_df.loc[j, 'OverlapPercentage']):
                    chr_df.loc[j, 'OverlapPercentage'] = overlap
    return chr_df


def summarize_overlapping_urs(gene_df, col_index):
    '''
    check per gene for overlapping unique regions and only keep the first instance
    this is done because not always the same ORFs overlap for overlapping unique regions
    '''
    if len(gene_df.index) > 1:
        gene_df_return = copy.deepcopy(gene_df)
        # search for pairwise overlaps of the OrfIDsOverlapping
        for index1 in gene_df.index:
            # compare 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
            index2 = index1
            while index2 < len(gene_df.index) - 1:
                index2 = index2 + 1
                orf_id_overlap_1 = gene_df.iloc[index1, col_index]
                orf_id_overlap_2 = gene_df.iloc[index2, col_index]
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
            if so_categorization_df[col].apply(lambda x: isinstance(x, list)).all():
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
    so_categorization_df = copy.deepcopy(all_predicted_so_orfs.groupby('OrfTransID').agg({
        'OrfID': list,
        'OrfStart': list,
        'geneID': 'first',
        'OrfPosition': list,
        'genomic_UR': list,
        'hasUR': list}).reset_index())

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
    # so_categorization_df['NrDistinctURs'] = so_categorization_df['genomic_UR'].apply(
    #     lambda x: len(set(ur for ur in x if ur == ur)))

    return so_categorization_df, all_predicted_so_orfs


def overlapping_ur_df_by_id(dna_ur_df, outdir, id, agg_col_index):
    '''
    group URs by gene, transcript or chromosomes and summarize if they overlap more than
    the indicated threshold: 0.2
    '''
    # get completely overlapping URs
    gene_dfs = {gene: copy.deepcopy(gene_df.reset_index(
        drop=True)) for gene, gene_df in dna_ur_df.groupby(id)}
    gene_dfs = {gene: get_max_overlap_of_regions_in_df(
        gene_df, 0.2) for gene, gene_df in gene_dfs.items()}
    dna_overlapping_ur_df = copy.deepcopy(pd.concat(
        gene_dfs.values()).reset_index(drop=True))

    dna_overlapping_ur_df['ORFs_sharing_region'] = dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: len(x))
    dna_overlapping_ur_df['shared_region_type'] = dna_overlapping_ur_df['OrfPositionsOverlapping'].apply(
        lambda x: len(x))
    # frozenset: order within the set does not matter!
    dna_overlapping_ur_df.loc[:, 'OrfIDsOverlapping'] = dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: frozenset(x))

    if id == 'OrfTransID':
        # aggregate ORFs that have overlapping URs with the exact same set of ORFs
        dna_distinct_ur_df = copy.deepcopy(dna_overlapping_ur_df.groupby('OrfIDsOverlapping').agg(
            {'genomic_UR': 'first',
             'ORFs_sharing_region': 'first',
             'OrfPosition': 'first',
             'ID': lambda x: ','.join(x),
             'OrfTransID': 'first',
             'OrfPositionsOverlapping': 'first',
             'OrfIDsOverlapping': 'first',
             'OverlapPercentage': 'max',
             'geneID': 'first',
             }).reset_index(drop=True))
    else:
        # aggregate ORFs that have overlapping URs with the exact same set of ORFs
        dna_distinct_ur_df = copy.deepcopy(dna_overlapping_ur_df.groupby('OrfIDsOverlapping').agg(
            {'genomic_UR': 'first',
             'ORFs_sharing_region': 'first',
             'OrfPosition': 'first',
             'ID': lambda x: ','.join(x),
             'OrfTransID': lambda x: ','.join(x),
             'OrfPositionsOverlapping': 'first',
             'OrfIDsOverlapping': 'first',
             'OverlapPercentage': 'max',
             'geneID': 'first',
             }).reset_index(drop=True))

    gene_dfs = {gene: gene_df.reset_index(drop=True).copy(
    ) for gene, gene_df in dna_distinct_ur_df.groupby(id)}
    gene_dfs = {gene: summarize_overlapping_urs(
        gene_df, agg_col_index) for gene, gene_df in gene_dfs.items()}
    dna_distinct_ur_df = copy.deepcopy(pd.concat(
        gene_dfs.values()).reset_index(drop=True))

    dna_distinct_ur_df['OrfPosition'].value_counts(
    ).reset_index().to_csv(os.path.join(outdir, f'distinct_URs_per_position_{id}.csv'))
    dna_distinct_ur_df.to_csv(os.path.join(
        outdir, f'dna_distinct_ur_df_{id}.csv'))

    return dna_distinct_ur_df, dna_overlapping_ur_df


def identify_overlapping_unique_regions(all_predicted_so_orfs, dna_ur_df, outdir):
    '''
    Categorize Split-ORF transcripts by their unique regions in distinct ORFs
    as well as the position of the ORFs (first, middle, last). Create a dataframe
    of distinct unique regions as well based on Split-ORf predictions.
    Write results to CSV files for later analyses and plotting.
    '''
    dna_ur_df = classify_ur_per_orf_position(dna_ur_df, all_predicted_so_orfs)

    # get completely overlapping URs on gene level
    dna_ur_df['geneID'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[0])

    dna_ur_df_unaltered = copy.deepcopy(dna_ur_df)

    dna_distinct_ur_df_trans, dna_overlapping_ur_df_trans = overlapping_ur_df_by_id(
        dna_ur_df_unaltered, outdir, 'OrfTransID', 6)
    dna_distinct_ur_df, dna_overlapping_ur_df = overlapping_ur_df_by_id(
        dna_ur_df, outdir, 'geneID', 6)

    so_categorization_df, all_predicted_so_orfs = so_transcript_categorization(
        dna_overlapping_ur_df_trans, all_predicted_so_orfs)

    so_categorization_df = format_categorization_df(so_categorization_df)

    # check that number of unique regions correspond!
    # if there is no unique region this is np.nan and np.nan != np.nan
    assert (so_categorization_df['genomic_UR'].apply(lambda x:
                                                     len([ur for ur in x.split(',') if ur != 'nan'])) == so_categorization_df['nrOrfsWithUR']).all()

    all_predicted_so_orfs.to_csv(os.path.join(
        outdir, 'all_predicted_so_orfs_position.csv'))

    return dna_distinct_ur_df, so_categorization_df, dna_distinct_ur_df_trans, dna_overlapping_ur_df_trans


def get_orfs_with_ur(row):
    has_ur_list = row['hasUR'].split(',')
    has_ur_list = [eval(has_ur) for has_ur in has_ur_list]
    orf_id_list = row['OrfID'].split(',')
    orfs_with_ur_list = [v for v, m in zip(orf_id_list, has_ur_list) if m]
    return orfs_with_ur_list


def nr_of_non_overlapping_urs(orfs_with_ur_list, dna_overlapping_ur_df, trans_id):
    '''
    get number of distinct UR ORFs
    '''
    # subset the unique df for the same gene
    orf_set = frozenset(orfs_with_ur_list)
    dna_overlapping_ur_df_sub = dna_overlapping_ur_df[dna_overlapping_ur_df['OrfTransID'] == trans_id]
    overlapping_ur_sets = dna_overlapping_ur_df_sub.apply(
        lambda x: x['OrfIDsOverlapping'].intersection(orf_set), axis=1)
    # reutrn all ORFs with UR, unless there are overlapping ORFs:
    # then at least some of the sets are > 1: subtract the number of overlapping ORFs
    # this works because each ORF is only listed once
    return len(orf_set) - sum(overlapping_ur_sets.apply(lambda x: len(x) - 1 if len(x) > 0 else 0))


def overlapping_orf_ids_within_trans(orfs_with_ur_list, dna_overlapping_ur_df, trans_id):
    '''
    get set of overlapping ORF IDs within transcript as frozen set
    '''
    # subset the unique df for the same gene
    orf_set = frozenset(orfs_with_ur_list)
    if len(orf_set) > 0:
        dna_overlapping_ur_df_sub = dna_overlapping_ur_df[
            dna_overlapping_ur_df['OrfTransID'] == trans_id]
        overlapping_ur_sets = dna_overlapping_ur_df_sub.apply(
            lambda x: x['OrfIDsOverlapping'].intersection(orf_set), axis=1).values
        overlapping_ur_sets = [
            overlap_set for overlap_set in overlapping_ur_sets if len(overlap_set) > 1]
        if len(overlapping_ur_sets) > 0:
            return overlapping_ur_sets
        else:
            return None
    else:
        return None


def check_for_overlapping_orfs_within_trans(row, dna_overlapping_ur_df):
    orfs_with_ur_list = get_orfs_with_ur(row)
    trans_id = row['geneID'] + '|' + row['OrfTransID']
    nr_distinct_urs = nr_of_non_overlapping_urs(
        orfs_with_ur_list, dna_overlapping_ur_df, trans_id)
    return nr_distinct_urs


def get_overlapping_ur_orfs_within_trans(row, dna_overlapping_ur_df):
    '''
    apply function to get the string of non-overlapping ORFs to consider 
    for ribo-seq coverage
    '''
    orfs_with_ur_list = get_orfs_with_ur(row)
    trans_id = row['geneID'] + '|' + row['OrfTransID']
    overlapping_ur_orfs_within_trans = overlapping_orf_ids_within_trans(
        orfs_with_ur_list, dna_overlapping_ur_df, trans_id)
    return overlapping_ur_orfs_within_trans


def orf_id_overlapping_first(row):
    '''
    Which ORF ID(s) are overlapping with the most 5' (first) ORF?
    '''
    if row['UROverlapWithinTrans'] == None:
        return None
    elif row['URInFirstORF'] == 0:
        return None
    else:
        # there is only one overlapping set where the first ORF is in!
        ov_set_first = [
            ov_set for ov_set in row['UROverlapWithinTrans'] if row['IDfirstORF'] in ov_set]
        if len(ov_set_first) > 0:
            ov_set_first = ov_set_first[0]
            return [orf for orf in ov_set_first if orf != row['IDfirstORF']]
        else:
            return None


def assign_has_distinct_ur(row):
    '''
    boolean comma-sep string of whether the repsective ORF has a distinct UR or not
    '''
    if row['nrOrfsWithUR'] == 0:
        return row['hasUR']
    elif row['nrOrfsWithUR'] == 1:
        return row['hasUR']
    else:
        overlapping_orf_sets = row['UROverlapWithinTrans']
        if overlapping_orf_sets != None:
            orf_dict = {}
            more_3_prime_orfs = []
            for overlapping_set in overlapping_orf_sets:
                for orf in overlapping_set:
                    orf_index = row['OrfID'].index(orf)
                    orf_start = row['OrfStart'][orf_index]
                    orf_dict[orf] = int(orf_start)

                most_5_prime = min(orf_dict, key=orf_dict.get)
                more_3_prime_orfs.extend(
                    [orf for orf in orf_dict.keys() if orf != most_5_prime])

            orfs_no_distinct_ur = [index for index, orf in enumerate(
                row['OrfID']) if orf in more_3_prime_orfs]
            has_distinct_ur_list = [
                x if i not in orfs_no_distinct_ur else False for i, x in enumerate(row['hasDistinctUR'])]
            return ','.join([str(has_ur) for has_ur in has_distinct_ur_list])

        else:
            return row['hasUR']


def orfs_for_which_ur_counts(row):
    '''
    assign each UR to one ORF, if overlap then the most
    5' ORF is assigned the UR
    only report these ORFs with distinct UR as a comma separated string
    '''
    if row['nrOrfsWithUR'] == 0:
        return None
    elif row['nrOrfsWithUR'] == 1:
        return row['OrfID'][row['hasDistinctUR'].index(True)]
    else:
        overlapping_orf_sets = row['UROverlapWithinTrans']
        if overlapping_orf_sets != None:
            orf_dict = {}
            # most_5_prime_orfs = []
            more_3_prime_orfs = []
            for overlapping_set in overlapping_orf_sets:
                for orf in overlapping_set:
                    orf_index = row['OrfID'].index(orf)
                    orf_start = row['OrfStart'][orf_index]
                    orf_dict[orf] = int(orf_start)

                most_5_prime = min(orf_dict, key=orf_dict.get)
                # most_5_prime_orfs.append(orf_dict[most_5_prime])
                more_3_prime_orfs.extend(
                    [orf for orf in orf_dict.keys() if orf != most_5_prime])

            all_orfs_with_ur_list = [v for v, m in zip(
                row['OrfID'], row['hasDistinctUR']) if m]

            orfs_with_distinct_ur = [
                orf for orf in all_orfs_with_ur_list if orf not in more_3_prime_orfs]

            orf_with_ur_string = ','.join(orfs_with_distinct_ur)

        else:
            orf_with_ur_string = ','.join(v for v, m in zip(
                row['OrfID'], row['hasDistinctUR']) if m)
        return orf_with_ur_string


def get_overlapping_info(dna_distinct_ur_df, so_categorization_df, dna_distinct_ur_df_trans):
    '''
    Add information to categorization_df of non-overlapping ORF URs to consider for 
    ribo-seq coverage in different formats
    '''
    so_categorization_df['UROverlapWithinTrans'] = ''
    so_categorization_df['UROverlapGeneral'] = ''
    so_categorization_df['NrDistinctURs'] = 0

    # subset for overlapping URs only
    dna_overlapping_ur_df = dna_distinct_ur_df_trans[dna_distinct_ur_df_trans['OrfIDsOverlapping'].apply(
        lambda x: len(x) > 1)]

    so_categorization_df['NrDistinctURs'] = so_categorization_df.apply(
        lambda x: check_for_overlapping_orfs_within_trans(x, dna_overlapping_ur_df), axis=1)
    so_categorization_df['UROverlapWithinTrans'] = so_categorization_df.apply(
        lambda x: get_overlapping_ur_orfs_within_trans(x, dna_overlapping_ur_df), axis=1)
    so_categorization_df['NrOverlapURsWithinTrans'] = so_categorization_df['nrOrfsWithUR'] - \
        so_categorization_df['NrDistinctURs']

    # check that overlapping URs within transcript numbers correspond
    assert sum(so_categorization_df['UROverlapWithinTrans'].apply(
        lambda x: x != None)) == sum(so_categorization_df['NrOverlapURsWithinTrans'] > 0)

    # assign this with the number of distinct URs in Second ORFs
    so_categorization_df['DistinctURInSecondORF'] = so_categorization_df['NrDistinctURs'] - \
        so_categorization_df['URInFirstORF']

    # keep the UR only for one ORF if overlap: always for the more 5' one!
    # have it as a True,False list etc
    # prepare the different cols accordingly
    so_categorization_df['hasDistinctUR'] = so_categorization_df['hasUR'].apply(
        lambda x: [eval(ur_indicator) for ur_indicator in x.split(',')])
    so_categorization_df['OrfPosition'] = so_categorization_df['OrfPosition'].apply(
        lambda x: x.split(','))
    so_categorization_df['OrfID'] = so_categorization_df['OrfID'].apply(
        lambda x: x.split(','))
    so_categorization_df['OrfStart'] = so_categorization_df['OrfStart'].apply(
        lambda x: x.split(','))
    so_categorization_df['IDfirstORF'] = so_categorization_df.apply(
        lambda x: x.loc['OrfID'][x['OrfPosition'].index('first')], axis=1)

    so_categorization_df['IDOverlapfirstORF'] = so_categorization_df.apply(
        lambda row: orf_id_overlapping_first(row), axis=1)

    so_categorization_df['OrfsWithDistinctURTrans'] = so_categorization_df.apply(
        lambda row: orfs_for_which_ur_counts(row), axis=1)
    so_categorization_df['hasDistinctUR'] = so_categorization_df.apply(
        lambda row: assign_has_distinct_ur(row), axis=1)

    # NR distinct URs == Nr ORFs with UR minus what overlaps and is not counted
    # here there is only one set of ORFs that overlaps per trnascript
    # this assertion might need to be adapted for other datasets
    print('UR test', (so_categorization_df['NrDistinctURs'] == so_categorization_df.apply(
        lambda x: x['nrOrfsWithUR'] if x['UROverlapWithinTrans'] == None else x['nrOrfsWithUR'] - len(x['UROverlapWithinTrans'][0]) + 1, axis=1)).all())

    # disntinct UR numbers need to correspond!
    assert (so_categorization_df['hasDistinctUR'].apply(lambda x: sum([eval(ur_bool) for ur_bool in x.split(',')])) ==
            so_categorization_df['NrDistinctURs']).all()

    # ORFs with distinct UR need to be same number as number of distinct URs, as each ORf is only considered for one UR!
    assert (so_categorization_df['OrfsWithDistinctURTrans'].apply(
        lambda x: len(x.split(',')) if isinstance(x, str) else 0) == so_categorization_df['NrDistinctURs']).all()

    so_categorization_df = format_categorization_df(so_categorization_df)

    return so_categorization_df


def write_genes_with_ur(so_categorization_df, outdir):
    '''
    write TXT file of gene names which can be found with ribocov method, ie genes that 
    have a unqiue region
    '''
    so_categorization_df_has_ur = so_categorization_df[so_categorization_df['nrOrfsWithUR'] > 0]
    pd.Series(so_categorization_df_has_ur['geneID'].unique()).to_csv(
        os.path.join(outdir, 'genes_with_unique_region.txt'), index=False)

# def write_categorization_table(so_categorization_df, outdir):
#     '''
#     write CSV file which indicates how many genes, transcripts, ORFs are there with
#     categorization by unique regions
#     '''


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

    dna_distinct_ur_df, so_categorization_df, dna_distinct_ur_df_trans, \
        dna_overlapping_ur_df_trans = identify_overlapping_unique_regions(
            all_predicted_so_orfs, dna_ur_df, outdir)
    so_categorization_df = get_overlapping_info(
        dna_distinct_ur_df, so_categorization_df, dna_distinct_ur_df_trans)
    write_genes_with_ur(so_categorization_df, outdir)
    # write_categorization_table(so_categorization_df, outdir)
    so_categorization_df.to_csv(os.path.join(
        outdir, 'so_categorization_df.csv'))


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    so_results = args.so_results
    ur_path = args.ur_path

    main(so_results, ur_path)
