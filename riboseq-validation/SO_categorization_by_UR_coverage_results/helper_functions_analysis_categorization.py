import re
import pandas as pd
import glob
from collections import Counter
import os
import os.path

from plotting import plot_three_category_pie


def read_riboseq_results(empirical_Ribo_findings_file, all_predicted_so_orfs, so_categorization_df):
    sample = os.path.basename(
        empirical_Ribo_findings_file).removesuffix('_unique_regions.csv')
    Ribo_results_df = pd.read_csv(
        empirical_Ribo_findings_file, header=0, index_col=0)
    Ribo_results_significant_df = Ribo_results_df[Ribo_results_df["significant"] == 1].copy(
    )
    Ribo_results_significant_df['ORF'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[1])
    Ribo_results_significant_df['OrfTransID'] = Ribo_results_significant_df['name'].apply(
        lambda x: x.split(':')[0].split('|')[1])
    all_predicted_so_orfs[sample] = 0
    all_predicted_so_orfs[sample] = all_predicted_so_orfs['OrfID'].isin(
        Ribo_results_significant_df['ORF'])
    return all_predicted_so_orfs, so_categorization_df


def validated_so_per_sample_analysis(ribo_coverage_path, all_predicted_so_orfs, so_categorization_df, sample_type):
    '''
    indicate for Split-ORFs whether they have Ribo-seq coverage in their URs per sample
    treat NMD inhibition and all other samples (control) in 2 different runs
    '''
    nr_samples = 0
    for empirical_Ribo_findings_file in glob.glob(f"{ribo_coverage_path}/*_unique_regions.csv"):
        if (sample_type == 'HEK_control' and 'HCT' and '0h' in empirical_Ribo_findings_file) or \
            (sample_type == 'control' and 'HCT' and '0h' in empirical_Ribo_findings_file) or \
            (sample_type == 'control' and 'SRR85' in empirical_Ribo_findings_file) or \
            (sample_type == 'control' and 'SRR10' in empirical_Ribo_findings_file) or \
            (sample_type == 'NMD_inhibition' and 'HCT' and '12h' in empirical_Ribo_findings_file) or \
            (sample_type == 'cancer' and 'SRR85' in empirical_Ribo_findings_file) or \
            (sample_type == 'cancer' and 'SRR10' in empirical_Ribo_findings_file) or \
            (sample_type == 'all' and 'HCT' in empirical_Ribo_findings_file) or \
            (sample_type == 'all' and 'SRR85' in empirical_Ribo_findings_file) or \
                (sample_type == 'all' and 'SRR10' in empirical_Ribo_findings_file):
            all_predicted_so_orfs, so_categorization_df = read_riboseq_results(
                empirical_Ribo_findings_file, all_predicted_so_orfs, so_categorization_df)
            nr_samples += 1
    return so_categorization_df, all_predicted_so_orfs, nr_samples


def explode_so_df(predicted_so_orfs):
    predicted_so_orfs = predicted_so_orfs[[
        'OrfTransID', 'OrfIDs', 'OrfStarts', 'geneID']].copy()
    predicted_so_orfs['OrfID'] = predicted_so_orfs.apply(
        lambda x: x['OrfIDs'].split(','), axis=1)
    predicted_so_orfs['OrfStart'] = predicted_so_orfs['OrfStarts']
    all_predicted_so_orfs = predicted_so_orfs.explode(
        ['OrfID', 'OrfStart'], ignore_index=True).copy()
    total_nr_so = len(all_predicted_so_orfs.index)
    total_nr_so_transcripts = len(
        all_predicted_so_orfs['OrfTransID'].unique())
    print('Total number of predicted Split-ORFs', total_nr_so)
    print('Total number of predicted Split-ORF transcripts',
          total_nr_so_transcripts)
    return all_predicted_so_orfs, predicted_so_orfs, total_nr_so


def subset_validated_sos_df(all_predicted_so_orfs, outdir, region_type):
    all_predicted_so_orfs['ValidationCount'] = all_predicted_so_orfs.iloc[:, range(
        6, len(all_predicted_so_orfs.columns))].sum(axis=1, numeric_only=True)
    validated_so_df = all_predicted_so_orfs[all_predicted_so_orfs['ValidationCount'] > 1].copy(
    )

    # write txt file with genes that have validated SOs
    pd.Series(validated_so_df['geneID'].unique()).to_csv(os.path.join(
        outdir, f'SO_validated_genes_{region_type}.txt'), header=False, index=False)
    nr_validated_so = len(validated_so_df.index)
    nr_validated_transcripts = len(validated_so_df['OrfTransID'].unique())
    print('Number of validated SOs:', nr_validated_so)
    print('Number of validated SO transcripts:', nr_validated_transcripts)
    return validated_so_df, nr_validated_so, nr_validated_transcripts


def subset_UR_for_expressed_genes(dna_ur_df, validated_so_df, ribo_coverage_path, outdir, region_type, sample_type):
    '''
        subset the URs for the genes that are included by TPM filter
    '''
    # get gene and transcript ID in dna_ur_df from name
    dna_ur_df['gene'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[0])
    dna_ur_df['tID'] = dna_ur_df['OrfTransID'].apply(
        lambda x: x.split('|')[1])

    genes_to_keep = filter_so_genes(ribo_coverage_path, sample_type)
    pd.Series(genes_to_keep).to_csv(os.path.join(
        outdir, f'all_considered_SO_genes_{region_type}.txt'))
    # subset and count ORFs with URs for expressed genes
    dna_ur_df_subset = dna_ur_df[(dna_ur_df['gene'].isin(genes_to_keep)) | (
        dna_ur_df['tID'].isin(validated_so_df['OrfTransID']))]
    nr_orfs_with_UR_expressed = len(dna_ur_df_subset['OrfID'].unique())

    nr_val_SO_trans = len(validated_so_df['OrfTransID'])
    nr_val_SO_trans_filtered = len(
        validated_so_df[validated_so_df['geneID'].isin(genes_to_keep)].index)
    nr_val_lost = nr_val_SO_trans - nr_val_SO_trans_filtered
    # there may not be URs validated of genes not included!
    assert nr_val_lost == 0
    return nr_orfs_with_UR_expressed, dna_ur_df_subset, genes_to_keep


def filter_so_genes(ribo_coverage_path, sample_type):
    '''
        filter genes for TPM threshold 
    '''
    genes_above_tpm_list = []
    for empirical_Ribo_findings_file in glob.glob(f"{ribo_coverage_path}/**/Unique_DNA_Regions_genomic_*_chrom_sorted.bed", recursive=True):
        if (sample_type == 'HEK_control' and 'HCT' and '0h' in empirical_Ribo_findings_file) or \
            (sample_type == 'control' and 'HCT' and '0h' in empirical_Ribo_findings_file) or \
            (sample_type == 'control' and 'SRR85' in empirical_Ribo_findings_file) or \
            (sample_type == 'control' and 'SRR10' in empirical_Ribo_findings_file) or \
            (sample_type == 'NMD_inhibition' and 'HCT' and '12h' in empirical_Ribo_findings_file) or \
            (sample_type == 'cancer' and 'SRR85' in empirical_Ribo_findings_file) or \
            (sample_type == 'cancer' and 'SRR10' in empirical_Ribo_findings_file) or \
            (sample_type == 'all' and 'HCT' in empirical_Ribo_findings_file) or \
            (sample_type == 'all' and 'SRR85' in empirical_Ribo_findings_file) or \
                (sample_type == 'all'  'SRR10' in empirical_Ribo_findings_file):
            unique_regions_tpm_tpm_filtered = pd.read_csv(
                empirical_Ribo_findings_file, header=None, index_col=None, sep='\t')
            unique_regions_tpm_tpm_filtered['gene_id'] = unique_regions_tpm_tpm_filtered.iloc[:, 3].apply(
                lambda x: x.split('|')[0])
            genes_above_tpm_list = genes_above_tpm_list + \
                unique_regions_tpm_tpm_filtered['gene_id'].to_list()
    counts_above_tpm = Counter(genes_above_tpm_list)
    genes_to_keep = []
    for gene, count in counts_above_tpm.items():
        if count > 1:
            genes_to_keep.append(gene)
    return genes_to_keep


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


def count_orfs_by_position(so_df):
    nr_first_orfs = sum(so_df['OrfPosition'] == 'first')
    nr_middle_orfs = sum(so_df['OrfPosition'] == 'middle')
    nr_last_orfs = sum(so_df['OrfPosition'] == 'last')
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def val_so_by_position(validated_so_df, nr_validated_so, outdir, region_type):
    validated_so_df = get_so_position_in_transcript(validated_so_df)
    nr_first_orfs, nr_middle_orfs, nr_last_orfs = count_orfs_by_position(
        validated_so_df)

    assert nr_first_orfs + nr_middle_orfs + nr_last_orfs == nr_validated_so

    plot_three_category_pie(nr_first_orfs,
                            nr_middle_orfs,
                            nr_last_orfs,
                            nr_validated_so,
                            ['# first ORFs', '# middle ORFs', '# last ORFs'],
                            'Pie chart of ORF positions of validated Split-ORFs',
                            outdir,
                            'positions_of_validated_SO_pie_chart',
                            region_type,
                            ['#75C1C5', '#FFC500', '#CC79A7']
                            )
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def all_URs_by_position(dna_ur_df, all_predicted_so_orfs, outdir, region_type):
    dna_ur_df = pd.merge(
        dna_ur_df, all_predicted_so_orfs.iloc[:, 1:6], on='OrfID', how='left')

    nr_DNA_URs = len(dna_ur_df.index)

    dna_ur_df = get_so_position_in_transcript(dna_ur_df)
    nr_first_orfs, nr_middle_orfs, nr_last_orfs = count_orfs_by_position(
        dna_ur_df)
    plot_three_category_pie(nr_first_orfs,
                            nr_middle_orfs,
                            nr_last_orfs,
                            nr_DNA_URs,
                            ['# first ORFs', '# middle ORFs', '# last ORFs'],
                            'Pie chart of ORF positions of sos with DNA UR',
                            outdir,
                            'positions_of_SO_with_UR_pie_chart',
                            region_type,
                            ['#75C1C5', '#FFC500', '#CC79A7']
                            )
    return nr_first_orfs, nr_middle_orfs, nr_last_orfs


def val_perc_first_middle_last_orfs_csv(nr_val_first_orfs,
                                        nr_val_middle_orfs,
                                        nr_val_last_orfs,
                                        nr_first_orfs_ur,
                                        nr_middle_orfs_ur,
                                        nr_last_orfs_ur,
                                        outdir,
                                        outname):
    perc_val_first_orfs = nr_val_first_orfs/nr_first_orfs_ur
    perc_val_middle_orfs = nr_val_middle_orfs/nr_middle_orfs_ur
    perc_val_last_orfs = nr_val_last_orfs/nr_last_orfs_ur
    perc_val_dict = {'perc_val_first_orfs': [perc_val_first_orfs],
                     'perc_val_middle_orfs': [perc_val_middle_orfs],
                     'perc_val_last_orfs': [perc_val_last_orfs]}
    perc_val_df = pd.DataFrame(perc_val_dict)
    perc_val_df.to_csv(os.path.join(outdir, outname))


def calculate_overlapping_region_percentage(start1, end1, start2, end2):
    if end1 <= start2 or end2 <= start1:
        return 0
    elif start1 < end2 and start2 < end1:
        overlap_start = max(start2, start1)
        overlap_end = min(end1, end2)
        nr_bp_overlap = overlap_end - overlap_start
        shorter_region = min(end2-start2, end1-start1)
        return nr_bp_overlap/shorter_region


def get_max_overlap_of_regions_in_df(chr_df, threshold=0.7):
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


def get_transcript_with_2_distinct_val_URs(val_dna_overlapping_ur_df):
    """this information needs to be seen, so printed or saved as a file, otherwise this function is useless"""
    val_dna_overlapping_ur_df['OrfTransID'] = val_dna_overlapping_ur_df['OrfTransID'].apply(
        lambda x: x.split(',')).copy()
    val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts()
    # how often are there two URs that do not overlap more than 50% validated
    sum(val_dna_overlapping_ur_df['OrfTransID'].explode(
    ).value_counts() > 1)
    val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts(
    )[val_dna_overlapping_ur_df['OrfTransID'].explode().value_counts() > 1]


def identify_overlapping_unique_regions(validated_so_df, dna_ur_df, outdir):
    # filter for validated ORFs
    val_dna_ur_df = dna_ur_df[dna_ur_df['OrfID'].isin(
        validated_so_df['OrfID'])].copy()

    # group several exonic URs per ORF together
    val_dna_ur_df = val_dna_ur_df.groupby('OrfID').agg({'start': 'min',
                                                        'stop': 'max',
                                                        'chr': 'first',
                                                        'ID': lambda x: ','.join(x),
                                                        'OrfTransID': 'first'}).reset_index().copy()

    # map ORF positions to IDs
    orf_id_position_map = validated_so_df.set_index('OrfID')['OrfPosition']
    val_dna_ur_df['OrfPosition'] = val_dna_ur_df['OrfID'].map(orf_id_position_map
                                                              )

    # concatenate genomic regions
    val_dna_ur_df['genomic_UR'] = val_dna_ur_df['chr'].astype(
        str) + '_' + val_dna_ur_df['start'].astype(str) + '_' + val_dna_ur_df['stop'].astype(str)

    val_dna_ur_df['OverlapPercentage'] = 0.0
    val_dna_ur_df['OrfPositionsOverlapping'] = val_dna_ur_df['OrfPosition'].apply(
        lambda x: set([x]))
    val_dna_ur_df['OrfIDsOverlapping'] = val_dna_ur_df['OrfID'].apply(
        lambda x: set([x]))

    # get completely overlapping URs
    chr_dfs = {chr: chr_df.reset_index(drop=True).copy(
    ) for chr, chr_df in val_dna_ur_df.groupby('chr')}
    chr_dfs = {chr: get_max_overlap_of_regions_in_df(
        chr_df, 0.7) for chr, chr_df in chr_dfs.items()}
    val_dna_overlapping_ur_df = pd.concat(
        chr_dfs.values()).reset_index(drop=True).copy()

    val_dna_overlapping_ur_df['ORFs_sharing_region'] = val_dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: len(x))
    val_dna_overlapping_ur_df['shared_region_type'] = val_dna_overlapping_ur_df['OrfPositionsOverlapping'].apply(
        lambda x: len(x))

    # subset for single ORF having the UR or ORFs from the same position (first, middle last)
    # val_dna_overlapping_ur_df = val_dna_overlapping_ur_df[(val_dna_overlapping_ur_df['ORFs_sharing_region'] == 1) | (
    #     val_dna_overlapping_ur_df['shared_region_type'] == 1)]

    # frozenset: order within the set does not matter!
    val_dna_overlapping_ur_df.loc[:, 'OrfIDsOverlapping'] = val_dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
        lambda x: frozenset(x))

    val_dna_overlapping_ur_df = val_dna_overlapping_ur_df.groupby('OrfIDsOverlapping').agg(
        {'genomic_UR': 'first',
            'ORFs_sharing_region': 'first',
            'OrfPosition': 'first',
            'ID': lambda x: ','.join(x),
            'OrfTransID': lambda x: ','.join(x),
            'OrfPositionsOverlapping': 'first',
            'OrfIDsOverlapping': 'first',
            'OverlapPercentage': 'max',
         }).reset_index(drop=True)

    val_dna_overlapping_ur_df['geneID'] = val_dna_overlapping_ur_df['ID'].apply(
        lambda x: x.split('|')[0])
    gene_dfs = {gene: gene_df.reset_index(drop=True).copy(
    ) for gene, gene_df in val_dna_overlapping_ur_df.groupby('geneID')}
    gene_dfs = {gene: summarize_overlapping_urs(
        gene_df) for gene, gene_df in gene_dfs.items()}
    val_dna_overlapping_ur_df = pd.concat(
        gene_dfs.values()).reset_index(drop=True).copy()

    get_transcript_with_2_distinct_val_URs(val_dna_overlapping_ur_df)

    val_dna_overlapping_ur_df['OrfPosition'].value_counts(
    ).reset_index().to_csv(os.path.join(outdir, 'distinct_URs_per_position.csv'))

    validated_so_df.to_csv(os.path.join(outdir, 'validated_so_df.csv'))

    return val_dna_overlapping_ur_df


def add_sample_info_ur_df(validated_so_df, val_dna_overlapping_ur_df, outdir, nr_samples):
    # add sample columns to UR frame
    for col in validated_so_df.columns[6:6+nr_samples]:
        val_dna_overlapping_ur_df[col] = False

    for index, row in validated_so_df.iterrows():
        # get row in ur df with the same ORF, this can only be one row!
        ORF_row = val_dna_overlapping_ur_df['OrfIDsOverlapping'].apply(
            lambda x: row['OrfID'] in x)
        # row index
        row_index = ORF_row.index[ORF_row == True]
        # same ORF may never appear twice, few ORFs are removed
        assert len(row_index) <= 1
        if len(row_index) == 1:
            for column in validated_so_df.columns[6:6+nr_samples]:
                val_dna_overlapping_ur_df.loc[row_index[0], column] = max(
                    val_dna_overlapping_ur_df.loc[row_index[0], column], row[column])

    val_dna_overlapping_ur_df.to_csv(os.path.join(
        outdir, "val_dna_overlapping_ur_df.csv"))

    return val_dna_overlapping_ur_df


def unique_regions_present(row):
    '''
    If URInFirstORF and DistinctURInSecondORF used to define what type of UR 
    transcript
    '''
    if row['URInFirstORF']:
        if row['DistinctURInSecondORF'] > 0:
            return 'UR 1. & 2. ORF'
        else:
            return 'UR 1. ORF'
    elif row['DistinctURInSecondORF'] > 0:
        return 'UR 2. ORF'
    else:
        return 'no UR'


def unique_regions_covered(row):
    '''
    covFirst and covSecond use to categorize the ribo-cov class of transcript
    '''
    if row['covFirst']:
        if row['covSecond']:
            return '1.& 2. ORF'
        else:
            return '1. ORF'
    elif row['covSecond']:
        return '2. ORF'
    else:
        return 'no UR cov'


def orf_has_coverage(row, covered_orfs):
    '''
    covFirst, Middle, Last and Second used together with OrfsWithDistinctURTrans
    only count ribo-cov for ORFs that are counted as distinctUR!
    '''
    if row['OrfID'] in covered_orfs and row['OrfID'] in row['OrfsWithDistinctURTrans']:
        if row['OrfPosition'] == 'first':
            row['covFirst'] = 1
        elif row['OrfPosition'] == 'middle':
            row['covMiddle'] = 1
            row['covSecond'] = 1
        elif row['OrfPosition'] == 'last':
            row['covLast'] = 1
            row['covSecond'] = 1
    return row


def split_orf_coverage_by_categorization(so_categorization_df, validated_so_df):
    so_categorization_df['covFirst'] = 0
    so_categorization_df['covMiddle'] = 0
    so_categorization_df['covLast'] = 0
    # coverage in middle or last is counted as coverage second
    so_categorization_df['covSecond'] = 0

    for column in ['OrfID', 'OrfStart', 'OrfPosition', 'genomic_UR', 'hasUR']:
        so_categorization_df[column] = so_categorization_df[column].apply(
            lambda x: x.split(','))

    so_categorization_df = so_categorization_df.explode(
        ['OrfID', 'OrfStart', 'OrfPosition', 'genomic_UR', 'hasUR'])

    so_categorization_df = so_categorization_df.apply(
        lambda x: orf_has_coverage(x, validated_so_df['OrfID'].tolist()), axis=1)

    # aggregate ORFs that have overlapping URs with the exact same set of ORFs
    so_categorization_df = so_categorization_df.groupby('OrfTransID').agg(
        {
            'OrfID': lambda x: ','.join(x),
            'OrfStart': lambda x: ','.join(x),
            'geneID': 'first',
            'OrfPosition': lambda x: ','.join(x),
            'genomic_UR': lambda x: ','.join(x),
            'hasUR': lambda x: ','.join(x),
            'nrOrfs': 'first',
            'nrOrfsWithUR': 'first',
            'URInFirstORF': 'first',
            'URInLastORF': 'first',
            'URInMiddleORF': 'first',
            'NrDistinctURs': 'first',
            'OrfsWithDistinctURTrans': 'first',
            'DistinctURInSecondORF': 'first',
            'covFirst': 'max',
            'covMiddle': 'max',
            'covLast': 'max',
            'covSecond': 'max',
        }).reset_index()

    so_categorization_df['covDistinctUr'] = so_categorization_df.apply(
        lambda x: min(x['NrDistinctURs'], x['covFirst'] + x['covSecond']), axis=1)

    so_categorization_df['UR present'] = so_categorization_df.apply(
        lambda x: unique_regions_present(x), axis=1)

    so_categorization_df['UR covered'] = so_categorization_df.apply(
        lambda x: unique_regions_covered(x), axis=1)

    return so_categorization_df


def get_ribocov_interesting_candidate_genes(so_categorization_df, outdir, sample_type, region_type):
    if sample_type == 'NMD_inhibition':
        so_categorization_df_covered = so_categorization_df[so_categorization_df['covSecond'] > 0]
    else:
        so_categorization_df_covered = so_categorization_df[
            so_categorization_df['UR covered'] != 'no UR cov']

    pd.Series(so_categorization_df_covered['geneID'].unique()).to_csv(os.path.join(
        outdir, f'{sample_type}_{region_type}_genes_interesting_candidates.txt'), index=False, header=False)
    so_categorization_df_covered.to_csv(os.path.join(
        outdir, f'{sample_type}_{region_type}_interesting_candidates.csv'), index=False)


def write_categorization_table(so_categorization_df, outdir, sample_type, region_type):
    covered_first_orfs = sum(so_categorization_df['UR covered'] == '1. ORF')
    covered_second_orfs = sum(so_categorization_df['UR covered'] == '2. ORF')
    covered_first_and_second_orfs = sum(
        so_categorization_df['UR covered'] == '1.& 2. ORF')
    has_ur_but_not_covered_orfs = len(so_categorization_df[(so_categorization_df['UR covered'] == 'no UR cov') & (
        so_categorization_df['ribocov_analyzable'] == 'ribocov analyzable')].index)

    category_dict = {'covered_first_orfs': [covered_first_orfs], 'covered_second_orfs': [covered_second_orfs],
                     'covered_first_and_second_orfs': [covered_first_and_second_orfs], 'has_ur_but_not_covered_orfs': [has_ur_but_not_covered_orfs]}
    pd.DataFrame.from_dict(category_dict).to_csv(os.path.join(
        outdir, f'{sample_type}_{region_type}_ribocov_category_table.csv'), index=False)
