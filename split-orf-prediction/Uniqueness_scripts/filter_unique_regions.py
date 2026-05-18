#


import sys
import pandas as pd
from pybedtools import BedTool


def main():
    def count_ORFs(unique_regions: pd.DataFrame) -> pd.DataFrame:
        '''helper function to count ORFs correctly, accounting for all ORFs found for the source
        and mapped to the respective target and for having several unique regions'''
        previous_transcript = unique_regions.iloc[0, 3]
        previous_alignment_start = unique_regions.iloc[0, 1]
        for i in range(1, len(unique_regions['NMD_transcript'])):
            current_transcript = unique_regions.iloc[i, 3]
            current_alignment_start = unique_regions.iloc[i, 1]
            # the following case happens if an ORF has several URs, then the ORF number is kept and not increased
            if current_transcript == previous_transcript and current_alignment_start == previous_alignment_start:
                # print(current_transcript)
                # print(current_alignment_start)
                unique_regions.iloc[i, 14] = unique_regions.iloc[i-1, 14]
            elif current_transcript == previous_transcript:
                unique_regions.iloc[i, 14] = unique_regions.iloc[i-1, 14] + 1
            previous_alignment_start = current_alignment_start
            previous_transcript = current_transcript
        return unique_regions

    valid_ORFs = pd.read_csv(sys.argv[1], sep=':|\t', header=None,
                             names=['target', 'protAlignPos_start', 'protAlignPos_end', 'NMD_transcript', 'ORF_nr'], engine='python')
    # UniqueProteinMatches.bed

    unique_regions = pd.read_csv(sys.argv[2], sep=':|\t',
                                 names=['gene_transcript', 'ORF_nr', 'ORF_start',
                                        'ORF_end', 'unique_start', 'unique_end'],
                                 engine='python')
    # Unique_DNA_Regions_merged.bed or Unique_Protein_Regions_merged_valid.bed

    # unique regions were filtered by all valid ORFs not the best protein match filtered ones
    # need to filter for the best matches
    unique_regions = pd.merge(
        valid_ORFs, unique_regions, on='ORF_nr', how='outer')

    unique_regions = unique_regions[unique_regions['protAlignPos_start'].notnull(
    )]

    unique_regions = unique_regions.sort_values(
        ['NMD_transcript', 'protAlignPos_start'])

    region_type = sys.argv[3]
    if region_type == 'protein':
        unique_regions['ORF_length'] = (
            unique_regions['ORF_end'] - unique_regions['ORF_start'])/3  # protein coordinates
    elif region_type == 'DNA':
        unique_regions['ORF_length'] = (
            unique_regions['ORF_end'] - unique_regions['ORF_start'])  # DNA coordinates
    else:
        raise ValueError('Please give an appropriate region type as the second argument.\
                         This can be either protein or DNA.')

    unique_regions['start_percent'] = (
        unique_regions['unique_start'])/unique_regions['ORF_length']
    unique_regions['stop_percent'] = (
        unique_regions['unique_end'])/unique_regions['ORF_length']
    unique_regions['unique_percent'] = unique_regions['stop_percent'] - \
        unique_regions['start_percent']

    # count number of ORFs per source transcript
    unique_regions['ORFcounter'] = 1
    unique_regions['nr_unique_regon'] = unique_regions.groupby(
        ['NMD_transcript', 'protAlignPos_start']).cumcount()+1

    unique_regions = count_ORFs(unique_regions)

    # indicate whether first ORF is True = 1, or if it is last ORF = -1 or intermediate = 0
    unique_regions['firstORF'] = 0
    unique_regions.loc[unique_regions['ORFcounter'] == 1, 'firstORF'] = 1
    indices = unique_regions.groupby(['NMD_transcript'])['ORFcounter'].transform(
        'max') == unique_regions['ORFcounter']
    unique_regions.loc[indices, 'firstORF'] = -1

    unique_regions['ORFcounter'] = unique_regions['ORFcounter'].astype(
        'category')
    unique_regions['firstORF'] = unique_regions['firstORF'].astype('category')

    # get dataframe of cases to inspect more closely
    weird_cases_first_orf = unique_regions[(unique_regions['firstORF'] == 1) & (
        unique_regions['start_percent'] < 0.01)]
    weird_cases_last_orf = unique_regions[(
        unique_regions['firstORF'] == -1) & (unique_regions['stop_percent'] > 0.99)]
    weird_cases_first_orf.reset_index(inplace=True)
    weird_cases_last_orf.reset_index(inplace=True)

    # get bedtool object of the unique regions
    unique_bed = BedTool(sys.argv[2])
    try:
        ORFs_first_orf = weird_cases_first_orf.apply(lambda row: ':'.join([row['ORF_nr'],
                                                                           str(int(row['unique_start'])), str(int(row['unique_end']))]), axis=1).to_list()
    except:
        ORFs_first_orf = []

    try:
        ORFs_last_orf = weird_cases_last_orf.apply(lambda row: ':'.join([row['ORF_nr'],
                                                                         str(int(row['unique_start'])), str(int(row['unique_end']))]), axis=1).to_list()
    except:
        ORFs_last_orf = []
    ORFs_to_filter = ORFs_first_orf + ORFs_last_orf

    unique_bed.filter(lambda unique_region:
                      unique_region.chrom.split(
                          ':')[1]+':'+str(unique_region.start)+':'+str(unique_region.end)
                      not in set(ORFs_to_filter)).saveas(f'{sys.argv[2][:-4]}_filtered.bed')


if __name__ == '__main__':
    main()
