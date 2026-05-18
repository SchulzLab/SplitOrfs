# ----- This script calculates unique genomic positions to unique ----- #
# ----- transcriptomic  positions   ----- #


import os
import pandas as pd
import argparse
# set_config(transform_output='pandas')


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('--ur_dna_regions_genomic',
                        help='Path to Unique_DNA_regions_genomic.bed')

    parser.add_argument('--exon_transcript_positions',
                        help='Path to ExonCoordsWIthChr110_transcript_positions.bed')

    parser.add_argument('--out_unique_dna_gen',
                        help='Name of output Unique DNA regions genomic')

    parser.add_argument('--out_unique_dna_trans',
                        help='Name of output Unique DNA regions transcriptomic')

    parser.add_argument('--out_unique_dna_for_fasta',
                        help='Name of output Unique DNA regions genomic for FASTA creation')

    return parser.parse_args()


def main(ur_dna_regions_genomic, exon_transcript_positions,
         out_unique_dna_gen, out_unique_dna_trans, out_unique_dna_for_fasta):

    ############################################################################
    # Load data                                                                #
    ############################################################################
    # read in exon transcript positions
    exon_transcript_positions_df = pd.read_csv(exon_transcript_positions,
                                               sep='\t',
                                               header=None,
                                               names=['trans_id',
                                                      'gen_start',
                                                      'gen_end',
                                                      'trans_start',
                                                      'trans_end',
                                                      'strand',
                                                      'chr'])

    # read in genomic uniwue regions
    ur_dna_regions_genomic_df = pd.read_csv(ur_dna_regions_genomic,
                                            sep='\t',
                                            header=None,
                                            names=['chr',
                                                   'gen_start',
                                                   'gen_end',
                                                   'ur_name',
                                                   'score',
                                                   'strand'])

    ############################################################################
    # Get necessary columns                                                    #
    ############################################################################
    # define new columns
    ur_dna_regions_genomic_df['trans_id'] = ur_dna_regions_genomic_df['ur_name'].apply(
        lambda x: x.split(':')[0])

    ur_dna_regions_genomic_df['orf'] = ur_dna_regions_genomic_df['ur_name'].apply(
        lambda x: x.split(':')[1])

    ur_dna_regions_genomic_df['orf_start'] = ur_dna_regions_genomic_df['ur_name'].apply(
        lambda x: x.split(':')[2]).astype(int)

    ur_dna_regions_genomic_df['orf_end'] = ur_dna_regions_genomic_df['ur_name'].apply(
        lambda x: x.split(':')[3]).astype(int)

    ur_dna_regions_genomic_df['prev_ur_start'] = ur_dna_regions_genomic_df['ur_name'].apply(
        lambda x: x.split(':')[4]).astype(int)

    ur_dna_regions_genomic_df['prev_ur_end'] = ur_dna_regions_genomic_df['ur_name'].apply(
        lambda x: x.split(':')[5]).astype(int)

    ur_dna_regions_genomic_df['ur_trans_start'] = 0

    ur_dna_regions_genomic_df['ur_trans_end'] = 0

    exon_transcript_positions_df = exon_transcript_positions_df[exon_transcript_positions_df['trans_id'].isin(
        ur_dna_regions_genomic_df['trans_id'])]

    ############################################################################
    # Transform genomic to transcriptomic positions                            #
    ############################################################################
    # calculate the transcriptomic coordinates from the genomic coordinates of the unique regions
    for ur in ur_dna_regions_genomic_df.itertuples():
        # 7th column in 'trans_id'
        trans_id = ur[7]
        exon_transcript_positions_df_ur = exon_transcript_positions_df[
            exon_transcript_positions_df['trans_id'] == trans_id]
        for row in exon_transcript_positions_df_ur.itertuples():
            trans_chr = row[7]
            gen_start_trans = row[2]
            gen_end_trans = row[3]
            trans_start = row[4]

            index = ur[0]
            ur_chr = ur[1]
            ur_gen_start = ur[2]
            ur_gen_end = ur[3]
            # 6th column is strand
            if ur[6] == '+':
                if ur_gen_start >= gen_start_trans and ur_gen_start <= gen_end_trans:
                    assert (trans_chr == ur_chr)
                    assert (ur_gen_end <= gen_end_trans)
                    start_offset = ur_gen_start - gen_start_trans
                    end_offset = ur_gen_end - gen_start_trans
                    trans_start_ur = trans_start + start_offset
                    trans_end_ur = trans_start + end_offset
                    ur_dna_regions_genomic_df.loc[index,
                                                  'ur_trans_start'] = trans_start_ur
                    ur_dna_regions_genomic_df.loc[index,
                                                  'ur_trans_end'] = trans_end_ur

            else:
                if ur_gen_start >= gen_start_trans and ur_gen_start <= gen_end_trans:
                    assert (trans_chr == ur_chr)
                    assert (ur_gen_end <= gen_end_trans)
                    end_offset = gen_end_trans - ur_gen_start
                    start_offset = gen_end_trans - ur_gen_end
                    trans_start_ur = trans_start + start_offset
                    trans_end_ur = trans_start + end_offset
                    ur_dna_regions_genomic_df.loc[index,
                                                  'ur_trans_start'] = trans_start_ur
                    ur_dna_regions_genomic_df.loc[index,
                                                  'ur_trans_end'] = trans_end_ur

    # group regions that span exons: in transcriptomic coordinates as one region
    ur_dna_regions_genomic_df_grouped = ur_dna_regions_genomic_df.groupby('ur_name').agg(
        {'trans_id': 'first', 'orf': 'first', 'orf_start': 'first',
         'orf_end': 'first', 'prev_ur_start': 'first', 'prev_ur_end': 'first',
         'ur_trans_start': 'min', 'ur_trans_end': 'max'})

    # create new orf names
    ur_dna_regions_genomic_df_grouped = ur_dna_regions_genomic_df_grouped.reset_index()
    ur_dna_regions_genomic_df_grouped['ur_name_old'] = ur_dna_regions_genomic_df_grouped['ur_name']
    # get a dict with the old and the new names
    ur_dna_regions_genomic_df_grouped['ur_name'] = ur_dna_regions_genomic_df_grouped['trans_id'] \
        + ':' + ur_dna_regions_genomic_df_grouped['orf'] + ':' + \
        ur_dna_regions_genomic_df_grouped['orf_start'].astype(
        str) + ':' + ur_dna_regions_genomic_df_grouped['orf_end'].astype(str) \
        + ':' + ur_dna_regions_genomic_df_grouped['ur_trans_start'].astype(str) \
        + ':' + ur_dna_regions_genomic_df_grouped['ur_trans_end'].astype(str)

    ur_dna_regions_genomic_df_grouped['orf_name'] = ur_dna_regions_genomic_df_grouped['orf'] \
        + ':' + ur_dna_regions_genomic_df_grouped['orf_start'].astype(
        str) + ':' + ur_dna_regions_genomic_df_grouped['orf_end'].astype(str)

    ur_name_dict = ur_dna_regions_genomic_df_grouped.set_index('ur_name_old')[
        'ur_name'].to_dict()

    ############################################################################
    # Write results                                                            #
    ############################################################################
    # write genomic coordinates with correct orf names
    ur_dna_regions_genomic_df['ur_name_new'] = ur_dna_regions_genomic_df['ur_name'].map(
        ur_name_dict)
    ur_dna_regions_genomic_df[['chr', 'gen_start', 'gen_end', 'ur_name_new', 'score', 'strand']].to_csv(
        out_unique_dna_gen, sep='\t', header=False, index=False)

    # write CDS corrected transcriptomic coordinates
    ur_dna_regions_genomic_df_grouped[['trans_id', 'ur_trans_start', 'ur_trans_end', 'orf_name']].to_csv(
        out_unique_dna_trans, sep='\t', header=False, index=False)

    ur_dna_regions_genomic_df_grouped['ur_trans_start'] = ur_dna_regions_genomic_df_grouped['ur_trans_start'] - \
        ur_dna_regions_genomic_df_grouped['orf_start']
    ur_dna_regions_genomic_df_grouped['ur_trans_end'] = ur_dna_regions_genomic_df_grouped['ur_trans_end'] - \
        ur_dna_regions_genomic_df_grouped['orf_start']
    ur_dna_regions_genomic_df_grouped['orf_id'] = ur_dna_regions_genomic_df_grouped['trans_id'] + \
        ':' + ur_dna_regions_genomic_df_grouped['orf_name']

    ur_dna_regions_genomic_df_grouped[['orf_id', 'ur_trans_start', 'ur_trans_end']].to_csv(
        out_unique_dna_for_fasta, sep='\t', header=False, index=False)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    ur_dna_regions_genomic = args.ur_dna_regions_genomic
    exon_transcript_positions = args.exon_transcript_positions
    out_unique_dna_gen = args.out_unique_dna_gen
    out_unique_dna_trans = args.out_unique_dna_trans
    out_unique_dna_for_fasta = args.out_unique_dna_for_fasta

    main(ur_dna_regions_genomic, exon_transcript_positions,
         out_unique_dna_gen, out_unique_dna_trans, out_unique_dna_for_fasta)
