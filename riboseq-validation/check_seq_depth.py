# ----- This script calculates the number of reads mapping to mRNA ----- #
# ----- if it is below a certain threshold then the sample will be filtered out ----- #

import pandas as pd
import argparse
# set_config(transform_output='pandas')


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('htseq_counts',
                        help='Path to HTSeq-counts')

    parser.add_argument('filter_treshold',
                        help='Threshold of reads mapping to genes to filter samples')

    return parser.parse_args()


def main(htseq_counts, filter_treshold):
    htseq_counts_df = pd.read_csv(
        htseq_counts, sep='\t', header=None, names=['Gene_name', 'count'])
    # filter out for the non-gene columns
    htseq_counts_df = htseq_counts_df[~htseq_counts_df['Gene_name'].apply(
        lambda x: x.startswith('__'))]

    gene_reads = htseq_counts_df['count'].sum()

    if gene_reads > filter_treshold * 1e6:
        print(True)
    else:
        print(False)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    htseq_counts = args.htseq_counts
    filter_treshold = int(args.filter_treshold)

    main(htseq_counts, filter_treshold)
