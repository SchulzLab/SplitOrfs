import os
import pandas as pd
import argparse
import seaborn as sbn
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('three_prime_coverage_file',
                        help='Path to three prime coverage file')

    parser.add_argument('cds_coverage_file',
                        help='Path to three cds coverage file')

    parser.add_argument('sample',
                        help='sample')

    parser.add_argument('three_prime_original_for_strand_info',
                        help='Path to 3 prime bed with strand information')

    parser.add_argument('alpha',
                        help='alpha')

    return parser.parse_args()


def main(three_prime_coverage_file, sample, cds_coverage_file, three_prime_original_for_strand_info, alpha):
    def plot_hist(x, bins, title, start, stop, outdir, outname, df):
        ax = sbn.histplot(
            df, x=x, bins=bins)
        ax.set_xlim(start, stop)
        ax.set_title(title)
        fig = ax.get_figure()
        fig.savefig(os.path.join(outdir,
                    outname), dpi=300, bbox_inches='tight')
        plt.close(fig)

    outdir = os.path.dirname(three_prime_coverage_file)

    three_prime_coverage_file_df = pd.read_csv(
        three_prime_coverage_file, sep='\t', header=None,
        names=['chr', 'start', 'stop', '3_prime_name', 'nr_overlap_reads',
               'nr_bases_covered', 'length_3_prime_UTR', 'covered_fraction'],
        low_memory=False)

    # length normalization of reads
    three_prime_coverage_file_df['RPK'] = three_prime_coverage_file_df['nr_overlap_reads'] / \
        (three_prime_coverage_file_df['length_3_prime_UTR']/1000)

    cds_coverage_file_df = pd.read_csv(
        cds_coverage_file, sep='\t', header=None,
        names=['chr', 'start', 'stop', '3_prime_name', 'nr_overlap_reads',
               'nr_bases_covered', 'length_cds', 'covered_fraction'],
        low_memory=False)

    # length normalization of reads
    cds_coverage_file_df['RPK'] = cds_coverage_file_df['nr_overlap_reads'] / \
        (cds_coverage_file_df['length_cds']/1000)

    # plot RPK distribution CDS hist
    plot_hist('RPK', 50000, f'CDS RPK for {sample}', 0, 3000, outdir,
              f'{sample}_CDS_RPK_hist_0_2000.png', cds_coverage_file_df)

    # plot RPK distribution 3primes hist
    plot_hist('RPK', 100, f'3 prime RPK for {sample}', 0, 1000, outdir,
              f'{sample}_three_primes_RPK_hist_0_1000.png',
              three_prime_coverage_file_df[three_prime_coverage_file_df['RPK'] < 2000])

    plot_hist('RPK', 100, f'3 prime RPK > 0 for {sample}', 0, 1000, outdir, f'{sample}_three_primes_RPK_hist_greater_0_1000.png',
              three_prime_coverage_file_df[(three_prime_coverage_file_df['RPK'] > 0) &
                                           (three_prime_coverage_file_df['RPK'] < 2000)])

    # FILTER OUT VALUES THAT ARE 0
    cds_coverage_file_df_filtered = cds_coverage_file_df[cds_coverage_file_df['RPK'] > 0]
    PI_lower = cds_coverage_file_df_filtered['RPK'].quantile(alpha)

    # here we can keep the 0-counts as we anyway want to have them
    # select all the values below the PI_lower
    three_prime_coverage_file_df_filtered = three_prime_coverage_file_df[
        three_prime_coverage_file_df['RPK'] < PI_lower].copy()

    # how does the three prime RPK distribution look like after the filtering?
    print('Fitlering threshold for RPK:', PI_lower)

    plot_hist('RPK', 50, f'CDS filtered 3 prime RPK for {sample}', 0, 300,
              outdir, f'{sample}_three_primes_RPK_after_CDS_filter_hist_greater_0_300.png',
              three_prime_coverage_file_df_filtered[three_prime_coverage_file_df_filtered['RPK'] > 0])

    # need to recover the strand information, based on gene
    strand_info_df = pd.read_csv(three_prime_original_for_strand_info, sep='\t',
                                 names=[
                                     'chr',
                                     'start',
                                     'stop',
                                     'name',
                                     'score',
                                     'strand'])

    strand_info_df['Gene'] = strand_info_df['name'].apply(
        lambda x: x.split('|')[0])
    strand_dict = dict(zip(strand_info_df['Gene'], strand_info_df['strand']))

    three_prime_coverage_file_df_filtered['Gene'] = \
        three_prime_coverage_file_df_filtered['3_prime_name'].apply(
        lambda x: x.split('|')[0])
    three_prime_coverage_file_df_filtered['strand'] = \
        three_prime_coverage_file_df_filtered['Gene'].map(
        strand_dict)
    three_prime_coverage_file_df_filtered['score'] = 0

    three_prime_coverage_file_df_filtered.loc[:,
                                              ['chr',
                                               'start',
                                               'stop',
                                               '3_prime_name',
                                               'score',
                                               'strand']
                                              ].to_csv(os.path.join(
                                                  outdir,
                                                  f'3_primes_filtered_for_CDS_distribution_{sample}.bed'),
        sep='\t', header=False, index=False)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    three_prime_coverage_file = args.three_prime_coverage_file
    sample = args.sample
    cds_coverage_file = args.cds_coverage_file
    three_prime_original_for_strand_info = args.three_prime_original_for_strand_info
    alpha = float(args.alpha)

    # three_prime_coverage_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367798/3_primes_genomic_merged_numbered_ERR3367798_windows_coverage.tsv'
    # cds_coverage_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367798/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered_ERR3367798_windows_coverage.tsv'
    # sample = 'ERR3367798'
    # three_prime_original_for_strand_info = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367798/3_primes_genomic_merged_numbered_ERR3367798.bed'
    # alpha = 0.05

    main(three_prime_coverage_file, sample, cds_coverage_file,
         three_prime_original_for_strand_info, alpha)
