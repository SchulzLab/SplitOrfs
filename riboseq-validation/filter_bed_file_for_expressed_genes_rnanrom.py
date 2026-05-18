import os
import pandas as pd
import argparse
from rnanorm import FPKM, TPM
from sklearn import set_config
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np
# set_config(transform_output='pandas')


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('bed_file_to_filter',
                        help='Path to BED file to filter')

    parser.add_argument('htseq_counts',
                        help='Path to HTSeq-counts')

    parser.add_argument('gtf_path',
                        help='Path to GTF')

    parser.add_argument('tpm_threshold',
                        help='TPM Threshold')

    return parser.parse_args()


def main(bed_file_to_filter, htseq_counts, gtf_path, tpm_threshold):
    outdir = os.path.dirname(htseq_counts)
    for file_end in ['_NMD_htseq_counts.tsv', '_RI_htseq_counts.tsv']:
        if os.path.basename(htseq_counts).endswith(file_end):
            sample = os.path.basename(htseq_counts).removesuffix(file_end)
    bed_file_name = os.path.basename(bed_file_to_filter).removesuffix('.bed')

    bed_file_to_filter_df = pd.read_csv(
        bed_file_to_filter, sep='\t', header=None, low_memory=False)
    bed_file_to_filter_df['Gene_name'] = bed_file_to_filter_df.iloc[:, 3].apply(
        lambda x: x.split('|')[0])

    htseq_counts_df = pd.read_csv(
        htseq_counts, sep='\t', header=None, names=['Gene_name', 'count'])
    # filter out for the non-gene columns
    htseq_counts_df = htseq_counts_df[~htseq_counts_df['Gene_name'].apply(
        lambda x: x.startswith('__'))]

    htseq_counts_df_transposed = htseq_counts_df.transpose()

    htseq_counts_df_transposed.columns = htseq_counts_df_transposed.loc['Gene_name', :]

    htseq_counts_df_transposed = htseq_counts_df_transposed.drop('Gene_name')

    tpm_values = TPM(gtf=gtf_path).fit_transform(htseq_counts_df_transposed)

    htseq_counts_df_transposed_tpm = pd.concat([htseq_counts_df_transposed, pd.DataFrame(
        [tpm_values.flatten(), np.log10(tpm_values.flatten() + 1)], columns=htseq_counts_df_transposed.columns)], ignore_index=True)

    htseq_counts_df_transposed_tpm.index = [
        'raw_counts', 'TPM', 'log10(TPM + 1)']

    tpm_greater_1_df = htseq_counts_df_transposed_tpm.transpose(
    )[(htseq_counts_df_transposed_tpm.transpose()['TPM'] >= 1) & (htseq_counts_df_transposed_tpm.transpose()['TPM'] < 1000)]

    ax = sbn.histplot(
        tpm_greater_1_df, x='TPM', bins=1000)
    ax.set_xlim(0, 200)
    ax.set_title(f'TPM values per gene (1-200) for {sample}')
    plt.show()

    fig = ax.get_figure()  # get the figure that contains the axes
    fig.savefig(os.path.join(outdir,
                f'{sample}_TPM_1_200.png'), dpi=300, bbox_inches='tight')  # save as PNG

    # filter the regions for the genes with TPM above TPM Treshold
    htseq_counts_df_fpkm = htseq_counts_df_transposed_tpm.transpose()
    genes_to_keep = htseq_counts_df_fpkm.index[htseq_counts_df_fpkm['TPM']
                                               >= tpm_threshold]

    bed_file_filtered = bed_file_to_filter_df[bed_file_to_filter_df['Gene_name'].isin(
        genes_to_keep)]
    bed_file_filtered.iloc[:, 0:6].to_csv(os.path.join(
        outdir, f'{bed_file_name}_{sample}.bed'), sep='\t', header=False, index=False)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    bed_file_to_filter = args.bed_file_to_filter
    htseq_counts = args.htseq_counts
    gtf_path = args.gtf_path
    tpm_threshold = int(args.tpm_threshold)

    # bed_file_to_filter = '/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged_numbered.bed'

    # htseq_counts = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/huvec_dnor_2/huvec_dnor_2_NMD_htseq_counts.tsv'

    # gtf_path = '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'
    main(bed_file_to_filter, htseq_counts, gtf_path, tpm_threshold)
