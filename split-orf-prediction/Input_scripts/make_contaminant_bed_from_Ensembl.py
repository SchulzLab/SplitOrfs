#### makes a BED file from Ens CDS coordinates download #####
#### filters for transcripts within a GTF and being protein coding #####

from pygtftk.gtf_interface import GTF
import argparse
import pandas as pd
from pybedtools import BedTool


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("--ensembl_genomic_coords_txt",
                        help="TXT file of Ensembl CDS coordinates")

    parser.add_argument("--output_bed_file",
                        help="Path to output BED file")

    return parser.parse_args()


def main(ensembl_genomic_coords_txt, output_bed_file):
    # read file with CDS cooridnates
    ensembl_genomic_cds_coords_df = pd.read_csv(
        ensembl_genomic_coords_txt, sep='\t')

    ensembl_genomic_cds_coords_df = ensembl_genomic_cds_coords_df[~ensembl_genomic_cds_coords_df['Exon region start (bp)'].isna(
    )]

    ensembl_genomic_cds_coords_df['ID'] = ensembl_genomic_cds_coords_df['Gene stable ID'] + \
        '|' + ensembl_genomic_cds_coords_df['Transcript stable ID']

    def change_strand(strand):
        if strand == 1:
            return '+'
        elif strand == -1:
            return '-'
        else:
            return '.'

    ensembl_genomic_cds_coords_df['Strand'] = ensembl_genomic_cds_coords_df['Strand'].apply(
        lambda x: change_strand(x))

    coordinate_string = ''
    for line in ensembl_genomic_cds_coords_df.index:
        coordinate_string += str(ensembl_genomic_cds_coords_df.loc[line, 'Chromosome/scaffold name'])+'\t' + str(int(ensembl_genomic_cds_coords_df.loc[line, 'Exon region start (bp)']) - 1)+'\t' +\
            str(int(ensembl_genomic_cds_coords_df.loc[line, 'Exon region end (bp)'])) + '\t' +\
            ensembl_genomic_cds_coords_df.loc[line, 'ID'] + '\t' +\
            str(0) + '\t' + \
            str(ensembl_genomic_cds_coords_df.loc[line, 'Strand']) + '\n'

    coordinate_bed = BedTool(
        coordinate_string, from_string=True).saveas(output_bed_file)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    ensembl_genomic_coords_txt = args.ensembl_genomic_coords_txt
    output_bed_file = args.output_bed_file

    main(ensembl_genomic_coords_txt, output_bed_file)
