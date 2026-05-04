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
    parser.add_argument("ensembl_genomic_cds_coords_txt",
                        help="TXT file of Ensembl CDS coordinates")

    parser.add_argument("gtf_file",
                        help="GTF file to use for filtering")

    parser.add_argument("output_bed_file",
                        help="Path to output BED file")

    return parser.parse_args()


def main(ensembl_genomic_cds_coords_txt, output_bed_file, gtf_file):
    # get transcripts to keep from the GTF
    filtered_gtf = GTF(gtf_file, check_ensembl_format=False)
    all_tids_gtf = filtered_gtf.get_tx_ids(nr=True)
    filtered_gtf_prot_cod = filtered_gtf.select_by_key(
        'transcript_biotype', 'protein_coding')
    tids_keep = filtered_gtf_prot_cod.get_tx_ids(nr=True)

    # read file with CDS cooridnates
    ensembl_genomic_cds_coords_df = pd.read_csv(
        ensembl_genomic_cds_coords_txt, sep='\t')

    ensembl_genomic_cds_coords_df = ensembl_genomic_cds_coords_df[~ensembl_genomic_cds_coords_df['Genomic coding start'].isna(
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

    # select transcripts from GTF
    ensembl_genomic_cds_coords_df = ensembl_genomic_cds_coords_df[ensembl_genomic_cds_coords_df['Transcript stable ID'].isin(
        tids_keep)]

    coordinate_string = ''
    for line in ensembl_genomic_cds_coords_df.index:
        coordinate_string += str(ensembl_genomic_cds_coords_df.loc[line, 'Chromosome/scaffold name'])+'\t' + str(int(ensembl_genomic_cds_coords_df.loc[line, 'Genomic coding start']) - 1)+'\t' +\
            str(int(ensembl_genomic_cds_coords_df.loc[line, 'Genomic coding end'])) + '\t' +\
            ensembl_genomic_cds_coords_df.loc[line, 'ID'] + '\t' +\
            str(0) + '\t' + \
            str(ensembl_genomic_cds_coords_df.loc[line, 'Strand']) + '\n'

    coordinate_bed = BedTool(
        coordinate_string, from_string=True).saveas(output_bed_file)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    ensembl_genomic_cds_coords_txt = args.ensembl_genomic_cds_coords_txt
    output_bed_file = args.output_bed_file
    gtf_file = args.gtf_file

    # gtf_file = '/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf'
    # ensembl_genomic_cds_coords_txt = '/projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_all.txt'

    main(ensembl_genomic_cds_coords_txt, output_bed_file, gtf_file)
