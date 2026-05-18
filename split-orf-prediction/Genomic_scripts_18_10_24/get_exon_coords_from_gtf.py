# Take any GTF file and extract the exon coordinates of the transcript
# in the format required by the Split-ORF pipeline

# usage get_exon_coords_from_gtf.py in.gtf out.txt


import argparse
import pandas as pd
from pygtftk.gtf_interface import GTF


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get exon coordinates of transcript for Split-ORF pipeline.")

    parser.add_argument(
        "gtf_file",
        type=str,
        help="Path to the input GTF file"
    )

    parser.add_argument(
        "output_txt_file",
        type=str,
        help="Path to the output file with the exon coords (txt)"
    )

    return parser.parse_args()


def main(path_to_custom_gtf, outfile):

    # load custom gtf
    custom_gtf = GTF(path_to_custom_gtf, check_ensembl_format=False)

    # get transcript info
    mando_transcript_info = custom_gtf.select_by_key(
        "feature", "transcript").extract_data("gene_id,transcript_id,start,end,strand,chr")

    mando_transcript_info_df = mando_transcript_info.as_data_frame()

    # get exon info
    mando_exon_info = custom_gtf.select_by_key(
        "feature", "exon").extract_data("transcript_id,start,end")

    mando_exon_info_df = mando_exon_info.as_data_frame()

    # merge dataframes based on transcript ID
    transcript_exon_df = pd.merge(mando_exon_info_df, mando_transcript_info_df,
                                  on='transcript_id', suffixes=('_exon', '_transcript'), how='outer')

    # names do not matter as the header is filtered out in the SO pipeline
    # entries are divided into belonging to the positive and negative strand

    # should add one to the start coordinate as the SO pipeline
    # expects Ensembl coords which are 1-based and adapts GTF files are also 1-based
    # ExonToTranscriptPositions.py script
    transcript_exon_df['start_exon'] = transcript_exon_df['start_exon'].astype(
        int)
    transcript_exon_df['start_transcript'] = transcript_exon_df['start_transcript'].astype(
        int)

    # reorder columns to the correct order
    transcript_exon_df = transcript_exon_df[['gene_id', 'transcript_id', 'start_exon',
                                            'end_exon', 'start_transcript', 'end_transcript', 'strand', 'seqid']]

    transcript_exon_df.to_csv(outfile, sep='\t', index=False)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.gtf_file, args.output_txt_file)
