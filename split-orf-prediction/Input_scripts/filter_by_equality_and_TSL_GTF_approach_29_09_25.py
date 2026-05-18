# We read in 3 GTF files
# the Ensembl GTF to filter, the Refseq equality and the TSL1 and 2 filtered ones
# we only keep transcripts in the Ensembl GTF that are also present in the
# at least one of the other GTF files
import sys
from typing import List
from pygtftk.gtf_interface import GTF


def get_transcript_string(transcript_ids: List[str]) -> str:
    '''get string of transcript or other ids for filtering 
    of a GTF class object from pygtftk'''
    transcript_string = ''
    if len(transcript_ids) > 1:
        for transcript_id in transcript_ids[:-1]:
            transcript_string += transcript_id+','
    transcript_string += transcript_ids[-1]
    return transcript_string

# select all transcript and exon entries from ensembl gtf and write to a new gtf file


def main():

    ensembl_gtf = GTF(sys.argv[1], check_ensembl_format=False)
    tsl_filtered_gtf = GTF(sys.argv[2], check_ensembl_format=False)
    equality_filtered_gtf = GTF(sys.argv[3], check_ensembl_format=False)

    tsl_filtered_genes = tsl_filtered_gtf.get_gn_ids(nr=True)
    tsl_filtered_transcripts = tsl_filtered_gtf.get_tx_ids(nr=True)

    equality_filtered_genes = equality_filtered_gtf.get_gn_ids(nr=True)
    equality_filtered_transcripts = equality_filtered_gtf.get_tx_ids(nr=True)

    transcripts_to_keep = list(
        set(tsl_filtered_transcripts + equality_filtered_transcripts))

    genes_to_keep = list(set(tsl_filtered_genes + equality_filtered_genes))

    transcript_string = get_transcript_string(transcripts_to_keep)

    gene_string = get_transcript_string(genes_to_keep)
    gene_gtf = ensembl_gtf.select_by_key(
        'feature', 'gene').select_by_key('gene_id', gene_string)

    transcript_gtf = ensembl_gtf.select_by_key(
        'feature', 'transcript,exon,CDS').select_by_key('transcript_id', transcript_string)

    print('Number of transcripts in clean gtf:', len(transcripts_to_keep))
    print("Number of genes in cleaned gtf:", len(genes_to_keep))

    transcript_gtf.write(sys.argv[4])
    gene_gtf.write(sys.argv[5])


if __name__ == "__main__":
    main()
