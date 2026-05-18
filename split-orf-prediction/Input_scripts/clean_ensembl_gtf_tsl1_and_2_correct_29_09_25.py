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
    ens_genes = ensembl_gtf.get_gn_ids(nr=True)
    ens_transcripts = ensembl_gtf.get_tx_ids(nr=True)
    filtered_gtf = ensembl_gtf.select_by_key('transcript_support_level', '1,2,1 (assigned to previous version 1),1 (assigned to previous version 2),1 (assigned to previous version 3),1 (assigned to previous version 4),1 (assigned to previous version 5),1 (assigned to previous version 6),1 (assigned to previous version 7),1 (assigned to previous version 8),1 (assigned to previous version 9),1 (assigned to previous version 10),1 (assigned to previous version 11),1 (assigned to previous version 12),1 (assigned to previous version 13),1 (assigned to previous version 14),1 (assigned to previous version 15),1 (assigned to previous version 16),2 (assigned to previous version 1),2 (assigned to previous version 2),2 (assigned to previous version 3),2 (assigned to previous version 4),2 (assigned to previous version 5),2 (assigned to previous version 6),2 (assigned to previous version 7),2 (assigned to previous version 8),2 (assigned to previous version 9),2 (assigned to previous version 10),2 (assigned to previous version 11),2 (assigned to previous version 12),2 (assigned to previous version 13),2 (assigned to previous version 14),2 (assigned to previous version 15)')
    genes = filtered_gtf.get_gn_ids(nr=True)
    transcripts = filtered_gtf.get_tx_ids(nr=True)
    gene_string = get_transcript_string(genes)
    gene_gtf = ensembl_gtf.select_by_key(
        'feature', 'gene').select_by_key('gene_id', gene_string)

    print('Number of transcripts in clean gtf:', len(transcripts))
    print("Number of genes in cleaned gtf:", len(genes))

    print('Number of transcripts in normal gtf:', len(ens_transcripts))
    print("Number of genes in normal gtf:", len(ens_genes))

    filtered_gtf.write(sys.argv[2])
    gene_gtf.write(sys.argv[3])


if __name__ == "__main__":
    main()
