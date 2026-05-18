# ----- Make txt file of all genes provdided in which SOs are searched ----- #
# ----- These are the correct genes for GO background ----- #

import argparse
from Bio import SeqIO


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get list of all genes with Input transcripts without duplicates.")

    parser.add_argument(
        "transcript_fasta",
        type=str,
        help="Path to the Transcript FASTA file for the Input transcripts"
    )

    parser.add_argument(
        "background_gene_file",
        type=str,
        help="Path to the output file with the background genes (txt)"
    )

    return parser.parse_args()


def main(background_gene_file, transcript_fasta):
    headers = []
    with open(transcript_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            headers.append(record.description)

    input_genes = list(set([header.split('|')[0] for header in headers]))

    with open(background_gene_file, 'w') as f:
        for gene in input_genes:
            f.write(f"{gene}\n")


if __name__ == "__main__":
    args = parse_arguments()
    background_gene_file = args.background_gene_file
    transcript_fasta = args.transcript_fasta

    main(background_gene_file, transcript_fasta)
