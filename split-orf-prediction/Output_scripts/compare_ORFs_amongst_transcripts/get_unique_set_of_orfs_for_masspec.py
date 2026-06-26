from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from collections import defaultdict
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get all unique ORF protein sequences for Masspec across 2 or more fasta files.")

    parser.add_argument(
        "fasta1",
        type=str,
        help="Path to first FASTA"
    )

    parser.add_argument(
        "fasta2",
        type=str,
        help="Path to second FASTA"
    )

    parser.add_argument(
        "outname",
        type=str,
        help="outname of FASTA with unique ORF sequences"
    )

    parser.add_argument('--fasta3', default='',
                        help='Path to third FASTA')

    parser.add_argument('--fasta4', default='',
                        help='Path to fourth FASTA')

    return parser.parse_args()


def main(fasta1, fasta2, outname, fasta3='', fasta4=''):
    dup_counter = 0
    total_counter = 0
    fasta_list = [fasta1, fasta2, fasta3, fasta4]
    dedup_orfs = defaultdict(list)
    for fasta in fasta_list:
        if fasta:
            sample_counter = 0
            for record in SeqIO.parse(fasta, "fasta"):
                dedup_orfs[str(record.seq)].append(
                    '|'.join(record.id.split('|')[1:3]))
                sample_counter += 1
                if len(dedup_orfs[str(record.seq)]) > 1:
                    dup_counter += 1
            total_counter += sample_counter
            print('sample counter', sample_counter)
    print("dup_counter", dup_counter)
    print('total counter', total_counter)
    orfs_list = [SeqRecord(Seq(seqi), id="sp" + "|" + gi[0],
                           name='', description='') for seqi, gi in dedup_orfs.items()]

    SeqIO.write(orfs_list, outname, 'fasta')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_arguments()

    fasta1 = args.fasta1
    fasta2 = args.fasta2
    outname = args.outname
    fasta3 = args.fasta3
    fasta4 = args.fasta4

    main(fasta1, fasta2, outname, fasta3, fasta4)
