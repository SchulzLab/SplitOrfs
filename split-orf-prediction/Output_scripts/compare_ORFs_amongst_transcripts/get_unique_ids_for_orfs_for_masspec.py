import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from collections import defaultdict
import pandas as pd


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

    parser.add_argument(
        "id_name",
        type=str,
        help="id_name of the newly assigned unique IDs"
    )

    parser.add_argument('--fasta3', default='',
                        help='Path to third FASTA')

    parser.add_argument('--fasta4', default='',
                        help='Path to fourth FASTA')

    return parser.parse_args()


def main(fasta1, fasta2, outname, id_name, fasta3='', fasta4=''):
    # count number of times a sequence has several identifiers
    dup_counter = 0
    # count total number of sequences in all FASTA files
    total_counter = 0
    fasta_list = [fasta1, fasta2, fasta3, fasta4]
    # create a dictionary that has the sequence as keys, and lists as values of the
    # respective SplitORF IDs
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

    # get IDs that belong together as list of lists
    so_id_list_of_lists = [ids for ids in dedup_orfs.values()]
    id_list = list(range(1, len(so_id_list_of_lists)+1))
    id_list = [id_name + '-' + str(id) for id in id_list]

    idx = 0
    so_id_to_unique_mapping_dict = {}
    for so_id_list in so_id_list_of_lists:
        for so_id in so_id_list:
            so_id_to_unique_mapping_dict[so_id] = id_list[idx]
        idx += 1

    so_id_df = pd.DataFrame(so_id_to_unique_mapping_dict,
                            index=['SO_unique_ID']).transpose()

    so_name = outname[:-6] + '_so_id_mapping.tsv'
    so_id_df.to_csv(so_name, sep='\t')

    orfs_list = [SeqRecord(Seq(seqi), id="sp" + "|" + so_id_df.loc[so_id_list[0], 'SO_unique_ID'],
                           name='', description='') for seqi, so_id_list in dedup_orfs.items()]

    SeqIO.write(orfs_list, outname, 'fasta')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_arguments()

    fasta1 = args.fasta1
    fasta2 = args.fasta2
    outname = args.outname
    id_name = args.id_name
    fasta3 = args.fasta3
    fasta4 = args.fasta4

    main(fasta1, fasta2, outname, id_name, fasta3, fasta4)
