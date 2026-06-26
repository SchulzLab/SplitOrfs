import os
import pandas as pd
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("reference_protein_fasta",
                        help="Path to reference protein FASTA")

    parser.add_argument("merged_splitorf_protein_fasta",
                        help="Path to merged Split-ORF protein FASTA")

    parser.add_argument("outfile",
                        help="Path to output file")

    return parser.parse_args()


def main(reference_protein_fasta, merged_splitorf_protein_fasta, outfile):

    # get output directory
    # outdir = os.path.dirname(outfile)

    # ref_seqs = SeqIO.parse(reference_protein_fasta, "fasta")
    # ref_seqs_set = set(str(record.seq) for record in ref_seqs)

    # so_seqs = SeqIO.parse(merged_splitorf_protein_fasta, "fasta")
    # so_seqs_set = set(str(record.seq) for record in so_seqs)

    # ref_seqs_set.isdisjoint(so_seqs_set)
    # FALSE

    # parse the reference and generate a unique mapping
    dedup_orfs = defaultdict(list)

    for record in SeqIO.parse(reference_protein_fasta, "fasta"):
        dedup_orfs[str(record.seq)].append(
            '|'.join(record.id.split('|')[1:3]))

    # get IDs that belong together as list of lists
    ref_id_list_of_lists = [ids for ids in dedup_orfs.values()]
    id_list = list(range(1, len(ref_id_list_of_lists)+1))
    id_list = ['ReferenceProtein' + '-' + str(id) for id in id_list]

    # create a dict of the original protein name as value and the
    # new reference name  ReferenceProtein-Nr as key
    idx = 0
    ref_id_to_unique_mapping_dict = {}
    for ref_id_list in ref_id_list_of_lists:
        for ref_id in ref_id_list:
            ref_id_to_unique_mapping_dict[ref_id] = id_list[idx]
        idx += 1

    ref_id_df = pd.DataFrame(ref_id_to_unique_mapping_dict,
                             index=['Reference_unique_ID']).transpose()
    ref_tsv_name = outfile[:-6] + '_ref_id_mapping.tsv'
    ref_id_df.to_csv(ref_tsv_name, sep='\t')

    # get SO dict of already deduplicated sequences: seq:SO-ID
    so_dedup_orfs = defaultdict(list)
    for record in SeqIO.parse(merged_splitorf_protein_fasta, "fasta"):
        so_dedup_orfs[str(record.seq)].append(
            '|'.join(record.id.split('|')[1:3]))

    # merge the df such that each SO-ID is the Index and the Reference_unique_ID
    merged_id_df = ref_id_df.copy()
    merged_seqs = dedup_orfs
    for sequence, so_list in so_dedup_orfs.items():
        if sequence not in dedup_orfs.keys():
            merged_seqs[sequence] = so_list
            merged_id_df.loc[so_list[0]] = [so_list[0]]

    orfs_list = [SeqRecord(Seq(seqi), id="sp" + "|" + merged_id_df.loc[id_list[0], 'Reference_unique_ID'],
                           name='', description='') for seqi, id_list in merged_seqs.items()]

    SeqIO.write(orfs_list, outfile, 'fasta')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    reference_protein_fasta = args.reference_protein_fasta
    merged_splitorf_protein_fasta = args.merged_splitorf_protein_fasta
    outfile = args.outfile

    main(reference_protein_fasta, merged_splitorf_protein_fasta, outfile)
