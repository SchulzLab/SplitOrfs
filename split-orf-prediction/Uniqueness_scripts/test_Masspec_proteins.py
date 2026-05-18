import sys
import re
from Bio import SeqIO

# this script implements the test of the completeness and correctness of the
# SO proteins used for Masspec by double checking with another file
# for MassSpec search
# input files are:
# 1: Valid_ORF_proteins.bed (bed file with all valid SOs: for double check)
# 2: Masspec FASTA


def get_SO_IDs_bed(file1):
    Valid_SOs = open(file1, 'r')
    Split_orf_ids = []
    for SplitORF in Valid_SOs:
        SplitORF = SplitORF.strip()
        SplitORF = SplitORF.split('\t')
        Split_orf_ids.append(SplitORF[0])
    Valid_SOs.close()
    return Split_orf_ids


def compare_SO_IDs_bed_Masspec_FASTA(Split_orf_ids, fasta_name):
    fin = open(fasta_name, 'r')
    SO_headers = []
    for record in SeqIO.parse(fin, 'fasta'):
        header_info = record.id[3:-1]
        SO_headers.append(header_info)
    fin.close()
    assert SO_headers == Split_orf_ids


if len(sys.argv) < 3:
    print("usage test_Masspec_proteins.py Valid_ORF_proteins.bed Proteins_with_unique_regions_for_masspec.fa")
else:
    Split_orf_ids = get_SO_IDs_bed(sys.argv[1])
    compare_SO_IDs_bed_Masspec_FASTA(Split_orf_ids, sys.argv[2])
