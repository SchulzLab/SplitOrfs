import sys
import re
from Bio import SeqIO

# this script extracts the protein sequences of the valid SplitORF transcripts
# for MassSpec search
# input files are:
# 1: file with the valid SO (.txt)
# 2: file with all possible ORFs (fasta)
# 3: Output file for protein sequences (fasta)


def extract_Split_ORF_ids(file1):
    Valid_SOs = open(file1, 'r')
    Split_orf_ids = {}
    for SplitORF in Valid_SOs:
        SplitORF = SplitORF.strip()
        SplitORF = SplitORF.split('\t')
        Split_orf_ids[SplitORF[2]] = SplitORF[4].split(',')
    Valid_SOs.close()
    return Split_orf_ids


def get_fasta_sequences(Split_orf_ids, fasta_name, outname):
    fin = open(fasta_name, 'r')
    fout = open(outname, 'w')
    for record in SeqIO.parse(fin, 'fasta'):
        header_info = re.split(r'[|:]', record.id)
        if header_info[1] in Split_orf_ids.keys():
            for transcript_id, list_ORFs in Split_orf_ids.items():
                if list_ORFs:
                    if transcript_id == header_info[1]:
                        # print(transcript_id,  header_info[1])
                        for ORF in list_ORFs:
                            if ORF == header_info[2]:
                                # print(transcript_id,  header_info[1], ORF, header_info[2])
                                record.id = 'sp|' + record.id + '|'
                                record.description = ''
                                SeqIO.write(record, fout, "fasta")
    fin.close()
    fout.close()


if len(sys.argv) < 4:
    print("usage get_protein_seqs_for_Masspec.py UniqueProteinORFPairs.txt OrfProteins.fa Proteins_with_unique_regions_for_masspec.fa")
else:
    Split_ORF_ids = extract_Split_ORF_ids(sys.argv[1])
    get_fasta_sequences(Split_ORF_ids, sys.argv[2], sys.argv[3])
