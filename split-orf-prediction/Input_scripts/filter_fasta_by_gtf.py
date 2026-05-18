####filters a FASTA file for transcripts present in a gtf file#####
import sys
from pygtftk.gtf_interface import GTF
from Bio import SeqIO

filtered_gtf = GTF(sys.argv[1], check_ensembl_format=False)
tids_keep = filtered_gtf.get_tx_ids(nr=True)

sequences_keep = []
for record in SeqIO.parse(sys.argv[2], "fasta"):
    if record.description.split('|')[1] in tids_keep:
        sequences_keep.append(record)


print('number of sequences to keep:', len(sequences_keep))

with open(sys.argv[3], "w") as output_handle:
    SeqIO.write(sequences_keep, output_handle, "fasta")