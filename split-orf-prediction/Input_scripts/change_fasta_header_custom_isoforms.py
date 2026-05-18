import sys
from pygtftk.gtf_interface import GTF
from Bio import SeqIO


ensembl_gtf = GTF(sys.argv[1], check_ensembl_format=False)
custom_fasta = sys.argv[2]
out_fasta = sys.argv[3]


def get_tid_gid_dict(ensembl_gtf):
    # this is a dict with genes as keys as transcripts as values
    gene_tid_dict = ensembl_gtf.get_gn_to_tx()

    # convert to tid-geneid dict
    tis_gene_dict = {}

    for gid, tids in gene_tid_dict.items():
        # are all lists, even if only one tid present
        if isinstance(tids, list):
            for tid in tids:
                tis_gene_dict[tid] = gid
    return tis_gene_dict


def modify_header(header, tis_gene_dict):
    if header in tis_gene_dict.keys():
        modified_header = tis_gene_dict[header] + '|' + header
    else:
        modified_header = '.'.join(header.split('.')[:-1]) + '|' + header
    return modified_header


tis_gene_dict = get_tid_gid_dict(ensembl_gtf)

# Read, modify, and write new FASTA
with open(out_fasta, "w") as out_fasta:
    for record in SeqIO.parse(custom_fasta, "fasta"):
        record.id = modify_header(record.id, tis_gene_dict)  # Modify header
        record.description = ""  # Remove additional descriptions if needed
        SeqIO.write(record, out_fasta, "fasta")  # Write modified record
