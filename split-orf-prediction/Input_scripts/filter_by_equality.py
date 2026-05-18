# This script takes a gtf of an assembly from stringtie that was merged of single sample assemblies
# and writes the ref_gene_id attribute as the gene_id
import sys
import regex as re
import pandas as pd
from pygtftk.gtf_interface import GTF

Ensembl_gtf = open(sys.argv[1], 'r')
Lines = Ensembl_gtf.readlines()
refseq_anno_file = pd.read_csv(sys.argv[2], sep="\t", comment='#', header=0)
refseq_equal_matches = refseq_anno_file[refseq_anno_file["class_code"] == "="]
print(refseq_equal_matches.head())
tids_matching_list = set(refseq_equal_matches['qry_id'].to_list())
print('Number of transcripts that have the = sign and will be kept in the Ensembl annotation',
      len(set(tids_matching_list)))
gtf_with_repalcement = open(sys.argv[3], 'w')


replaced_lines = []

for line in Lines:
    if not line.startswith('#'):
        trans_id_match = re.search(r'transcript_id "([^"]+)"', line)
        if trans_id_match:
            ref_trans_id = trans_id_match.group(1)
            if ref_trans_id in tids_matching_list:
                replaced_lines.append(line)
    else:
        replaced_lines.append(line)

gtf_with_repalcement.writelines(replaced_lines)
