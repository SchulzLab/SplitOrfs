#This script takes a gtf of an assembly from stringtie that was merged of single sample assemblies
#and writes the ref_gene_id attribute as the gene_id
import sys
import regex as re
import pandas as pd
from pygtftk.gtf_interface import GTF

refseq_gtf = open(sys.argv[1], 'r')
Lines = refseq_gtf.readlines()
refseq_anno_file = pd.read_csv(sys.argv[2], sep="\t", comment='#', header = 0)
#print(refseq_anno_file.head(50))

#create a dictionary for the chromosome mapping
chr_dict= refseq_anno_file.set_index('RefSeq-Accn')['Sequence-Name'].to_dict()
#print(chr_dict)

gtf_with_repalcement = open(sys.argv[3], 'w')


replaced_lines = []

for line in Lines:
    if not line.startswith('#'):
        columns = re.split(r'\s+', line.strip())
        line = re.sub(columns[0], chr_dict[columns[0]], line)
        #print(line)
    replaced_lines.append(line)

gtf_with_repalcement.writelines(replaced_lines)