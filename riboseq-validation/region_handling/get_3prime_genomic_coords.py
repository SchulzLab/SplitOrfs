# get_CDS_coords.py is given structures information from Biomart on the CDS coordiantes
# required fields are: Gene ID, transcript ID, cDNA coding end, cDNA coding start and
# Transcript length (including UTRs and CDS)
# the exon wise information is combined such that the CDS start and end is extracted
# in transcript coordinates
# from this a bed file containing the 3' UTR regions is extracted for the random
# intersection in the Riboseq pipeline

import sys
import pandas as pd
from pybedtools import BedTool
from pygtftk.gtf_interface import GTF

# read in gtf file Ensembl that is tsl 1 and 2 filtered
Ensembl_tsl_filtered_gtf = GTF(sys.argv[1])
tsl_filtered_prot_cod_tids = Ensembl_tsl_filtered_gtf.select_by_key(
    'transcript_support_level', '1,1 (assigned to previous version 1),1 (assigned to previous version 2),1 (assigned to previous version 3),1 (assigned to previous version 4),1 (assigned to previous version 5),1 (assigned to previous version 6),1 (assigned to previous version 7),1 (assigned to previous version 8),1 (assigned to previous version 9),1 (assigned to previous version 10),1 (assigned to previous version 11),1 (assigned to previous version 12),1 (assigned to previous version 13),1 (assigned to previous version 14),1 (assigned to previous version 15),1 (assigned to previous version 16)').select_by_key('transcript_biotype', 'protein_coding').get_tx_ids(nr=True)


# read file with CDS cooridnates
Ensembl_genomic_three_primes = pd.read_csv(sys.argv[2], sep='\t')
Ensembl_genomic_three_primes = Ensembl_genomic_three_primes[~Ensembl_genomic_three_primes["3' UTR start"].isna(
)]
Ensembl_genomic_three_primes = Ensembl_genomic_three_primes[Ensembl_genomic_three_primes['Transcript stable ID'].isin(
    tsl_filtered_prot_cod_tids)]
Ensembl_genomic_three_primes['ID'] = Ensembl_genomic_three_primes['Gene stable ID'] + \
    '|' + Ensembl_genomic_three_primes['Transcript stable ID']
print(Ensembl_genomic_three_primes.columns)
Ensembl_genomic_three_primes = Ensembl_genomic_three_primes[[
    'Chromosome/scaffold name',  "3' UTR start",  "3' UTR end", 'Strand', 'ID']]

coordinate_string = ''
for line in Ensembl_genomic_three_primes.index:
    coordinate_string += str(Ensembl_genomic_three_primes.loc[line, 'Chromosome/scaffold name'])+'\t' + str(int(Ensembl_genomic_three_primes.loc[line, "3' UTR start"]))+'\t' +\
        str(int(Ensembl_genomic_three_primes.loc[line, "3' UTR end"])) + '\t' +\
        Ensembl_genomic_three_primes.loc[line, 'ID'] + '\t' +\
        str(0) + '\t' + \
        str(int(Ensembl_genomic_three_primes.loc[line, 'Strand'])) + '\n'


coordinate_bed = BedTool(
    coordinate_string, from_string=True).saveas(sys.argv[3])
