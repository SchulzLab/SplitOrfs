import sys
import pandas as pd

# NMD test cases
# ENSG00000164082|ENST00000296479	1565	1662	ORF-440:1565:2993
# end of region is end of exon: need to put end coordinate as the first base of the
# neighboring intron
# ENST00000309608 ORF-789
# unique region crossing exon borders and on minus strand
# ENST00000306450:ORF-907
# positive strand covering 3 exons
# ENSG00000007341|ENST00000361846:ORF-1506
# minus strand with the start coordinate being the start of the exon (end in transcript universe)

# RI test cases


transcript_coords_test = pd.read_csv(
    sys.argv[1], sep='\t', header=None, names=['chr', 'gen_start', 'gen_stop', 'ID',  'phase', 'strand'])

genomic_pos_calculated_all = pd.read_csv(
    sys.argv[2], sep='\t', header=None, names=['chr', 'gen_start', 'gen_stop', 'ID',  'phase', 'strand'])

IDs_to_compare = transcript_coords_test['ID'].to_list()
genomic_pos_subset = genomic_pos_calculated_all[genomic_pos_calculated_all['ID'].isin(
    IDs_to_compare)]
IDs_used = genomic_pos_subset['ID'].to_list()
transcript_coords_test = transcript_coords_test[transcript_coords_test['ID'].isin(
    IDs_used)]

genomic_pos_subset = genomic_pos_subset.reset_index(drop=True)
transcript_coords_test = transcript_coords_test.reset_index(drop=True)

genomic_pos_subset = genomic_pos_subset.astype(transcript_coords_test.dtypes)
# print(genomic_pos_subset)
# print(transcript_coords_test)
assert genomic_pos_subset.equals(transcript_coords_test)
