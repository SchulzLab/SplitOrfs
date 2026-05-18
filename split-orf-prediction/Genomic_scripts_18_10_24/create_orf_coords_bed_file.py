"""
create_ORF_coords_bed_file.py - create a bed file of all valid ORFs in transcript coordinates


This script takes hte UniqueProteinORFPairs.txt as input.
Usage:
    python create_ORF_coords_bed_file.py UniqueProteinORFPairs.txt ORF_transcript_coords.bed
"""


import sys
import pandas as pd

ORFpairs = sys.argv[1]
bed_outfile = sys.argv[2]


ORFpairs_df = pd.read_csv(ORFpairs, sep='\t', header=0)

ORFpairs_df['ID'] = ORFpairs_df['geneID'] + '|' + ORFpairs_df['OrfTransID']
ORFpairs_df['OrfPos'] = ORFpairs_df['OrfPos'].apply(lambda x: x.split(','))
ORFpairs_df['OrfIDs'] = ORFpairs_df['OrfIDs'].apply(lambda x: x.split(','))

ORFpairs_df = ORFpairs_df[['ID', 'OrfPos', 'OrfIDs']]
ORFpairs_df = ORFpairs_df.explode(["OrfIDs", "OrfPos"], ignore_index=True)

ORFpairs_df['start'] = ORFpairs_df['OrfPos'].apply(
    lambda x: int(x.split('-')[0]))
ORFpairs_df['stop'] = ORFpairs_df['OrfPos'].apply(
    lambda x: int(x.split('-')[1]))

ORFpairs_df = ORFpairs_df[['ID', 'start', 'stop', 'OrfIDs']]

ORFpairs_df.to_csv(bed_outfile, sep="\t", index=False, header=False)
