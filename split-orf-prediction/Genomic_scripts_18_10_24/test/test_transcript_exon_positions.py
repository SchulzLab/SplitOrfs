import sys
import pandas as pd

transcript_coords_test = pd.read_csv(
    sys.argv[1], sep='\t', header=None, names=['ID', 'gen_start', 'gen_stop', 'transcript_start', 'transcript_stop', 'strand', 'chr'])

transcript_coords_test['transcript_length'] = transcript_coords_test['transcript_stop'] - \
    transcript_coords_test['transcript_start']
transcript_coords_test['genomic_length'] = transcript_coords_test['gen_stop'] - \
    transcript_coords_test['gen_start']

# first check: genomic and transcriptomic regions need to have the same length
assert (transcript_coords_test['transcript_length']
        == transcript_coords_test['genomic_length']).all(), 'transcript and genomic lengths differ'


whole_transcript_length = transcript_coords_test.groupby(
    'ID')['transcript_length'].apply(lambda x: x.sum())


# second test: from ENsmebl archives v110, take the transcript lengths and compare
assert whole_transcript_length['ENSG00000163138|ENST00000508952'] == 987, 'transcript length of ENST00000508952 incorrect'
assert whole_transcript_length['ENSG00000126524|ENST00000697861'] == 1339, 'transcript length of ENST00000697861 incorrect'
assert whole_transcript_length['ENSG00000126524|ENST00000697862'] == 1532, 'transcript length of ENST00000697862 incorrect'


# for the minus strand take one transcript and check coordiantes with Ensembl directly
assert transcript_coords_test[transcript_coords_test['ID'] ==
                              'ENSG00000126524|ENST00000697862']['gen_start'].min() == 66987700 - 1
assert transcript_coords_test[transcript_coords_test['ID'] ==
                              'ENSG00000126524|ENST00000697862']['gen_stop'].max() == 66995604


transcript_coords_full = pd.read_csv(
    sys.argv[3], sep='\t', header=None, names=['ID', 'gen_start', 'gen_stop', 'transcript_start', 'transcript_stop', 'strand', 'chr'])

exon_coords = pd.read_csv(
    sys.argv[2], sep='\t', header=None, names=['gID', 'tID', 'gen_start', 'gen_stop', 'overall_start', 'overall_stop', 'strand', 'chr'])

# make the chr column as string, to be sure
transcript_coords_full['chr'] = transcript_coords_full['chr'].astype(str)
exon_coords['chr'] = exon_coords['chr'].astype(str)

# sort both dataframes by start coordinate to ensure the same ordering for the
# following comparison of start coordinates
transcript_coords_full = transcript_coords_full.sort_values(
    by=['chr', 'gen_start'])
exon_coords = exon_coords.sort_values(by=['chr', 'gen_start'])

assert (transcript_coords_full['gen_start']
        == exon_coords['gen_start'] - 1).all()

# Compare the start and end position of all calculated coordinates to the given ones
calculated_start_coords = transcript_coords_full.groupby('ID')[
    'gen_start'].min()
calculated_end_coords = transcript_coords_full.groupby('ID')['gen_stop'].max()

given_start_coords = exon_coords.groupby(
    'tID')['gen_start'].min().apply(lambda x: x-1)
given_end_coords = exon_coords.groupby('tID')['gen_stop'].max()

calculated_start_coords = calculated_start_coords.sort_values().reset_index(drop=True)
calculated_end_coords = calculated_end_coords.sort_values().reset_index(drop=True)

given_start_coords = given_start_coords.sort_values().reset_index(drop=True)
given_end_coords = given_end_coords.sort_values().reset_index(drop=True)


# assert equlity of start coordinates
assert (calculated_start_coords ==
        given_start_coords).all()

# assert equality of end coordiantes
assert (calculated_end_coords == given_end_coords).all()
