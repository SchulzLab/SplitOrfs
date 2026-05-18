import pandas as pd


def load_so_results(so_results):
    predicted_so_orfs = pd.read_csv(so_results, header=0, sep='\t')
    so_transcripts = predicted_so_orfs['OrfTransID'].to_list()
    predicted_so_orfs['OrfPos'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: x.split(','))
    predicted_so_orfs['OrfStarts'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: [y.split('-')[0] for y in x])
    predicted_so_orfs['nr_SO_starts'] = predicted_so_orfs['OrfPos'].apply(
        lambda x: len(x))
    return predicted_so_orfs, so_transcripts


def load_so_categorization_df(so_categorization_df):
    so_trans_categorized_df = pd.read_csv(
        so_categorization_df, header=0, index_col=0)
    so_trans_categorized_df['OrfsWithDistinctURTrans'] = so_trans_categorized_df['OrfsWithDistinctURTrans'].apply(
        lambda x: x.split(',') if not isinstance(x, float) else [])
    return so_trans_categorized_df


def load_dna_ur_df(UR_path):
    dna_ur_df = pd.read_csv(UR_path, sep='\t', header=None, names=[
                            'chr', 'start', 'stop', 'ID', 'score', 'strand'])
    dna_ur_df['OrfID'] = dna_ur_df['ID'].str.split(':').apply(lambda x: x[1])
    dna_ur_df['OrfTransID'] = dna_ur_df['ID'].str.split(
        ':').apply(lambda x: x[0])

    nr_orfs_with_UR = len(dna_ur_df['OrfID'].unique())
    nr_transcripts_with_UR = len(dna_ur_df['OrfTransID'].unique())
    print('Number of ORFs with unique region', nr_orfs_with_UR)
    print('Number of transcripts with unique region', nr_transcripts_with_UR)

    return dna_ur_df, nr_orfs_with_UR, nr_transcripts_with_UR
