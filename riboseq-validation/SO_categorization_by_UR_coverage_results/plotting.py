import os
import os.path
import matplotlib.pyplot as plt
import plotly.express as px


def plot_three_category_pie(cat1, cat2, cat3, nr_total_so, names_list, title, out_dir, figname, region_type, color_order):
    # plot pie chart with numbers: first, middle, last ORF
    fig, ax = plt.subplots()
    ax.pie([cat1, cat2, cat3],
           labels=names_list,
           autopct=lambda p: '{:.0f}'.format(p * nr_total_so / 100),
           colors=color_order)
    plt.title(f'{title} - {region_type}')

    plt.savefig(
        os.path.join(out_dir, f'{region_type}_{figname}.svg'), bbox_inches='tight')
    plt.close()


def plot_val_so_sets(nr_orfs_with_UR, nr_validated_so, total_nr_so, outdir, region_type, category_names=['# validated SO', '# SO with UR not validated',
                                                                                                         '# SO without UR']):
    nr_so_not_validated = nr_orfs_with_UR - nr_validated_so
    nr_so_no_UR = total_nr_so - nr_orfs_with_UR

    assert nr_validated_so + nr_so_not_validated + nr_so_no_UR == total_nr_so

    plot_three_category_pie(nr_validated_so,
                            nr_so_not_validated,
                            nr_so_no_UR,
                            total_nr_so,
                            category_names,
                            'Split-ORF validation pie chart',
                            outdir,
                            'SO_validation_pie_chart',
                            region_type,
                            ['#CC79A7', '#FFC500', '#75C1C5']
                            )


def plot_sunburst_ribo_cov_orf(so_categorization_df, outdir, region_type, sample_type):
    '''
    sunburst plots of Split-ORF transcript categorization by unique region presence
    and coverage
    '''
    df_agg = (
        so_categorization_df.groupby(['UR present', 'expressed', 'UR covered'])
        .size()
        .reset_index(name='count')
    )
    fig = px.sunburst(
        df_agg, path=['UR present', 'expressed', 'UR covered'], values='count')
    fig.write_image(os.path.join(
        outdir, f'sunburst_expression_info_{region_type}_{sample_type}.svg'))

    df_agg = (
        so_categorization_df.groupby(['UR present', 'UR covered'])
        .size()
        .reset_index(name='count')
    )
    fig = px.sunburst(
        df_agg, path=['UR present', 'UR covered'], values='count')
    fig.write_image(os.path.join(
        outdir, f'sunburst_no_expression_info_{region_type}_{sample_type}.svg'))

    so_categorization_df_filtered = so_categorization_df[
        so_categorization_df['expressed'] == 'ribocov gene']
    df_agg = (
        so_categorization_df_filtered.groupby(['UR present', 'UR covered'])
        .size()
        .reset_index(name='count')
    )
    fig = px.sunburst(
        df_agg, path=['UR present', 'UR covered'], values='count')
    fig.write_image(os.path.join(
        outdir, f'sunburst_expression_filtered_{region_type}_{sample_type}.svg'))

    so_categorization_df_filtered_2 = so_categorization_df[
        so_categorization_df['UR present'] != 'no UR']
    df_agg = (
        so_categorization_df_filtered_2.groupby(
            ['UR present', 'expressed', 'UR covered'])
        .size()
        .reset_index(name="count")
    )

    # color_map = {
    #     "no UR cov": "#0f6fb4",
    #     "cov UR in first ORF": "#d66b9b",
    #     "cov Ur in first and second ORF": "#dd3a55",
    #     "cov UR in second ORF": "#da7a2b",
    # }

    fig = px.sunburst(df_agg, path=[
                      'UR present', 'expressed',  'UR covered'], values="count", color='UR present')
    fig.write_image(os.path.join(
        outdir, f'sunburst_UR_present_filtered_{region_type}_{sample_type}.svg'))

    so_categorization_df_filtered = so_categorization_df[
        (so_categorization_df['expressed'] == 'ribocov gene') & (so_categorization_df['UR present'] != 'no UR')]
    df_agg = (
        so_categorization_df_filtered.groupby(['UR present', 'UR covered'])
        .size()
        .reset_index(name='count')
    )
    fig = px.sunburst(
        df_agg, path=['UR present', 'UR covered'], values='count', color='UR covered')
    fig.write_image(os.path.join(
        outdir, f'sunburst_expression_and_UR_present_filtered_{region_type}_{sample_type}.svg'))


def plot_possible_ribocov_information_pie(so_categorization_df, outdir, region_type, sample_type):
    so_categorization_df['ribocov_analyzable'] = 'ribocov analyzable'
    so_categorization_df.loc[so_categorization_df['expressed'] ==
                             'not ribocov gene', 'ribocov_analyzable'] = 'not ribocov gene'
    so_categorization_df.loc[so_categorization_df['UR present']
                             == 'no UR', 'ribocov_analyzable'] = 'no UR'

    nr_ribocov_analyzable = sum(
        so_categorization_df['ribocov_analyzable'] == 'ribocov analyzable')
    nr_not_ribocov = sum(
        so_categorization_df['ribocov_analyzable'] == 'not ribocov gene')
    nr_no_ur_ribocov = sum(
        so_categorization_df['ribocov_analyzable'] == 'no UR')

    total_nr_so_trans = len(so_categorization_df.index)

    assert nr_ribocov_analyzable + nr_not_ribocov + \
        nr_no_ur_ribocov == total_nr_so_trans
    category_names = ['ribocov analyzable', 'not ribocov gene', 'no UR']

    plot_three_category_pie(nr_ribocov_analyzable,
                            nr_not_ribocov,
                            nr_no_ur_ribocov,
                            total_nr_so_trans,
                            category_names,
                            'Split-ORF ribocov pie chart',
                            outdir,
                            'SO_ribocov_pie_chart',
                            f'{region_type}_{sample_type}',
                            ['#CC79A7', '#FFC500', '#75C1C5']
                            )
