import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def combine_tcrex(path, sample_list, file_in, tcrex_in='TCRex_results',
                  delim_in=',', head='infer'):
    """
    For (un)processed TCRex_results, loads files one by one and concatenates them into one combined df
    (which is saved as a new file)
    :param path: str, path to a parent folder (not ending with /)
    :param sample_list: list, all sample names (the same as folder names where TCRex results are stored);
    all files should have the same column names
    :param file_in: str, name of the input files (the same for all samples, folder names are sample-specific)
    :param tcrex_in: str, folder name within the parent directory with TCRex results per sample
    (just a name, no / at the end, ref. to default value)
    :param delim_out: str, delimeter in the output files (as specified in sep argument of pd.read_csv)
    :param delim_in: str, delimeter in the input files (as specified in sep argument of pd.read_csv)
    :param head: int, length of the file header (the number of lines in the file before column names)
    :return: df, dataframe (and saved file as tsv) with TCRex results combined from different samples
    """
    # create empty dataframe where to combine all separate files to
    col_names = pd.read_csv(f'{path}/{tcrex_in}/' + sample_list[0] +
                            f'/{file_in}',
                            sep=delim_in, header=head).columns.to_list()
    # column names are the same as in the individual dataframes after processing (tcrex_seq files)
    # from tcrex_seq.columns.to_list()
    # ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr', 'sample_id', '
    # total_count', 'N_nt_seq', 'count', 'freq']
    df = pd.DataFrame(columns=col_names)

    for one_sample in sample_list:
        print(f'Processing file {one_sample}:')
        file_name = f'{path}/{tcrex_in}/{one_sample}/{file_in}'
        one_tcrex = pd.read_csv(file_name, sep=delim_in, header=head)
        assert all(elem in one_tcrex.columns.to_list() for elem in col_names), \
            f"The file {one_sample} doesn't contain required columns, please, check the spelling and/or case. " \
            f"The following columns are expected: {col_names} based on the first input file {sample_list[0]}."
        one_tcrex = one_tcrex[col_names]
        df = df.append(one_tcrex, ignore_index=True, sort=False)
        print(f"{one_sample} has been added!")
    return df


def tcr_counter_split_ds(all_data, one_id):
    """

    :param all_data: df, a data frame with all TCRex predictions for all patients, cell, types, time points
    :param one_id: str, id of one  patient
    :return: df
    """
    # tmp - data of one type of T cells of one patient
    tmp = all_data[all_data['new_id'] == one_id].copy()
    if tmp.empty:
        print(f"{one_id} has no TCRs that occured at least in two repeats and fulfilled the filtering criteria")
        return

    # drop unnecessary columns
    tmp.drop(columns=['score', 'bpr', 'N_nt_seq'], inplace=True)

    # choosing mode (the most frequent) values from the repeats (so basically the repeat that has the majority of TCRs)
    # or just keep the highest since I'm keeping the freqs from the highest?
    tmp['total_count'] = tmp['total_count'].mode()[0]
    tmp['unique_TCRs'] = tmp['unique_TCRs'].mode()[0]
    tmp['unique_CDR3s'] = tmp['unique_CDR3s'].mode()[0]

    # renew column orders
    tmp = tmp[['new_id', 'patient_id', 'T-cell_type', 'study_patient_id', 'disease_severity', 'day', 'week',
               'total_count', 'unique_TCRs', 'unique_CDR3s', 'species',
               'count', 'freq', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology']]

    # number of unique cdr3-epitope pairs
    # tmp['Total_TCRs'] = len(tmp[['CDR3_beta', 'epitope']].drop_duplicates()['CDR3_beta'].unique())

    # Epitopes
    # number of different epitope groups: unique to SARS-CoV-2 / common for other coronaviruses / cr with other viruses
    tmp['N_SARS-CoV-2_Eps'] = len(tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'epitope'].unique())
    tmp['N_SARS-CoV-2_only_Eps'] = len(tmp.loc[(tmp['pathology'] == 'SARS-CoV-2') & (tmp['species'] == 1),
                                               'epitope'].unique())
    tmp['N_CoV_Eps'] = len(tmp.loc[(tmp['pathology'] == 'SARS-CoV-2') & (tmp['species'] > 1), 'epitope'].unique())
    tmp['Frac_CoV_in_SARS-CoV-2_Eps'] = tmp['N_CoV_Eps'] / tmp['N_SARS-CoV-2_Eps']
    tmp['Frac_SARS-CoV-2_only_Eps'] = tmp['N_SARS-CoV-2_only_Eps'] / tmp['N_SARS-CoV-2_Eps']
    tmp['N_viral_Eps'] = len(tmp.loc[tmp['pathology'] != 'SARS-CoV-2', 'epitope'].unique())

    # new addition, works correctly
    ep_list = tmp['epitope'].unique().tolist()
    covid_ep_list = tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'epitope'].unique().tolist()
    ep_list.sort()
    covid_ep_list.sort()

    # TCRs
    #
    # Covid-specific TCRs
    #
    tmp['N_SARS-CoV-2_TCRs'] = tmp[tmp['pathology'] ==
                                   'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                         'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_SARS-CoV-2_TCRs'] = tmp[tmp['pathology'] ==
                                      'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                            'TRBJ_gene'])['freq'].sum()
    #
    # Covid-specific ONLY TCRs (NOT recognizing other coronaviruses)
    #
    tmp['N_SARS-CoV-2_only_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                                        (tmp['species'] == 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                       'TRBJ_gene'])[
        'CDR3_beta'].count()
    # fraction out of the entire repertoire
    tmp['Frac_SARS-CoV-2_only_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                                           (tmp['species'] == 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                          'TRBJ_gene'])['freq'].sum()
    #
    # CoV TCRs: Covid-specific TCRs recognizing epitopes that occur in other coronaviruses
    #
    tmp['N_CoV_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                            (tmp['species'] > 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                          'TRBJ_gene'])['CDR3_beta'].count()
    # fraction out of the entire repertoire
    tmp['Frac_CoV_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                               (tmp['species'] > 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                             'TRBJ_gene'])['freq'].sum()

    # fraction out of the SARS-CoV-2-specific repertoire
    if tmp['Frac_SARS-CoV-2_TCRs'].all() > 0:
        tmp['Frac_CoV_in_SARS-CoV-2_TCRs'] = tmp['Frac_CoV_TCRs'] / tmp['Frac_SARS-CoV-2_TCRs']
    else:
        print(f"{one_id}: no SARS-CoV-2_TCRs")
        # writer.write(f"{one_id}: no SARS-CoV-2_TCRs")
        tmp['Frac_CoV_in_SARS-CoV-2_TCRs'] = 0

    #
    # TCRs specific to other viruses
    #
    tmp['N_viral_TCRs'] = tmp[tmp['pathology'] !=
                              'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                    'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_viral_TCRs'] = tmp[tmp['pathology'] !=
                                 'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                       'TRBJ_gene'])['freq'].sum()
    # CROSS-REACTIVITY
    #
    # SARS-CoV-2 crTCRs: tcrs that recognize several covid epitopes
    #
    tmp['N_SARS-CoV-2_crTCRs'] = len(tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'CDR3_beta']) - \
                                 tmp[tmp['pathology'] == 'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                               'TRBJ_gene'])[
                                     'CDR3_beta'].count()

    # fraction of SARS-CoV-2 crTCRs out of the entire repertoire
    subset = tmp[tmp['pathology'] == 'SARS-CoV-2']
    tmp['Frac_SARS-CoV-2_crTCRs'] = subset.loc[subset.duplicated(subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'],
                                                                 keep='last'), 'freq'].sum()

    # fraction of SARS-CoV-2 crTCRs out of the SARS-CoV-2-specific repertoire
    if tmp['Frac_SARS-CoV-2_TCRs'].all() > 0:
        tmp['Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs'] = tmp['Frac_SARS-CoV-2_crTCRs'] / tmp['Frac_SARS-CoV-2_TCRs']
    else:
        print(f"{one_id}: no SARS-CoV-2_TCRs")
        # writer.write(f"{one_id}: no SARS-CoV-2_TCRs")
        tmp['Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs'] = 0

        # fraction of CoV TCRs out of SARS-CoV-2 crTCRs
    subset = subset[subset.duplicated(subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'], keep='last')]
    if tmp['Frac_SARS-CoV-2_crTCRs'].all() > 0:
        tmp['Frac_CoV_in_crTCRs'] = subset.loc[subset['species'] > 1, 'freq'].sum() / tmp['Frac_SARS-CoV-2_crTCRs']
    else:
        print(f"{one_id}: no TCRs recognizing multiple SARS-CoV-2 epitopes")
        # writer.write(f"{one_id}: no TCRs recognizing multiple SARS-CoV-2 epitopes")
        tmp['Frac_CoV_in_crTCRs'] = 0
    del subset

    #
    # viral_crTCRs: covid tcrs cross-reactive with other viruses
    #
    # used to be .count() instead of .nunique() -> suspect that duplicates of pathogens were counted
    tmp2 = tmp.groupby('CDR3_beta', as_index=False)['pathology'].nunique()
    tmp2.rename(columns={'pathology': 'cr'}, inplace=True)
    tmp = pd.merge(tmp, tmp2, on='CDR3_beta')
    virus_list = tmp.loc[tmp['cr'] > 1, 'pathology'].unique().tolist()
    virus_list.sort()
    # fraction of viral_crTCRs out of the entire repertoire where is drop duplicated????
    # tmp['N_SARS-CoV-2_viral_crTCRs'] = tmp.loc[(tmp['cr'] > 1) & (tmp['pathology'] == 'SARS-CoV-2'),
    #                                            'CDR3_beta'].count()
    # tmp['N_SARS-CoV-2_viral_crTCRs'] = tmp.loc[(tmp['cr'] > 1) & (tmp['pathology'] == 'SARS-CoV-2'),
    #                                            'freq'].sum()

    tmp['N_SARS-CoV-2_viral_crTCRs'] = tmp[(tmp['cr'] > 1) &
                                           (tmp['pathology'] ==
                                            'SARS-CoV-2')].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                   'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_SARS-CoV-2_viral_crTCRs'] = tmp[(tmp['cr'] > 1) &
                                              (tmp['pathology'] ==
                                               'SARS-CoV-2')].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                      'TRBJ_gene'])['freq'].sum()

    # fraction of viral_crTCRs out of the SARS-CoV-2-specific repertoire WRONG!
    # tmp['Frac_viral_crTCRs_in_SARS-CoV-2'] = tmp[(tmp['cr'] > 1) &
    #                                              (tmp['pathology'] ==
    #                                               'SARS-CoV-2')].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
    #                                                                                      'TRBJ_gene'])['freq'].sum() /
    #                                          tmp['Frac_SARS-CoV-2_TCRs']

    # fraction of CoV TCRs out of viral_crTCRs CHECK!
    # instead of try-except ZeroDivisionError - pandas handles it to NaN and doesn't through exception
    if tmp['Frac_SARS-CoV-2_viral_crTCRs'].all() > 0:
        tmp['Frac_CoV_in_viral_crTCRs'] = tmp[(tmp['cr'] > 1) &
                                              (tmp['species'] > 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                            'TRBJ_gene'])[
                                              'freq'].sum() / \
                                          tmp['Frac_SARS-CoV-2_viral_crTCRs']
    else:
        print(f"{one_id}: no TCRs recognizing SARS-CoV-2 and other viral epitopes")
        # doesn't write to a file ;(
        # print(f"{one_id}: no TCRs recognizing SARS-CoV-2 and other viral epitopes", file=writer)
        tmp['Frac_CoV_in_viral_crTCRs'] = 0

    # calculate diff diversity measures
    # from skbio.diversity.alpha import shannon

    # delete columns
    tmp.drop(columns=['count', 'freq', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'species', 'cr'],
             inplace=True)
    tmp.drop_duplicates(inplace=True)

    # add columns listing viruses and epitopes recognized
    tmp['cr_viruses'] = [virus_list]
    tmp['all_epitopes'] = [ep_list]
    tmp['all_SARS-CoV-2_epitopes'] = [covid_ep_list]

    return tmp


def tcr_counter_mixed_ds(all_data, one_id):
    """

    :param all_data: df, a data frame with all TCRex predictions for all patients, cell, types, time points
    :param one_id: str
    :return: df IR-Binder
    """
    # tmp - data of one patient
    tmp = all_data[all_data['patient_id'] == one_id].copy()
    if tmp.empty:
        print("\n" + f"{one_id} has no TCRs with predicted specificity" + "\n")
        tmp.drop(columns=['count', 'freq', 'score', 'bpr', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope',
                          'pathology', 'species'], inplace=True)
        tmp.drop_duplicates(inplace=True)
        return tmp

    # drop unnecessary columns
    tmp.drop(columns=['score', 'bpr', 'N_nt_seq'], inplace=True)

    # substitute 0 with nan in counts and freqs and then nan with min values
    tmp['count'] = tmp['count'].replace(to_replace=0.0, value=np.nan)
    tmp['freq'] = tmp['freq'].replace(to_replace=0.0, value=np.nan)
    min_cnt = tmp['count'].min()
    min_freq = tmp['freq'].min()
    values = {'count': min_cnt, 'freq': min_freq}
    tmp.fillna(value=values, inplace=True)
    # filter out low freq clones
    before = tmp.shape[0]
    tmp = tmp[tmp['freq'] > 1 / 100000]
    print(f'{one_id}: {before - tmp.shape[0]} clones were removed due to low frequency (< 1/100000), '
          f'{tmp.shape[0]} clones remained')

    # Epitopes
    # number of different epitope groups: unique to SARS-CoV-2 / common for other coronaviruses
    tmp['N_SARS-CoV-2_Eps'] = len(tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'epitope'].unique())
    tmp['N_SARS-CoV-2_only_Eps'] = len(tmp.loc[(tmp['pathology'] == 'SARS-CoV-2') & (tmp['species'] == 1),
                                               'epitope'].unique())
    tmp['N_CoV_Eps'] = len(tmp.loc[(tmp['pathology'] == 'SARS-CoV-2') & (tmp['species'] > 1), 'epitope'].unique())

    if tmp['N_SARS-CoV-2_Eps'].all() > 0:
        tmp['Frac_CoV_in_SARS-CoV-2_Eps'] = tmp['N_CoV_Eps'] / tmp['N_SARS-CoV-2_Eps']
        tmp['Frac_SARS-CoV-2_only_Eps'] = tmp['N_SARS-CoV-2_only_Eps'] / tmp['N_SARS-CoV-2_Eps']
    else:
        print(f"{one_id}: no SARS-CoV-2_Eps")
        tmp['Frac_CoV_in_SARS-CoV-2_Eps'] = 0
        tmp['Frac_SARS-CoV-2_only_Eps'] = 0
    # tmp['N_viral_Eps'] = len(tmp.loc[tmp['pathology'] != 'SARS-CoV-2', 'epitope'].unique())

    # new addition, works correctly
    ep_list = tmp['epitope'].unique().tolist()
    covid_ep_list = tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'epitope'].unique().tolist()
    ep_list.sort()
    covid_ep_list.sort()

    # TCRs
    #
    # Covid-specific TCRs
    #
    tmp['N_SARS-CoV-2_TCRs'] = tmp[tmp['pathology'] ==
                                   'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                         'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_SARS-CoV-2_TCRs'] = tmp[tmp['pathology'] ==
                                      'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                            'TRBJ_gene'])['freq'].sum()
    #
    # Covid-specific ONLY TCRs (NOT recognizing other coronaviruses)
    #
    tmp['N_SARS-CoV-2_only_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                                        (tmp['species'] == 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                       'TRBJ_gene'])[
        'CDR3_beta'].count()
    # fraction out of the entire repertoire
    tmp['Frac_SARS-CoV-2_only_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                                           (tmp['species'] == 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                          'TRBJ_gene'])['freq'].sum()
    #
    # CoV TCRs: Covid-specific TCRs recognizing epitopes that occur in other coronaviruses
    #
    tmp['N_CoV_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                            (tmp['species'] > 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                          'TRBJ_gene'])['CDR3_beta'].count()
    # fraction out of the entire repertoire
    tmp['Frac_CoV_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                               (tmp['species'] > 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                             'TRBJ_gene'])['freq'].sum()
    # fraction out of the SARS-CoV-2-specific repertoire
    if tmp['Frac_SARS-CoV-2_TCRs'].all() > 0:
        tmp['Frac_CoV_in_SARS-CoV-2_TCRs'] = tmp['Frac_CoV_TCRs'] / tmp['Frac_SARS-CoV-2_TCRs']
    else:
        print(f"{one_id}: no SARS-CoV-2_TCRs")
        tmp['Frac_CoV_in_SARS-CoV-2_TCRs'] = 0

    #
    # TCRs specific to other viruses
    #
    tmp['N_viral_TCRs'] = tmp[tmp['pathology'] !=
                              'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                    'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_viral_TCRs'] = tmp[tmp['pathology'] !=
                                 'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                       'TRBJ_gene'])['freq'].sum()
    # CROSS-REACTIVITY
    #
    # SARS-CoV-2 crTCRs: tcrs that recognize several covid epitopes
    #
    tmp['N_SARS-CoV-2_crTCRs'] = len(tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'CDR3_beta']) - \
                                 tmp[tmp['pathology'] == 'SARS-CoV-2'].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                               'TRBJ_gene'])[
                                     'CDR3_beta'].count()
    # fraction of SARS-CoV-2 crTCRs out of the entire repertoire
    subset = tmp[tmp['pathology'] == 'SARS-CoV-2']
    tmp['Frac_SARS-CoV-2_crTCRs'] = subset.loc[subset.duplicated(subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'],
                                                                 keep='last'), 'freq'].sum()
    # fraction of SARS-CoV-2 crTCRs out of the SARS-CoV-2-specific repertoire
    if tmp['Frac_SARS-CoV-2_TCRs'].all() > 0:
        tmp['Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs'] = tmp['Frac_SARS-CoV-2_crTCRs'] / tmp['Frac_SARS-CoV-2_TCRs']
    else:
        print(f"{one_id}: no SARS-CoV-2_TCRs")
        # writer.write(f"{one_id}: no SARS-CoV-2_TCRs")
        tmp['Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs'] = 0

        # fraction of CoV TCRs out of SARS-CoV-2 crTCRs
    subset = subset[subset.duplicated(subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'], keep='last')]
    if tmp['Frac_SARS-CoV-2_crTCRs'].all() > 0:
        tmp['Frac_CoV_in_crTCRs'] = subset.loc[subset['species'] > 1, 'freq'].sum() / tmp['Frac_SARS-CoV-2_crTCRs']
    else:
        print(f"{one_id}: no TCRs recognizing multiple SARS-CoV-2 epitopes")
        tmp['Frac_CoV_in_crTCRs'] = 0
    del subset

    #
    # viral_crTCRs: covid tcrs cross-reactive with other viruses
    #
    # used to be .count() instead of .nunique() -> suspect that duplicates of pathogens were counted
    tmp2 = tmp.groupby('CDR3_beta', as_index=False)['pathology'].nunique()
    tmp2.rename(columns={'pathology': 'cr'}, inplace=True)
    tmp = pd.merge(tmp, tmp2, on='CDR3_beta')
    virus_list = tmp.loc[tmp['cr'] > 1, 'pathology'].unique().tolist()
    virus_list.sort()
    # fraction of viral_crTCRs out of the entire repertoire where is drop duplicated????
    # tmp['N_SARS-CoV-2_viral_crTCRs'] = tmp.loc[(tmp['cr'] > 1) & (tmp['pathology'] == 'SARS-CoV-2'),
    #                                            'CDR3_beta'].count()
    # tmp['N_SARS-CoV-2_viral_crTCRs'] = tmp.loc[(tmp['cr'] > 1) & (tmp['pathology'] == 'SARS-CoV-2'),
    #                                            'freq'].sum()

    tmp['N_SARS-CoV-2_viral_crTCRs'] = tmp[(tmp['cr'] > 1) &
                                           (tmp['pathology'] ==
                                            'SARS-CoV-2')].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                   'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_SARS-CoV-2_viral_crTCRs'] = tmp[(tmp['cr'] > 1) &
                                              (tmp['pathology'] ==
                                               'SARS-CoV-2')].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                      'TRBJ_gene'])['freq'].sum()

    # fraction of viral_crTCRs out of the SARS-CoV-2-specific repertoire WRONG!
    # tmp['Frac_viral_crTCRs_in_SARS-CoV-2'] = tmp[(tmp['cr'] > 1) &
    #                                              (tmp['pathology'] ==
    #                                               'SARS-CoV-2')].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
    #                                                                                      'TRBJ_gene'])['freq'].sum() /
    #                                          tmp['Frac_SARS-CoV-2_TCRs']

    # fraction of CoV TCRs out of viral_crTCRs CHECK!
    # instead of try-except ZeroDivisionError - pandas handles it to NaN and doesn't through exception
    if tmp['Frac_SARS-CoV-2_viral_crTCRs'].all() > 0:
        tmp['Frac_CoV_in_viral_crTCRs'] = tmp[(tmp['cr'] > 1) &
                                              (tmp['species'] > 1)].drop_duplicates(subset=['TRBV_gene', 'CDR3_beta',
                                                                                            'TRBJ_gene'])[
                                              'freq'].sum() / \
                                          tmp['Frac_SARS-CoV-2_viral_crTCRs']
    else:
        print(f"{one_id}: no TCRs recognizing SARS-CoV-2 and other viral epitopes")
        # doesn't write to a file ;(
        # print(f"{one_id}: no TCRs recognizing SARS-CoV-2 and other viral epitopes", file=writer)
        tmp['Frac_CoV_in_viral_crTCRs'] = 0

    # calculate diff diversity measures
    # from skbio.diversity.alpha import shannon

    # delete columns
    tmp.drop(columns=['count', 'freq', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'species', 'cr'],
             inplace=True)
    tmp.drop_duplicates(inplace=True)

    # add columns listing viruses and epitopes recognized
    tmp['cr_viruses'] = [virus_list]
    tmp['all_epitopes'] = [ep_list]
    tmp['all_SARS-CoV-2_epitopes'] = [covid_ep_list]

    return tmp


def tcr_counter_sc_ds(all_data, one_id):
    # tmp - data of one patient
    tmp = all_data[all_data['patient_id'] == one_id].copy()
    if tmp.empty:
        print("\n" + f"{one_id} has no TCRs with predicted specificity" + "\n")
        tmp.drop(
            columns=['cloneCount', 'cloneFreq', 'score', 'bpr',
                     'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope',
                     'pathology', 'species'], inplace=True)
        tmp.drop_duplicates(inplace=True)
        return tmp

    # drop unnecessary columns
    tmp.drop(columns=['score', 'bpr'], inplace=True)

    # substitute 0 with nan in counts and freqs and then nan with min values
    tmp['cloneCount'] = tmp['cloneCount'].replace(to_replace=0.0, value=np.nan)
    tmp['cloneFreq'] = tmp['cloneFreq'].replace(to_replace=0.0, value=np.nan)
    min_cnt = tmp['cloneCount'].min()
    min_freq = tmp['cloneFreq'].min()
    values = {'cloneCount': min_cnt, 'cloneFreq': min_freq}
    tmp.fillna(value=values, inplace=True)

    #
    # Epitopes
    #
    # number of different epitope groups: SC2-unique / CoV-common
    tmp['N_SARS-CoV-2_Eps'] = len(
        tmp.loc[tmp['pathology'] == 'SARS-CoV-2', 'epitope'].unique())
    tmp['N_SARS-CoV-2_only_Eps'] = len(
        tmp.loc[(tmp['pathology'] == 'SARS-CoV-2') & (tmp['species'] == 1),
                'epitope'].unique())
    tmp['N_CoV_Eps'] = len(tmp.loc[(tmp['pathology'] == 'SARS-CoV-2') & (
                tmp['species'] > 1), 'epitope'].unique())

    # all recognized epitopes
    covid_ep_list = tmp.loc[
        tmp['pathology'] == 'SARS-CoV-2', 'epitope'].unique().tolist()
    covid_ep_list.sort()

    # TCRs
    #
    # Covid-specific TCRs
    #
    tmp['N_SARS-CoV-2_TCRs'] = tmp[tmp['pathology'] ==
                                   'SARS-CoV-2'].drop_duplicates(
        subset=['TRBV_gene', 'CDR3_beta',
                'TRBJ_gene'])['CDR3_beta'].count()
    tmp['Frac_SARS-CoV-2_TCRs'] = tmp[tmp['pathology'] ==
                                      'SARS-CoV-2'].drop_duplicates(
        subset=['TRBV_gene', 'CDR3_beta',
                'TRBJ_gene'])['cloneFreq'].sum()
    #
    # Covid-specific ONLY TCRs (NOT recognizing other coronaviruses)
    #
    tmp['N_SARS-CoV-2_only_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                                        (tmp['species'] == 1)].drop_duplicates(
        subset=['TRBV_gene', 'CDR3_beta',
                'TRBJ_gene'])[
        'CDR3_beta'].count()
    # fraction out of the entire repertoire
    tmp['Frac_SARS-CoV-2_only_TCRs'] = \
        tmp[(tmp['pathology'] == 'SARS-CoV-2')
            & (tmp['species'] == 1)].drop_duplicates(
            subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])['cloneFreq'].sum()
    #
    # CoV-common TCRs
    #
    tmp['N_CoV_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                            (tmp['species'] > 1)].drop_duplicates(
        subset=['TRBV_gene', 'CDR3_beta',
                'TRBJ_gene'])['CDR3_beta'].count()
    # fraction out of the entire repertoire
    tmp['Frac_CoV_TCRs'] = tmp[(tmp['pathology'] == 'SARS-CoV-2') &
                               (tmp['species'] > 1)].drop_duplicates(
        subset=['TRBV_gene', 'CDR3_beta',
                'TRBJ_gene'])['cloneFreq'].sum()

    # delete columns
    tmp.drop(columns=['cloneCount', 'cloneFreq',
                      'TRBV_gene', 'CDR3_beta', 'TRBJ_gene',
                      'epitope', 'pathology', 'species'],
             inplace=True)
    tmp.drop_duplicates(inplace=True)

    # add columns listing viruses and epitopes recognized
    tmp['all_SARS-CoV-2_epitopes'] = [covid_ep_list]

    return tmp


# # PLOTS: Dynamics per patient and Trends per group
#
#
# def plot_person_dynamics(data_in, x_in, y_in, hue_order_in, hue_in, palette_in, style_in, saveas):
#     """Plots personal dynamics (connected dots)"""
#     plt.figure(figsize=(9, 7))
#     sns.lineplot(x=x_in, y=y_in, data=data_in, hue_order=hue_order_in, hue=hue_in,
#                  palette=palette_in, style=style_in, markers=True, dashes=False)
#     plt.ylabel(y_in)
#     plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=9)
#     plt.tight_layout()
#     plt.savefig(f'{saveas}/PersonDyn_{y_in}.png')
#     plt.close()
#
#
# def plot_dots_regression(data_in, x_in, y_in, hue_order_in, hue_in, palette_in, saveas):
#     """Plots dots with regression line"""
#     g = sns.lmplot(x=x_in, y=y_in, data=data_in, scatter_kws={'s': 10},
#                    hue_order=hue_order_in, hue=hue_in, palette=palette_in, ci=95, legend=False)
#     # g.set(ylim=(0, None))
#     plt.ylabel(y_in)
#     plt.legend(title=hue_in)
#     plt.tight_layout()
#     plt.savefig(f'{saveas}/RegressDots_{y_in}.png')
#     plt.close()
#
#
# def plot_group_trend(data_in, x_in, y_in, hue_order_in, hue_in, palette_in, saveas):
#     """Plots dots with trend line"""
#     plt.figure(figsize=(9, 7))
#     sns.lineplot(x=x_in, y=y_in, data=data_in, markers=True, dashes=False, style=hue_in,
#                  estimator=np.median, ci=95, n_boot=1000, seed=None, sort=True, err_style='band',
#                  hue_order=hue_order_in, hue=hue_in, palette=palette_in)
#     plt.ylabel(y_in)
#     plt.legend(title=hue_in)  # bbox_to_anchor=(1, 1) redundant to loc loc='upper right',
#     plt.tight_layout()
#     plt.savefig(f'{saveas}/TrendDots_{y_in}.png')
#     plt.close()
