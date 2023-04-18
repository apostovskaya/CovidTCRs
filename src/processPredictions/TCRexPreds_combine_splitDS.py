import os
import pandas as pd
import datetime
from Functions import combine_tcrex

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine into one file TCRex predictions for every individual in split dataset
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder = './data/splitDS_IMSEQ'
# specify path to where TCRex output files are stored
directory = f'{folder}/TCRex_results'
# where to save processed TCRex output files
if not os.path.exists(f'{folder}/TCRex_processed'):
    os.makedirs(f'{folder}/TCRex_processed')
    print(f'Processed TCRex files will be saved in {folder}/TCRex_processed')
directory_new = f'{folder}/TCRex_processed'

# from CovidTCRS in PhD_2020
# list of all IDs of patients (folder names in the TCRex_results folder)
id_list = [one_id for one_id in os.listdir(directory)
           if os.path.isdir(os.path.join(directory, one_id))]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add seq counts and frequencies to tcrex output files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for sample in id_list:
    # tcrex output (predictions) file
    print(str(datetime.datetime.now()) +
          f':\nStarted combining TCRex predictions with '
          f'sequencing counts data for {sample}')
    file_tcrex = f'{directory}/{sample}/filtered_results.csv'
    tcrex = pd.read_csv(file_tcrex, sep='\t', header=6)
    tcrex['sample_id'] = sample

    # sequencing files with counts
    file_seq = f'{folder}/TCRex_format_meta_TRB/{sample}.txt'
    seq = pd.read_csv(file_seq, sep='\t', index_col=0)
    seq.drop(columns=['d', 'VEnd', 'DStart', 'DEnd', 'JStart'], inplace=True)
    # Total number of sequenced TCRs
    seq['total_count'] = seq['count'].sum()
    # Total number of different TCRs (V-CDR3-J)
    seq['unique_TCRs'] = len(
        seq[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene']].drop_duplicates()[
            'CDR3_beta'])
    # Total number of different CDR3s
    seq['unique_CDR3s'] = len(seq['CDR3_beta'].unique())

    # change V-J gene notation to IMGT (TCRex also has the same, IMGT)
    # remove white spaces V-J genes
    seq['TRBV_gene'] = seq['TRBV_gene'].apply(str).str.replace(' ', '')
    seq['TRBJ_gene'] = seq['TRBJ_gene'].apply(str).str.replace(' ', '')

    # Convert V-J genes to IMGT format
    # seq['TRBV_gene'] = seq['TRBV_gene'].str.replace('\*00.+$', '', regex=True)
    # change TRBV/J1- into TRBV/J01-
    seq.loc[seq['TRBV_gene'].str.contains(pat='TRBV[1-9]-',
                                          regex=True), 'TRBV_gene'] = \
        seq.loc[seq['TRBV_gene'].str.contains(pat='TRBV[1-9]-',
                                              regex=True), 'TRBV_gene'].str.replace(
            'TRBV', 'TRBV0',
            regex=False)
    seq.loc[seq['TRBJ_gene'].str.contains(pat='TRBJ[1-9]-',
                                          regex=True), 'TRBJ_gene'] = \
        seq.loc[seq['TRBJ_gene'].str.contains(pat='TRBJ[1-9]-',
                                              regex=True), 'TRBJ_gene'].str.replace(
            'TRBJ', 'TRBJ0',
            regex=False)
    # change TRBV/J1 into TRBV/J01 (no - )
    seq.loc[seq['TRBV_gene'].str.contains(pat='TRBV[1-9]$',
                                          regex=True), 'TRBV_gene'] = \
        seq.loc[seq['TRBV_gene'].str.contains(pat='TRBV[1-9]$',
                                              regex=True), 'TRBV_gene'].str.replace(
            'TRBV', 'TRBV0',
            regex=False)
    seq.loc[seq['TRBJ_gene'].str.contains(pat='TRBJ[1-9]$',
                                          regex=True), 'TRBJ_gene'] = \
        seq.loc[seq['TRBJ_gene'].str.contains(pat='TRBJ[1-9]$',
                                              regex=True), 'TRBJ_gene'].str.replace(
            'TRBJ', 'TRBJ0',
            regex=False)
    # change -1 into -01
    seq.loc[seq['TRBV_gene'].str.contains(pat='-[1-9]$',
                                          regex=True), 'TRBV_gene'] = \
        seq.loc[seq['TRBV_gene'].str.contains(pat='-[1-9]$',
                                              regex=True), 'TRBV_gene'].str.replace(
            '-', '-0',
            regex=False)
    seq.loc[seq['TRBJ_gene'].str.contains(pat='-[1-9]$',
                                          regex=True), 'TRBJ_gene'] = \
        seq.loc[seq['TRBJ_gene'].str.contains(pat='-[1-9]$',
                                              regex=True), 'TRBJ_gene'].str.replace(
            '-', '-0',
            regex=False)

    # combine tcrex and sequencing files
    tcrex_seq = pd.merge(tcrex, seq, how="left",
                         on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])
    tcrex_seq.drop_duplicates(inplace=True)

    # number of different nt sequences per CDR3 with the same V-J genes, epitope, etc.
    nt_seq = tcrex_seq.groupby(
        ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
         'score', 'bpr'],
        as_index=False, sort=False)['count'].count().rename(
        columns={'count': 'N_nt_seq'})
    # combined counts for CDR3 with the same V-J genes, epitope, etc. but different nt sequences
    new_count = tcrex_seq.groupby(
        ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
         'score', 'bpr'],
        as_index=False, sort=False)['count'].sum()
    # combined frequencies for CDR3 with the same V-J genes, epitope, etc. but different nt sequences
    new_freq = tcrex_seq.groupby(
        ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
         'score', 'bpr'],
        as_index=False, sort=False)['freq'].sum()

    # step-by-step merging to update count and freq columns and substitude cdr3nt with N_nt_seq
    interim = pd.merge(nt_seq, new_count, how='inner',
                       on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope',
                           'pathology',
                           'score', 'bpr'])
    interim = pd.merge(interim, new_freq, how='inner',
                       on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope',
                           'pathology',
                           'score', 'bpr'])
    tcrex_seq = pd.merge(tcrex_seq[
                             ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope',
                              'pathology', 'score', 'bpr',
                              'sample_id', 'total_count', 'unique_TCRs',
                              'unique_CDR3s']], interim, how='inner',
                         on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope',
                             'pathology', 'score', 'bpr'])

    tcrex_seq.drop_duplicates(inplace=True)
    tcrex_seq.sort_values(by=['bpr', 'pathology', 'epitope', 'CDR3_beta'],
                          inplace=True)
    tcrex_seq.to_csv(f'{directory}/{sample}/results_with_seq_counts.csv',
                     index=False)
    print(f"{sample} files have been processed!")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine individual sample files into one; add metadata; keep only viral TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# metadata file
meta = pd.read_csv(folder + '/metadata_TRB_10umi_filter_MiXCR_to_VDJtools.txt',
                   sep='\t')
meta = meta[
    ['sample_id', 'patient_id', 'T-cell_type', 'Repeat', 's_id', 'Gene',
     'study_patient_id', 'disease_severity', 'time']]

tcrex_data = combine_tcrex(folder, id_list, 'results_with_seq_counts.csv')
# All predicted epitopes: ['SARS-CoV-2', 'HIV', 'EBV', 'CMV', 'HCV', 'DENV1', 'Melanoma', 'Influenza',
# 'DENV2', 'DENV3/4', 'HTLV1', 'HSV2', 'Multiple Myeloma', 'YellowFeverVirus']

tcrex_viral = tcrex_data.loc[~((tcrex_data.pathology == 'Melanoma') | (
            tcrex_data.pathology == 'Multiple Myeloma'))]
tcrex_viral = tcrex_viral.sort_values(
    by=['sample_id', 'bpr', 'pathology', 'epitope', 'CDR3_beta'])

# add metadata
tcrex_viral_meta = pd.merge(tcrex_viral, meta, on='sample_id')
tcrex_viral_meta = tcrex_viral_meta.sort_values(
    by=['study_patient_id', 'patient_id', 's_id', 'Repeat',
        'bpr', 'pathology', 'epitope', 'CDR3_beta'])

# convert time from d10 (str) to day: 10 (int)
tcrex_viral_meta['time'] = tcrex_viral_meta['time'].str.replace(pat='d',
                                                                repl='',
                                                                regex=False).map(
    lambda x: int(x))
tcrex_viral_meta.rename(columns={'time': 'day'}, inplace=True)

# creating week column
tcrex_viral_meta.loc[tcrex_viral_meta['day'] <= 7, 'week'] = 1
tcrex_viral_meta.loc[(tcrex_viral_meta['day'] > 7) & (
            tcrex_viral_meta['day'] <= 14), 'week'] = 2
tcrex_viral_meta.loc[(tcrex_viral_meta['day'] > 14) & (
            tcrex_viral_meta['day'] <= 21), 'week'] = 3
tcrex_viral_meta.loc[(tcrex_viral_meta['day'] > 21) & (
            tcrex_viral_meta['day'] <= 28), 'week'] = 4
tcrex_viral_meta['week'] = tcrex_viral_meta['week'].map(lambda x: int(x))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MERGING (combined or intersected) REPEATS TOGETHER; PATIENTS (TP and cell type) are still SEPARATE
# ? update frequencies to the median ones when merging (adj.Freq and adj.Total_counts as an alternative)
# ? add CDR3 diversity column
# ? add N_unique V, J genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new patient id will contain patient_id and cell type only (united/intersected repeats -> no s_id)
new_pat_id = []
for patient in tcrex_viral_meta['patient_id'].unique().tolist():
    for cell in ['cd4', 'cd8']:
        new_pat_id.append(f'{str(patient)}_{cell}')

# combine 3 sequencing repeats into one TCRex output file per patient
# TP and cell type are still separate samples
all_combined = pd.DataFrame()
all_intersect2 = pd.DataFrame()  # occurs at least twice, in two repeats
for new_patient in new_pat_id:
    # to drop columns 'sample_id', 'Repeat', 's_id', 'Gene'
    df1 = tcrex_viral_meta.loc[
        tcrex_viral_meta['sample_id'].str.startswith(pat=(new_patient + '_1')),
        ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
         'score', 'bpr',
         'total_count', 'unique_TCRs', 'unique_CDR3s', 'N_nt_seq', 'count',
         'freq',
         'patient_id', 'T-cell_type', 'study_patient_id', 'disease_severity',
         'day', 'week']]
    df2 = tcrex_viral_meta.loc[
        tcrex_viral_meta['sample_id'].str.startswith(pat=(new_patient + '_2')),
        ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
         'score', 'bpr',
         'total_count', 'unique_TCRs', 'unique_CDR3s', 'N_nt_seq', 'count',
         'freq',
         'patient_id', 'T-cell_type', 'study_patient_id', 'disease_severity',
         'day', 'week']]
    df3 = tcrex_viral_meta.loc[
        tcrex_viral_meta['sample_id'].str.startswith(pat=(new_patient + '_3')),
        ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
         'score', 'bpr',
         'total_count', 'unique_TCRs', 'unique_CDR3s', 'N_nt_seq', 'count',
         'freq',
         'patient_id', 'T-cell_type', 'study_patient_id', 'disease_severity',
         'day', 'week']]
    # merge all entries of all three repeats (filter later)
    df_combined = pd.concat([df1, df2, df3], ignore_index=True, sort=False)
    df_combined['new_id'] = new_patient

    # keep CDR3s that occur at least in 2 repeats out of 3
    # (many same CDR3s have different predicted specificity and V/J genes)
    # all repeated CDR3s
    # confusing keep parameter, but every repeat gets assigned True
    df_intersect2 = df_combined[
        df_combined.duplicated(subset=['CDR3_beta'], keep=False)].copy()
    df_intersect2.sort_values('CDR3_beta', inplace=True)
    try:
        # V-J-CDR3 that occur at least in two repeats
        df_intersect2 = pd.concat(g for _, g in df_intersect2.groupby(
            ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])
                                  if len(g['total_count'].unique()) > 1)
        df_intersect2.sort_values(
            by=['total_count', 'bpr', 'pathology', 'CDR3_beta'],
            inplace=True)  # ascending
        # keep only the last one, with the highest counts (gets assigned False)
        df_intersect2 = df_intersect2[~df_intersect2.duplicated(
            subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene',
                    'epitope'], keep='last')]
        df_intersect2.sort_values(by=['bpr', 'pathology', 'CDR3_beta'],
                                  inplace=True)
    except ValueError as err:
        print(
            f"{err}: {new_patient} does not have sequences with the same V-J-CDR3 occurring in at least two repeats. "
            f"Looking for the same CDR3 sequences only (disregarding V-J genes.)")
        try:
            # CDR3 that occur at least in two repeats
            df_intersect2 = pd.concat(g for _, g in df_intersect2.groupby(
               'CDR3_beta') if len(g['total_count'].unique()) > 1)
            df_intersect2.sort_values(
                by=['total_count', 'bpr', 'pathology', 'CDR3_beta'],
                inplace=True)  # ascending
            # keep only the last one, with the highest counts (gets assigned False)
            df_intersect2 = df_intersect2[~df_intersect2.duplicated(
                subset=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene',
                        'epitope'], keep='last')]
            df_intersect2.sort_values(by=['bpr', 'pathology', 'CDR3_beta'],
                                      inplace=True)
        except ValueError as err2:
            print(
                f"{err2}: {new_patient} does not have sequences with the same V-CDR3-J or just CDR3 occurring "
                f"in at least two repeats. Skipping to the next file.")
        else:
            print(
                f"File with the same CDR3 only occurring in at least two repeats has been created successfully "
                f"for {new_patient}.")
    else:
        print(
            f"File with the same V-CDR3-J occurring in at least two repeats has been created successfully "
            f"for {new_patient}.")
    # assemble all patients together
    all_intersect2 = all_intersect2.append(df_intersect2, ignore_index=True,
                                           sort=False)
all_intersect2.to_csv(f'{directory_new}/viralTCRexPreds_withCounts_splitDS_intersected2x_meta.tsv',
                      sep='\t', index=False)
# all_intersected_double_viral_meta
