import os
import re
import datetime
import pandas as pd
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For Public data, S-paper (IR-Binder)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder = '/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/IR-Binder_PRJEB38339/'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make metafile of MiXCR data for VDJtools
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta = {'file_name': [], 'sample_id': [], 's_id': [], 'type': [], 'patient_id': [], 'patient': [], 'Gene': []}
for f_name in os.listdir(folder + 'MiXCR/MiXCR_outputs'):
    if f_name.endswith('.TRB.txt'):
        if f_name.startswith('TRB-Pt'):
            # TRB-Pt-1-1-500ng-03-04-2020-gDNA_S27_mixcrOut.clonotypes.TRB or TRB-Pt-1-5_S94_mixcrOut.clonotypes.TRB
            meta['file_name'].append('./MiXCR/MiXCR_outputs/' + f_name)
            meta['sample_id'].append(f_name.split('_mixcr')[0])
            meta['s_id'].append(f_name.split('_')[1])
            meta['type'].append('Pt')
            match = re.search(r'Pt-[0-9]+-[0-9][\-\_]', f_name)[0][:-1]  # the last element is _ or -
            meta['patient_id'].append(match)
            match = match.split('Pt-')[1]
            meta['patient'].append(match)
            meta['Gene'].append(f_name.split('.')[2])
        else:
            # TCRb-HD39-PB-gDNA_S59_mixcrOut.clonotypes.TRB or TCRb-HD29_S78_mixcrOut.clonotypes.TRB
            # HD18-Mar2015-TCRb_S36_mixcrOut.clonotypes.TRB
            # TRB-HD01_S102_mixcrOut.clonotypes.TRB or TRB-HD04-09-03-2017-gDNA_S79_mixcrOut.clonotypes.TRB
            meta['file_name'].append('./MiXCR/MiXCR_outputs/' + f_name)
            meta['sample_id'].append(f_name.split('_mixcr')[0])
            meta['s_id'].append(f_name.split('_')[1])
            meta['type'].append('HD')
            match = re.search(r'HD[0-9][0-9]', f_name)[0]
            meta['patient_id'].append(match)
            match = match.split('HD')[1]
            meta['patient'].append(match)
            meta['Gene'].append(f_name.split('.')[2])
metadata = pd.DataFrame(meta)
metadata.sort_values(by=['patient_id'], inplace=True, ignore_index=True)
metadata.to_csv(path_or_buf=(folder + 'metadata_TRB_IR_Binder_MiXCR.txt'), sep='\t', index=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# transform MiXCR files to use without VDJtools (semiTCRex format: more columns and no predictions yet)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not os.path.exists(folder + 'MiXCR_parsed'):
    os.makedirs(folder + 'MiXCR_parsed')

# the number of files in the folder
total = 0
for f_name in os.listdir(folder + 'MiXCR/MiXCR_outputs'):
    if f_name.endswith('.TRB.txt'):
        total += 1

current = 0
error = 0
for f_name in os.listdir(folder + 'MiXCR/MiXCR_outputs'):
    if f_name.endswith('.TRB.txt'):
        current += 1
        if f_name.startswith('TRB-Pt'):
            new_name = re.search(r'Pt-[0-9]+-[0-9][\-\_]', f_name)[0][:-1]
            health = 'Pt'
            patient = new_name.split('Pt-')[1]
        else:
            new_name = re.search(r'HD[0-9][0-9]', f_name)[0]
            health = 'HD'
            patient = new_name.split('HD')[1]

        df = pd.read_csv(folder + f'MiXCR/MiXCR_outputs/{f_name}', sep='\t')
        print(str(datetime.datetime.now()) + ':\n' + f'Parsing MiXCR files of {new_name}...')
        # rename columns to VDJtools/TCRex style
        df.rename(columns={'aaSeqCDR3': 'CDR3_beta', 'bestVGene': 'TRBV_gene', 'bestJGene': 'TRBJ_gene',
                           'cloneCount': 'count', 'cloneFraction': 'freq'}, inplace=True)
        df = df[['count', 'freq', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'nSeqCDR3']]
        df['count'] = df['count'].astype('int32')
        df['total_count'] = df['count'].sum()
        df['unique_TCRs'] = len(df[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene']].drop_duplicates()['CDR3_beta'])
        df['unique_CDR3s'] = len(df['CDR3_beta'].unique())
        # remove rows with non-canonical CDR3 amino acid notation
        df = df.loc[~df['CDR3_beta'].str.contains(pat=r'[^A-Z]', regex=True)]
        # remove CDR3s not starting with C / not ending with F/W
        df = df.loc[df['CDR3_beta'].str.contains(pat=r'^C', regex=True)]
        df = df.loc[df['CDR3_beta'].str.contains(pat=r'F|W$', regex=True)]
        # adding patient id PMe compatible style
        df['patient'] = patient
        df['type'] = health
        df['study'] = 'IR-Binder'

        # change V-J gene notation to IMGT (TCRex also has the same, IMGT)
        # remove white spaces V-J genes
        df['TRBV_gene'] = df['TRBV_gene'].apply(str).str.replace(' ', '')
        df['TRBJ_gene'] = df['TRBJ_gene'].apply(str).str.replace(' ', '')

        # Convert V-J genes to IMGT format
        # df['TRBV_gene'] = df['TRBV_gene'].str.replace('\*00.+$', '', regex=True)
        # change TRBV/J1- into TRBV/J01-
        df.loc[df['TRBV_gene'].str.contains(pat='TRBV[1-9]-', regex=True), 'TRBV_gene'] = \
            df.loc[df['TRBV_gene'].str.contains(pat='TRBV[1-9]-',
                                                regex=True), 'TRBV_gene'].str.replace('TRBV', 'TRBV0', regex=False)
        df.loc[df['TRBJ_gene'].str.contains(pat='TRBJ[1-9]-', regex=True), 'TRBJ_gene'] = \
            df.loc[df['TRBJ_gene'].str.contains(pat='TRBJ[1-9]-',
                                                regex=True), 'TRBJ_gene'].str.replace('TRBJ', 'TRBJ0', regex=False)
        # change TRBV/J1 into TRBV/J01 (no - )
        df.loc[df['TRBV_gene'].str.contains(pat='TRBV[1-9]$', regex=True), 'TRBV_gene'] = \
            df.loc[df['TRBV_gene'].str.contains(pat='TRBV[1-9]$',
                                                regex=True), 'TRBV_gene'].str.replace('TRBV', 'TRBV0', regex=False)
        df.loc[df['TRBJ_gene'].str.contains(pat='TRBJ[1-9]$', regex=True), 'TRBJ_gene'] = \
            df.loc[df['TRBJ_gene'].str.contains(pat='TRBJ[1-9]$',
                                                regex=True), 'TRBJ_gene'].str.replace('TRBJ', 'TRBJ0', regex=False)
        # change -1 into -01
        df.loc[df['TRBV_gene'].str.contains(pat='-[1-9]$', regex=True), 'TRBV_gene'] = \
            df.loc[df['TRBV_gene'].str.contains(pat='-[1-9]$',
                                                regex=True), 'TRBV_gene'].str.replace('-', '-0', regex=False)
        df.loc[df['TRBJ_gene'].str.contains(pat='-[1-9]$', regex=True), 'TRBJ_gene'] = \
            df.loc[df['TRBJ_gene'].str.contains(pat='-[1-9]$',
                                                regex=True), 'TRBJ_gene'].str.replace('-', '-0', regex=False)
        if df.empty:
            error += 1
            print('\n' + f'{current}/{total}: file is empty')
            print(f'{error} files were empty after parsing and hence not saved' + '\n')
        else:
            df.to_csv(path_or_buf=(folder + f'MiXCR_parsed/{new_name}.tsv'), sep='\t', index=False)
            print(f'{current}/{total}: {new_name} is parsed and saved!' + '\n')
print(f'{error} files skipped')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REMAKE from here on with proper TCRex predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge those parsed MiXCR files with TCRex predictions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rename tcrex files so it will be easier to join
tcrex_path = f'/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/IR-Binder-000001_SARS-CoV-2_TCRex_predictions/'
for file in os.listdir(tcrex_path):
    new_file = file.split('_')[0] + '_tcrex.tsv'
    os.rename(os.path.join(tcrex_path, file), os.path.join(tcrex_path, new_file))

# list of all IDs of patients (filenames match the start of TCRex prediction file names)
id_list = [one_id.split('.')[0] for one_id in os.listdir(folder + 'MiXCR_parsed/')
           if one_id.startswith('Pt')]
id_list.sort()

# output location
if not os.path.exists(folder + 'MiXCR_TCRex_short'):
    os.makedirs(folder + 'MiXCR_TCRex_short')

# create empty dataframe where to combine all separate files to
tcrex_seq_all = pd.DataFrame()

current = 0
error = 0
for one_id in id_list:
    # sequencing file with counts
    seq = pd.read_csv(folder + f'MiXCR_parsed/{one_id}.tsv', sep='\t')
    # tcrex output (predictions) file (SarsCov2-only)
    try:
        tcrex = pd.read_csv(tcrex_path + f'{one_id}_tcrex.tsv', sep='\t', header=5)
        if tcrex.empty:
            print('\n' + f'No TCRex predctions for file {one_id}_tcrex.tsv' + '\n')
        # combine tcrex and sequencing files
        tcrex_seq = pd.merge(tcrex, seq, how="left", on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])


        # number of different nt sequences per CDR3 with the same V-J genes, epitope, etc.
        nt_seq = tcrex_seq.groupby(['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr'],
                                   as_index=False, sort=False)['nSeqCDR3'].count().rename(columns={'nSeqCDR3': 'N_nt_seq'})
        # combined counts for CDR3 with the same V-J genes, epitope, etc. but different nt sequences
        new_count = tcrex_seq.groupby(['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr'],
                                      as_index=False, sort=False)['count'].sum()
        # combined frequencies for CDR3 with the same V-J genes, epitope, etc. but different nt sequences
        new_freq = tcrex_seq.groupby(['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr'],
                                     as_index=False, sort=False)['freq'].sum()

        # step-by-step merging to update count and freq columns and substitude cdr3nt with N_nt_seq
        interim = pd.merge(nt_seq, new_count, how='inner',
                           on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
                               'score', 'bpr'])
        interim = pd.merge(interim, new_freq, how='inner',
                           on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology',
                               'score', 'bpr'])
        tcrex_seq = pd.merge(tcrex_seq[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr',
                                        'total_count', 'unique_TCRs', 'unique_CDR3s',
                                        'patient', 'type', 'study']], interim,
                             how='inner',
                             on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr'])

        # to have empty tcrex entry but not empty seq entry
        tcrex_seq[['total_count', 'unique_TCRs', 'unique_CDR3s',
                   'patient', 'type', 'study']] = seq[['total_count', 'unique_TCRs', 'unique_CDR3s',
                                                       'patient', 'type', 'study']]
        tcrex_seq['patient_id'] = one_id

        tcrex_seq.drop_duplicates(inplace=True)
        tcrex_seq.to_csv(path_or_buf=(folder + f'MiXCR_TCRex_short/{one_id}_mixcr_tcrex.tsv'), sep='\t', index=False)
        # add to one combined file for all patients
        tcrex_seq_all = tcrex_seq_all.append(tcrex_seq, ignore_index=True, sort=False)
        print(str(datetime.datetime.now()) + ':\n' +
              f'{current}/{len(id_list)}: Joined MiXCR and TCRex files of {one_id}' + '\n')
    except Exception as E:
        error += 1
        print('\n' + f'There was a problem with {error} tcrex files')
        print(f'{E}' + '\n' + f'for file {one_id}_tcrex.tsv' + '\n')
print(f'{error} tcrex files skipped')

# number of rows with nan after merging (no seq data for TCRex predictions)
# tcrex_seq_all.shape[0] - tcrex_seq_all.dropna().shape[0]  # 168 out of 789!!
zeros = tcrex_seq_all[tcrex_seq_all['N_nt_seq'] == 0].shape[0]
print(f'{zeros} out of {tcrex_seq_all.shape[0]} rows have no seq data for tcrex predictions')
# 167 out of 769 rows have no seq data for tcrex predictions; Pt-23-1 has nan

# PMe epitope data: common CV vs SARS-CoV-2 eps
ep_data = pd.read_csv('/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/epitope_species.txt', sep='\t')
tcrex_seq_all = pd.merge(tcrex_seq_all, ep_data, on='epitope', how='left')
tcrex_seq_all.fillna(value={'species': 0}, inplace=True)
tcrex_seq_all['species'] = tcrex_seq_all['species'].astype('int32')

tcrex_seq_all.to_csv(path_or_buf=(folder + 'MiXCR_shortTCRex_all.tsv'), sep='\t', index=False)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -> code below is in TCRexout_countSpec_IR_Binder for REMAKE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tcrex_seq_all = pd.read_csv(folder + 'MiXCR_shortTCRex_all.tsv', sep='\t')
# tcrex_seq_all = tcrex_seq_all[tcrex_seq_all['freq'] > 1 / 100000]


def tcr_counter_short(all_data, one_id):
    """

    :param all_data: df, a data frame with all TCRex predictions for all patients, cell, types, time points
    :param one_id: str
    :return: df
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
    tmp.drop(columns=['score', 'bpr'], inplace=True)

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

    # delete columns
    tmp.drop(columns=['count', 'freq', 'TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'species'],
             inplace=True)
    tmp.drop_duplicates(inplace=True)

    return tmp


# list of all IDs of patients (filenames match the start of TCRex prediction file names)
id_list = [one_id.split('.')[0] for one_id in os.listdir(folder + 'MiXCR_parsed/')
           if one_id.startswith('Pt')]
id_list.sort()

tcr_per_virus_short = pd.DataFrame()
# tmp df with TCR counts per covid/virus per TP per person per cell type
n = 0
for one_patient in id_list:
    n += 1
    counted = tcr_counter_short(tcrex_seq_all, one_patient)
    # append to the df with tcr count info for all individuals at all time points
    tcr_per_virus_short = tcr_per_virus_short.append(counted, ignore_index=True, sort=False)
    print(f"{n}/{len(id_list)}: " + one_patient + " TCR specificities have been counted!")

# add meta data from PMe file
meta_pme = pd.read_csv('/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/PMedata/SARShits09.csv', index_col=0)
tcr_per_virus_short = pd.merge(tcr_per_virus_short, meta_pme, on='patient')

# rename and add week column
tcr_per_virus_short.rename(columns={'days': 'day', 'disease': 'disease_severity'}, inplace=True)
tcr_per_virus_short.loc[tcr_per_virus_short['day'] <= 7, 'week'] = 1
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 7) & (tcr_per_virus_short['day'] <= 14), 'week'] = 2
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 14) & (tcr_per_virus_short['day'] <= 21), 'week'] = 3
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 21) & (tcr_per_virus_short['day'] <= 28), 'week'] = 4
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 28) & (tcr_per_virus_short['day'] <= 35), 'week'] = 5
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 35) & (tcr_per_virus_short['day'] <= 42), 'week'] = 6
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 42) & (tcr_per_virus_short['day'] <= 49), 'week'] = 7
tcr_per_virus_short.loc[(tcr_per_virus_short['day'] > 49) & (tcr_per_virus_short['day'] <= 56), 'week'] = 8
tcr_per_virus_short['week'] = tcr_per_virus_short['week'].map(lambda x: int(x))

tcr_per_virus_short.sort_values(by=['patient_id', 'day'], inplace=True)
tcr_per_virus_short = tcr_per_virus_short[['patient_id', 'total_count', 'unique_TCRs', 'unique_CDR3s',
                                           'disease_severity', 'day', 'week', 'N_SARS-CoV-2_Eps',
                                           'N_SARS-CoV-2_only_Eps',
                                           'N_CoV_Eps', 'Frac_CoV_in_SARS-CoV-2_Eps', 'Frac_SARS-CoV-2_only_Eps',
                                           'N_SARS-CoV-2_TCRs', 'Frac_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_only_TCRs',
                                           'Frac_SARS-CoV-2_only_TCRs', 'N_CoV_TCRs', 'Frac_CoV_TCRs',
                                           'Frac_CoV_in_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_crTCRs',
                                           'Frac_SARS-CoV-2_crTCRs', 'Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs',
                                           'Frac_CoV_in_crTCRs', 'SARShits', 'SARSsig', 'ORF1hits', 'ORF2hits',
                                           'uniquehits', 'otherhits', 'patient', 'type', 'study']]
tcr_per_virus_short.drop_duplicates(inplace=True)
tcr_per_virus_short.rename(columns={"patient": 'patient_sample'}, inplace=True)
tcr_per_virus_short['patient'] = tcr_per_virus_short['patient_sample'].str.split('-', expand=True)[0]
# save one file with all TCR counts
tcr_per_virus_short.to_csv(path_or_buf=(folder + 'TCRcounts_per_virus_short.tsv'), index=False, sep='\t')
