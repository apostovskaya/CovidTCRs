import os
import pandas as pd
import datetime
from Functions import combine_tcrex

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                   File processing: load in, combine together, save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IR-Binder
folder = './data/mixedDS_IR_Binder'
# specify path to where TCRex output files are stored
directory = f'{folder}/TCRex_results'
# where to save processed TCRex output files
if not os.path.exists(f'{folder}/TCRex_processed'):
    os.makedirs(f'{folder}/TCRex_processed')
    print(f'Processed TCRex files will be saved in {folder}/TCRex_processed')
directory_new = f'{folder}/TCRex_processed'

# from ViralTCRS in PhD_2020
# list of all IDs of patients (folder names in the TCRex_results folder)
id_list = [one_id for one_id in os.listdir(directory)
           if os.path.isdir(os.path.join(directory, one_id))]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add seq counts and frequencies to each of individual tcrex output files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for sample in id_list:
    # tcrex output (predictions) file
    print(str(datetime.datetime.now()) + ':\n' + f'Started combining TCRex predictions with sequencing counts data '
                                                 f'for {sample}')
    file_tcrex = f'{directory}/{sample}/filtered_results.csv'
    tcrex = pd.read_csv(file_tcrex, sep='\t', header=6)
    tcrex['patient_id'] = sample

    # sequencing files with counts
    file_seq = f'{folder}/MiXCR_parsed/{sample}.tsv'
    seq = pd.read_csv(file_seq, sep='\t')

    # combine tcrex and sequencing files
    tcrex_seq = pd.merge(tcrex, seq, how="outer", on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene'])
    tcrex_seq.drop_duplicates(inplace=True)

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
                                    'patient_id', 'total_count', 'unique_TCRs', 'unique_CDR3s',
                                    'patient', 'type', 'study']], interim, how='inner',
                         on=['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'epitope', 'pathology', 'score', 'bpr'])

    tcrex_seq.drop_duplicates(inplace=True)
    tcrex_seq.sort_values(by=['bpr', 'pathology', 'epitope', 'CDR3_beta'], inplace=True)
    tcrex_seq.to_csv(f'{directory}/{sample}/results_with_seq_counts.csv', index=False)
    print(f"{sample} files have been processed!")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine individual sample files into one; keep only viral TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tcrex_data = combine_tcrex(folder, id_list, 'results_with_seq_counts.csv')
# All predicted epitopes: ['SARS-CoV-2', 'HIV', 'EBV', 'CMV', 'HCV', 'DENV1', 'Melanoma', 'Influenza',
# 'DENV2', 'DENV3/4', 'HTLV1', 'HSV2', 'Multiple Myeloma', 'YellowFeverVirus']

tcrex_viral = tcrex_data.loc[~((tcrex_data.pathology == 'Melanoma') | (tcrex_data.pathology == 'Multiple Myeloma'))]
tcrex_viral = tcrex_viral.sort_values(by=['patient_id', 'bpr', 'pathology', 'epitope', 'CDR3_beta'])

# PMe epitope data: common CV vs SARS-CoV-2 eps
ep_data = pd.read_csv('./data/epitope_species.txt', sep='\t')
tcrex_viral = pd.merge(tcrex_viral, ep_data, on='epitope', how='left')
tcrex_viral.fillna(value={'species': 0}, inplace=True)
tcrex_viral['species'] = tcrex_viral['species'].astype('int32')

# save one file with all TCRex predictions for viruses
tcrex_viral.to_csv(f'{directory_new}/viralTCRexPreds_withCounts_mixedDS.tsv', sep='\t', index=False)