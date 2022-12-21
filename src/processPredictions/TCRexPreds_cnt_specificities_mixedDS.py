import os
import pandas as pd
import numpy as np
from Functions import tcr_counter_mixed_ds

folder = './data/mixedDS_IR_Binder'

# one file with all TCRex predictions for viruses
tcrex_seq_all = pd.read_csv(f'{folder}/TCRex_processed/viralTCRexPreds_withCounts_mixedDS.tsv', sep='\t')  # 5458, 18

# filter out CDR3s with frequency less than 1 in 100.000 to compensate for different seq depth
tcrex_seq_all = tcrex_seq_all[tcrex_seq_all['freq'] > 1 / 100000]  # 5458 -> 4525
# filter out CDR3s with a score < 0.9
tcrex_seq_all = tcrex_seq_all[tcrex_seq_all['score'] >= 0.9]  # 4525 -> 840

# create files for Patients and HD separately, just in case
tcrex_viral_Pt = tcrex_seq_all[tcrex_seq_all['type'] == 'Pt'].copy()
tcrex_viral_HD = tcrex_seq_all[tcrex_seq_all['type'] == 'HD'].copy()
# RUN ONLY ONCE, DON'T REPEAT: save files for Patients and HD separately, just in case
tcrex_viral_Pt.to_csv(path_or_buf=(folder + '/TCRex_processed/viralTCRexPreds_withCounts_mixedDS_filtered_Pt.tsv'),
                      index=False, sep='\t')
tcrex_viral_HD.to_csv(path_or_buf=(folder + '/TCRex_processed/viralTCRexPreds_withCounts_mixedDS_filtered_HD.tsv'),
                      index=False, sep='\t')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PATIENTS ONLY
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# list of all IDs of patients (filenames match the start of TCRex prediction file names)
# id_list = [one_id.split('.')[0] for one_id in os.listdir(folder + 'MiXCR_parsed/')
#            if one_id.startswith('Pt')]
id_list_pt = tcrex_viral_Pt['patient_id'].unique().tolist()
id_list_pt.sort()

tcr_per_virus = pd.DataFrame()
# tmp df with TCR counts per covid/virus per TP per person
n = 0
for one_patient in id_list_pt:
    n += 1
    counted = tcr_counter_mixed_ds(tcrex_seq_all, one_patient)
    # append to the df with tcr count info for all individuals at all time points
    tcr_per_virus = tcr_per_virus.append(counted, ignore_index=True, sort=False)
    print(f"{n}/{len(id_list_pt)}: " + one_patient + " TCR specificities have been counted!")

# hits = # specific unique TCR / # unique TCRs (to compare with PMe pre-print)
tcr_per_virus['Hits_SARS-CoV-2'] = tcr_per_virus['N_SARS-CoV-2_only_TCRs'] / tcr_per_virus['unique_TCRs']
tcr_per_virus['Hits_CoV'] = tcr_per_virus['N_CoV_TCRs'] / tcr_per_virus['unique_TCRs']
tcr_per_virus['Hits'] = (tcr_per_virus['N_SARS-CoV-2_only_TCRs'] +
                         tcr_per_virus['N_CoV_TCRs']) / tcr_per_virus['unique_TCRs']

# add meta data from PMe file
meta_pme = pd.read_csv('/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/PMedata/SARShits09.csv', index_col=0)
tcr_per_virus = pd.merge(tcr_per_virus, meta_pme, on='patient')

# rename and add week column
tcr_per_virus.rename(columns={'days': 'day', 'disease': 'disease_severity'}, inplace=True)
tcr_per_virus.loc[tcr_per_virus['day'] <= 7, 'week'] = 1
tcr_per_virus.loc[(tcr_per_virus['day'] > 7) & (tcr_per_virus['day'] <= 14), 'week'] = 2
tcr_per_virus.loc[(tcr_per_virus['day'] > 14) & (tcr_per_virus['day'] <= 21), 'week'] = 3
tcr_per_virus.loc[(tcr_per_virus['day'] > 21) & (tcr_per_virus['day'] <= 28), 'week'] = 4
tcr_per_virus.loc[(tcr_per_virus['day'] > 28) & (tcr_per_virus['day'] <= 35), 'week'] = 5
tcr_per_virus.loc[(tcr_per_virus['day'] > 35) & (tcr_per_virus['day'] <= 42), 'week'] = 6
tcr_per_virus.loc[(tcr_per_virus['day'] > 42) & (tcr_per_virus['day'] <= 49), 'week'] = 7
tcr_per_virus.loc[(tcr_per_virus['day'] > 49) & (tcr_per_virus['day'] <= 56), 'week'] = 8
tcr_per_virus['week'] = tcr_per_virus['week'].map(lambda x: int(x))

tcr_per_virus.sort_values(by=['patient_id', 'day'], inplace=True)
# reorder all exiting columns
tcr_per_virus = tcr_per_virus[['patient_id', 'total_count', 'unique_TCRs', 'unique_CDR3s',
                               'disease_severity', 'day', 'week', 'N_SARS-CoV-2_Eps',
                               'N_SARS-CoV-2_only_Eps',
                               'N_CoV_Eps', 'Frac_CoV_in_SARS-CoV-2_Eps', 'Frac_SARS-CoV-2_only_Eps',
                               'N_SARS-CoV-2_TCRs', 'Frac_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_only_TCRs',
                               'Frac_SARS-CoV-2_only_TCRs', 'N_CoV_TCRs', 'Frac_CoV_TCRs',
                               'Frac_CoV_in_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_crTCRs',
                               'Frac_SARS-CoV-2_crTCRs', 'Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs',
                               'Frac_CoV_in_crTCRs', 'Hits_SARS-CoV-2', 'Hits_CoV', 'Hits',
                               'SARShits', 'SARSsig', 'ORF1hits', 'ORF2hits', 'uniquehits', 'otherhits',
                               'N_viral_TCRs', 'Frac_viral_TCRs', 'N_SARS-CoV-2_viral_crTCRs',
                               'Frac_SARS-CoV-2_viral_crTCRs', 'Frac_CoV_in_viral_crTCRs', 'cr_viruses',
                               'all_epitopes', 'all_SARS-CoV-2_epitopes', 'patient', 'type', 'study']]
# tcr_per_virus.drop_duplicates(inplace=True)  # 85, 35 -> 63, 35; now 63, 43
tcr_per_virus.rename(columns={"patient": 'patient_sample'}, inplace=True)
# different days of the same patient are in patient_sample (1-1, 1-2,...); patient - 1, 2, ...
tcr_per_virus['patient'] = tcr_per_virus['patient_sample'].str.split('-', expand=True)[0]
# save one file with all TCR counts
tcr_per_virus.to_csv(path_or_buf=(folder + '/TCRex_processed/viralTCRexPreds_withCounts_mixedDS_perSpecificity_Pt.tsv'), index=False, sep='\t')