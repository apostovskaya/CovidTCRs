import pandas as pd

# import
folder = './data'
folder_out = './results/Supplementary'
# PROCESSED DATA
file_in1 = 'mixedDS_IR_Binder/TCRex_processed/MetaForPlotting_viralTCRexPreds_withCounts_mixedDS_perSpecificity_Pt.tsv'
file_in2 = 'splitDS_IMSEQ/TCRex_processed/MetaForPlotting_viralTCRexPreds_withCounts_splitDS_intersected2x_meta_perSpecificity.tsv'

# IR-Binder public mixed dataset, patient data
tcr_per_virus_pt = pd.read_csv(f'{folder}/{file_in1}', sep='\t')
# IMSEQ in-house split dataset, only TCRs that occurred in at least 2 repeats
tcr_per_virus_inter = pd.read_csv(f'{folder}/{file_in2}', sep='\t')
# all patients are active in this dataset, adding a column to match mixedDS
tcr_per_virus_inter['disease_status'] = 'active'
# separating CD8 TCRs which will undergo further analysis
imseq_cd8 = tcr_per_virus_inter[tcr_per_virus_inter['T-cell_type'] == 'cd8']

# merging all data (for the supplementary table)
df1 = tcr_per_virus_inter.copy()
df1.drop(columns='patient_id', inplace=True)
df1.rename(columns={'study_patient_id': 'patient_id'}, inplace=True)
df1['source'] = 'IMSEQ'
df2 = tcr_per_virus_pt.copy()
# severe in this paper corresponds to critical in WHO grading
# (wasn't available yet when they were classifying patients)
df2.loc[df2['disease_severity'] == 'severe', 'disease_severity'] = 'critical'
df2['source'] = 'IR'
df2['T-cell_type'] = 'cd8'
df2 = df2[df2['disease_severity'] != 'asymptomatic']
common_cols = list(set(df1.columns.tolist()) & set(df2.columns.tolist()))
df1 = df1[common_cols]
df2 = df2[common_cols]
all_data = df1.append(df2)
all_data.sort_values(by=['critical_disease', 'disease_severity', 'day', 'patient_id'],
                     ignore_index=True, inplace=True)
# clean up epitope column to remove [] and '' and transform epitopes per person
# to actual lists instead of list-looking strings
all_data['all_SARS-CoV-2_epitopes'] = all_data['all_SARS-CoV-2_epitopes'].str.strip('[]').str.replace("'", "").str.split(', ')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Supplementary table S2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppT2 = all_data[['patient_id', 'disease_status',
                   'critical_disease', 'disease_severity', 'day', 'week',
                   'T-cell_type', 'total_count', 'unique_TCRs',
                   'all_SARS-CoV-2_epitopes', 'source']].copy()
suppT2.sort_values(by=['critical_disease', 'disease_severity', 'patient_id', 'day', 'T-cell_type'], inplace=True)
suppT2.reset_index(inplace=True, drop=True)
suppT2.loc[suppT2['critical_disease'] == 'Yes', 'critical_disease'] = 'critical'
suppT2.loc[suppT2['critical_disease'] == 'No', 'critical_disease'] = 'non-critical'
suppT2.loc[suppT2['source'] == 'IR', 'dataset'] = '"mixed"'
suppT2.loc[suppT2['source'] == 'IMSEQ', 'dataset'] = '"split"'
suppT2.loc[suppT2['source'] == 'IR', 'source'] = 'Schultheiss et al.'
suppT2.loc[suppT2['source'] == 'IMSEQ', 'source'] = 'in-house'
suppT2.rename(columns={'patient_id': 'patient ID', 'disease_severity': 'COVID-19 severity',
                       'critical_disease': 'Patient group',
                       'T-cell_type': 'T-cell population',
                       'total_count': 'total No. of TCRs',
                       'unique_TCRs': 'No. of unique TCRs',
                       'disease_status': 'COVID19 status',
                       'all_SARS-CoV-2_epitopes': 'recognized SC2 epitopes'},
              inplace=True)
suppT2 = suppT2[['patient ID', 'COVID19 status', 'COVID-19 severity',
                 'Patient group', 'day', 'week', 'T-cell population',
                 'total No. of TCRs', 'No. of unique TCRs',
                 'recognized SC2 epitopes', 'dataset', 'source']]
suppT2.to_csv(path_or_buf=f'{folder_out}/TableS2_PatientInfo.tsv', sep='\t', index=False)
