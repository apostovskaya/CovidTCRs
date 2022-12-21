import sys
import os
import pandas as pd

assert len(sys.argv) == 5, "The input is incorrect. Expected: (1) script name; (2) main folder; " \
                           "(3) input file in .tsv; (4) input metadata file in .tsv; (5) output file prefix"

print(f'\n"{sys.argv[0].split("/")[-1]}" script is being executed')

# input/output folder path
folder = os.path.abspath(sys.argv[1])  # 'IMSEQ'
# folder = '/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/'
file_in = sys.argv[2]  # 'TCRex_processed/per_patient/new_intersected_double_TCRcounts_per_virus.tsv'
file_in_dir = '/'.join(file_in.split('/')[0:-1])
filename_in = file_in.split('/')[-1].split('.')[0]
meta_file = sys.argv[3]  # 'Metadata/clean_IMSEQ_age.tsv'
prefix_out = sys.argv[4]  # 'MetaForPlotting'

orig_stdout = sys.stdout
log_f = open(f'{folder}/{file_in_dir}/LogFile_addMeta_{filename_in}.txt', 'w+')
sys.stdout = log_f

print("~ Executed script: " + sys.argv[0])
print("~ Main folder: " + folder)
print("~ Input file: " + file_in)
print("~ Provided metafile: " + meta_file)
print("~ Provided prefix for the output file name: " + prefix_out)
print("~ The output file will be saved in the same directory as the input file: " + file_in_dir)

tcr_per_virus_inter = pd.read_csv(f'{folder}/{file_in}', sep='\t')

# only 1 ICU patient, and profile is very similar to critical, so will unite those two groups into critical
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'INZO', 'disease_severity'] = 'critical'
tcr_per_virus_inter.sort_values(by=['disease_severity', 'day', 'new_id', 'T-cell_type'], inplace=True)
print('\n1 ICU patient is renamed to be in critical category.')

# want to have the same sorting, but not alphabetical disease_severity, but asympt-moder-crit
# Create the dictionary that defines the order for sorting
sorter = ['asymptomatic', 'moderate', 'critical']
sorterIndex = dict(zip(sorter, range(len(sorter))))  # {'asymptomatic': 0, 'moderate': 1, 'critical': 2}
# Generate a rank column that will be used to sort the dataframe numerically
tcr_per_virus_inter['group_Rank'] = tcr_per_virus_inter['disease_severity'].map(sorterIndex)
tcr_per_virus_inter.sort_values(by=['group_Rank', 'day', 'new_id', 'T-cell_type'], inplace=True)
tcr_per_virus_inter.drop('group_Rank', 1, inplace=True)

# hits = N specific unique TCRs / N unique TCRs (to compare with PMe pre-print)
tcr_per_virus_inter['Hits_SARS-CoV-2'] = tcr_per_virus_inter['N_SARS-CoV-2_only_TCRs'] / tcr_per_virus_inter[
    'unique_TCRs']
tcr_per_virus_inter['Hits_CoV'] = tcr_per_virus_inter['N_CoV_TCRs'] / tcr_per_virus_inter['unique_TCRs']
tcr_per_virus_inter['Hits'] = (tcr_per_virus_inter['N_SARS-CoV-2_only_TCRs'] +
                               tcr_per_virus_inter['N_CoV_TCRs']) / tcr_per_virus_inter['unique_TCRs']
print('\nHits column is created, where hits = N specific unique TCRs / N unique TCRs')

# to add age
age = pd.read_csv(f'{folder}/{meta_file}', sep='\t')
age = age[['study_patient_id', 'age', 'comorbidities', 'disease_severity']]
tcr_per_virus_inter = pd.merge(tcr_per_virus_inter, age, on=['study_patient_id', 'disease_severity'], how='left')
# no age data for patients 'imseq27', 'imseq25', 'imseq28'
no_age_info = tcr_per_virus_inter.loc[tcr_per_virus_inter['age'].isna(), 'study_patient_id'].unique()
print(f'\nNo age information is available for the following patients: \n{no_age_info}')
# no comorbidities data for 'covidam12', 'covidam10', 'imseq27', 'imseq25', 'imseq28'
no_comorb_info = tcr_per_virus_inter.loc[tcr_per_virus_inter['comorbidities'].isna(), 'study_patient_id'].unique()
print(f'No information about comorbidities is available for the following patients: \n{no_comorb_info}')
# tcr_per_virus_inter.loc[tcr_per_virus_inter['age'].notna(), 'age'] = \
#     tcr_per_virus_inter.loc[tcr_per_virus_inter['age'].notna(), 'age'].astype('int')

# UPD on patient imseq28 - severe, not critical (mistake in the original excel file)
tcr_per_virus_inter.loc[tcr_per_virus_inter['study_patient_id'] == 'imseq28', 'disease_severity'] = 'severe'

# add disease severity category label to compare with IR-Binder data
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'asymptomatic', 'disease_burden'] = 'low'
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'moderate', 'disease_burden'] = 'medium'
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'severe', 'disease_burden'] = 'medium'
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'critical', 'disease_burden'] = 'high'

# add critical disease label to compare with IR-Binder data
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'asymptomatic', 'critical_disease'] = 'No'
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'moderate', 'critical_disease'] = 'No'
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'severe', 'critical_disease'] = 'No'
tcr_per_virus_inter.loc[tcr_per_virus_inter['disease_severity'] == 'critical', 'critical_disease'] = 'Yes'

tcr_per_virus_inter = tcr_per_virus_inter[['new_id', 'patient_id', 'T-cell_type', 'study_patient_id',
                                           'disease_severity', 'disease_burden', 'critical_disease', 'day', 'week',
                                           'total_count', 'unique_TCRs', 'unique_CDR3s',
                                           'N_SARS-CoV-2_Eps', 'N_SARS-CoV-2_only_Eps',
                                           'N_CoV_Eps', 'Frac_CoV_in_SARS-CoV-2_Eps', 'Frac_SARS-CoV-2_only_Eps',
                                           'N_viral_Eps', 'N_SARS-CoV-2_TCRs', 'Frac_SARS-CoV-2_TCRs',
                                           'N_SARS-CoV-2_only_TCRs', 'Frac_SARS-CoV-2_only_TCRs', 'N_CoV_TCRs',
                                           'Frac_CoV_TCRs', 'Frac_CoV_in_SARS-CoV-2_TCRs', 'N_viral_TCRs',
                                           'Frac_viral_TCRs', 'N_SARS-CoV-2_crTCRs', 'Frac_SARS-CoV-2_crTCRs',
                                           'Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs', 'Frac_CoV_in_crTCRs',
                                           'N_SARS-CoV-2_viral_crTCRs', 'Frac_SARS-CoV-2_viral_crTCRs',
                                           'Frac_CoV_in_viral_crTCRs', 'cr_viruses', 'all_epitopes',
                                           'all_SARS-CoV-2_epitopes', 'Hits_SARS-CoV-2', 'Hits_CoV', 'Hits', 'age',
                                           'comorbidities']]

df = tcr_per_virus_inter.drop_duplicates(subset='study_patient_id', ignore_index=True)

# age per group
ages = {}
ages['critical'] = df.loc[df['disease_severity'] == 'critical', 'age']
ages['moderate'] = df.loc[df['disease_severity'] == 'moderate', 'age']
ages['asymptomatic'] = df.loc[df['disease_severity'] == 'asymptomatic', 'age']
ages['not_critical'] = df.loc[df['disease_severity'] != 'critical', 'age']
print(ages)
# {'critical': 10     NaN
#  11     NaN
#  12    60.0  # imseq21 inrease of SC2-hits & Fracs as in moderate
#  13    50.0
#  Name: age, dtype: float64,
#  'moderate': 3     NaN
#  4    53.0
#  5    57.0
#  6    36.0  # imseq12 from  df.loc[df.age == 36, 'study_patient_id']; obesity
#  7    61.0
#  8    56.0
#  9    65.0
#  Name: age, dtype: float64,
#  'asymptomatic': 0    26.0
#  1    24.0
#  2    33.0  # imseq23

# median age per group
ages_median = {}
ages_median['high'] = df.loc[df['disease_burden'] == 'high', 'age'].median()
ages_median['medium'] = df.loc[df['disease_burden'] == 'medium', 'age'].median()
ages_median['low'] = df.loc[df['disease_burden'] == 'low', 'age'].median()
ages_median['critical'] = df.loc[df['critical_disease'] == 'Yes', 'age'].median()
ages_median['not_critical'] = df.loc[df['critical_disease'] == 'No', 'age'].median()

# {'critical': 55.0, 'moderate': 56.5, 'asymptomatic': 26.0}
print(f'\nMedian age per group: \n{ages_median}')

sys.stdout = orig_stdout
log_f.close()

print(f'\nMedian age per group: \n{ages_median}')
print(f'\nMetadata has been added to the {filename_in}')

# exclude asymptomatic as we are not going to analyse this group
tcr_per_virus_inter = tcr_per_virus_inter[tcr_per_virus_inter['disease_severity'] != 'asymptomatic']
tcr_per_virus_inter.to_csv(path_or_buf=f'{folder}/{file_in_dir}/{prefix_out}_{filename_in}.tsv', sep='\t', index=False)

# python src/processPredictions/addMetadata_splitDS.py 'data/splitDS_IMSEQ'
# 'TCRex_processed/viralTCRexPreds_withCounts_splitDS_intersected2x_meta_perSpecificity.tsv' 'clean_IMSEQ_age.tsv' 'MetaForPlotting'