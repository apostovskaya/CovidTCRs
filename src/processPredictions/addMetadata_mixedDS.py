import sys
import os
import pandas as pd

assert len(sys.argv) == 5, "The input is incorrect. Expected: (1) script name; (2) main folder; " \
                           "(3) input file in .tsv; (4) input metadata file in .tsv; (5) output file prefix"

print(f'\n"{sys.argv[0].split("/")[-1]}" script is being executed\n')

# input/output folder path
folder = os.path.abspath(sys.argv[1])  # 'IMSEQ'
# folder = '/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/'
file_in = sys.argv[2]  # 'IR-Binder_PRJEB38339/full_TCRcounts_per_virus_Pt.tsv'
file_in_dir = '/'.join(file_in.split('/')[0:-1])
filename_in = file_in.split('/')[-1].split('.')[0]
meta_file = sys.argv[3]  # 'Metadata/IR_Binder_Meta_manual.csv'
prefix_out = sys.argv[4]  # 'MetaForPlotting'

print("~ Executed script: " + sys.argv[0])
print("~ Main folder: " + folder)
print("~ Input file: " + file_in)
print("~ Provided metafile: " + meta_file)
print("~ Provided prefix for the output file name: " + prefix_out)
print("~ The output file will be saved in the same directory as the input file: " + file_in_dir)

tcr_per_virus_pt = pd.read_csv(f'{folder}/{file_in}', sep='\t')
meta = pd.read_csv(f'{folder}/{meta_file}', sep=';')

# rename columns
tcr_per_virus_pt.rename(columns={'patient_id': 'sample_id', 'patient': 'patient_id',
                                 'disease_severity': 'disease_status'}, inplace=True)
# add disease severity from the meta file
tcr_per_virus_pt = tcr_per_virus_pt.merge(meta[['patient_id', 'disease_severity']], how='left', on='patient_id')

# add disease severity category label to compare with IMSEQ data
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'asymptomatic', 'disease_burden'] = 'low'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'mild', 'disease_burden'] = 'low'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'moderate', 'disease_burden'] = 'medium'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'severe', 'disease_burden'] = 'high'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'fatal', 'disease_burden'] = 'high'

# add critical disease label to compare with IMSEQ data
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'asymptomatic', 'critical_disease'] = 'No'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'mild', 'critical_disease'] = 'No'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'moderate', 'critical_disease'] = 'No'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'severe', 'critical_disease'] = 'Yes'
tcr_per_virus_pt.loc[tcr_per_virus_pt['disease_severity'] == 'fatal', 'critical_disease'] = 'Yes'

# median age per group
# -> add values from pdf first
# df = tcr_per_virus_pt.drop_duplicates(subset='patient_id', ignore_index=True)
# ages_median = {}
# ages_median['high'] = df.loc[df['disease_burden'] == 'high', 'age_range'].median()
# ages_median['medium'] = df.loc[df['disease_burden'] == 'medium', 'age_range'].median()
# ages_median['low'] = df.loc[df['disease_burden'] == 'low', 'age_range'].median()
# print(f'\nMedian age per group: \n{ages_median}')

tcr_per_virus_pt = tcr_per_virus_pt[['patient_id', 'sample_id', 'total_count', 'unique_TCRs', 'unique_CDR3s',
                                     'disease_status', 'disease_severity', 'disease_burden', 'critical_disease',
                                     'day', 'week', 'N_SARS-CoV-2_Eps', 'N_SARS-CoV-2_only_Eps', 'N_CoV_Eps',
                                     'Frac_CoV_in_SARS-CoV-2_Eps', 'Frac_SARS-CoV-2_only_Eps', 'N_SARS-CoV-2_TCRs',
                                     'Frac_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_only_TCRs',
                                     'Frac_SARS-CoV-2_only_TCRs', 'N_CoV_TCRs',
                                     'Frac_CoV_TCRs', 'Frac_CoV_in_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_crTCRs',
                                     'Frac_SARS-CoV-2_crTCRs', 'Frac_SARS-CoV-2_crTCRs_in_SARS-CoV-2_TCRs',
                                     'Frac_CoV_in_crTCRs', 'Hits_SARS-CoV-2', 'Hits_CoV', 'Hits', 'SARShits',
                                     'SARSsig', 'ORF1hits', 'ORF2hits', 'uniquehits', 'otherhits',
                                     'N_viral_TCRs', 'Frac_viral_TCRs', 'N_SARS-CoV-2_viral_crTCRs',
                                     'Frac_SARS-CoV-2_viral_crTCRs', 'Frac_CoV_in_viral_crTCRs',
                                     'cr_viruses', 'all_epitopes', 'all_SARS-CoV-2_epitopes',
                                     'patient_sample', 'type', 'study']]

tcr_per_virus_pt.to_csv(path_or_buf=f'{folder}/{file_in_dir}/{prefix_out}_{filename_in}.tsv', sep='\t', index=False)

print(f'\nMetadata has been added to the {filename_in}')

# python src/processPredictions/addMetadata_mixedDS.py 'data/mixedDS_IR_Binder' 'TCRex_processed/viralTCRexPreds_withCounts_mixedDS_perSpecificity_Pt.tsv' 'IR_Binder_Meta_manual.csv' 'MetaForPlotting'