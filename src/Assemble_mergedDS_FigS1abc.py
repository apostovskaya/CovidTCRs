import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merging CD8+ split and mixed datasets,
# keeping only relevant for further analysis columns
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df1 = imseq_cd8[['study_patient_id', 'disease_severity', 'disease_status',
                 'critical_disease', 'day', 'week',
                 'total_count', 'unique_TCRs', 'Frac_CoV_TCRs',
                 'Frac_SARS-CoV-2_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
                 'N_CoV_TCRs', 'N_SARS-CoV-2_only_TCRs', 'N_SARS-CoV-2_TCRs',
                 'N_SARS-CoV-2_Eps', 'N_SARS-CoV-2_only_Eps', 'N_CoV_Eps',
                 'all_SARS-CoV-2_epitopes']].copy()
df1.rename(columns={'study_patient_id': 'patient_id'}, inplace=True)
df1['source'] = 'IMSEQ'
df2 = tcr_per_virus_pt[
    ['patient_id', 'disease_severity', 'disease_status', 'critical_disease',
     'day', 'week',
     'total_count', 'unique_TCRs',
     'Frac_CoV_TCRs', 'Frac_SARS-CoV-2_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
     'N_CoV_TCRs', 'N_SARS-CoV-2_only_TCRs', 'N_SARS-CoV-2_TCRs',
     'N_SARS-CoV-2_Eps', 'N_SARS-CoV-2_only_Eps',
     'N_CoV_Eps', 'all_SARS-CoV-2_epitopes']].copy()

# severe in this paper corresponds to critical in WHO grading
# (wasn't available yet when they were classifying patients)
df2.loc[df2['disease_severity'] == 'severe', 'disease_severity'] = 'critical'
df2['source'] = 'IR'

all_data = df1.append(df2)
# we are not interested in asymptomatic patients
all_data = all_data[all_data['disease_severity'] != 'asymptomatic']
all_data.sort_values(
    by=['critical_disease', 'disease_severity', 'day', 'patient_id'],
    inplace=True)
# clean up epitope column to remove [] and '' and transform epitopes per person
# to actual lists instead of list-looking strings
all_data['all_SARS-CoV-2_epitopes'] = all_data[
    'all_SARS-CoV-2_epitopes'].str.strip('[]').str.replace("'", "").str.split(
    ', ')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rescaling to [0:1] by subtracting min(value) and dividing by max(value)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# discussion about ways of rescaling
# https://stackoverflow.com/questions/703924/correct-way-to-standardize-scale-normalize-multiple-variables-following-power-la
df1['Frac_CoV_TCRs'] = (df1['Frac_CoV_TCRs'] - df1['Frac_CoV_TCRs'].min()) / \
                       df1['Frac_CoV_TCRs'].max()
df1['Frac_SARS-CoV-2_TCRs'] = (df1['Frac_SARS-CoV-2_TCRs'] - df1[
    'Frac_SARS-CoV-2_TCRs'].min()) / \
                              df1['Frac_SARS-CoV-2_TCRs'].max()
df1['Frac_SARS-CoV-2_only_TCRs'] = (df1['Frac_SARS-CoV-2_only_TCRs'] - df1[
    'Frac_SARS-CoV-2_only_TCRs'].min()) / \
                                   df1['Frac_SARS-CoV-2_only_TCRs'].max()

df2['Frac_CoV_TCRs'] = (df2['Frac_CoV_TCRs'] - df2['Frac_CoV_TCRs'].min()) / \
                       df2['Frac_CoV_TCRs'].max()
df2['Frac_SARS-CoV-2_TCRs'] = (df2['Frac_SARS-CoV-2_TCRs'] - df2[
    'Frac_SARS-CoV-2_TCRs'].min()) / \
                              df2['Frac_SARS-CoV-2_TCRs'].max()
df2['Frac_SARS-CoV-2_only_TCRs'] = (df2['Frac_SARS-CoV-2_only_TCRs'] - df2[
    'Frac_SARS-CoV-2_only_TCRs'].min()) / \
                                   df2['Frac_SARS-CoV-2_only_TCRs'].max()
all_data_scaled = df1.append(df2)

# we are not interested in asymptomatic patients (-1 row from 'mixed' DS)
all_data_scaled = all_data_scaled[
    all_data_scaled['disease_severity'] != 'asymptomatic']
all_data_scaled.sort_values(by=['critical_disease', 'patient_id', 'day'],
                            inplace=True)

all_data_scaled['%unique_TCRs'] = (all_data_scaled['unique_TCRs'] /
                                   all_data_scaled['total_count']) * 100

# calculating Epitope proportion out of uniqueTCRs
all_data_scaled['%CoV-common_Eps'] = 100 * all_data_scaled['N_CoV_Eps'] / \
                                     all_data_scaled['unique_TCRs']
all_data_scaled['%SC2-unique_Eps'] = 100 * all_data_scaled[
    'N_SARS-CoV-2_only_Eps'] / all_data_scaled['unique_TCRs']

# rescaling of N of epitopes to median repertoire size
median_TCRs = all_data_scaled['unique_TCRs'].median()  # 6394
all_data_scaled['scaled_N_CoV-common_Eps'] = all_data_scaled[
                                                 'N_CoV_Eps'] * median_TCRs / \
                                             all_data_scaled['unique_TCRs']
all_data_scaled['scaled_N_SC2-unique_Eps'] = all_data_scaled[
                                                 'N_SARS-CoV-2_only_Eps'] * median_TCRs \
                                             / all_data_scaled['unique_TCRs']

#
all_data_scaled['%CoV-common_TCRs'] = 100 * all_data_scaled['N_CoV_TCRs'] / \
                                      all_data_scaled['unique_TCRs']
all_data_scaled['%SC2-unique_TCRs'] = 100 * all_data_scaled[
    'N_SARS-CoV-2_only_TCRs'] / all_data_scaled['unique_TCRs']
all_data_scaled['%allSC2_TCRs'] = 100 * all_data_scaled['N_SARS-CoV-2_TCRs'] / \
                                  all_data_scaled['unique_TCRs']

all_data_scaled = all_data_scaled[['patient_id', 'disease_status',
                                   'disease_severity', 'critical_disease',
                                   'day', 'week',
                                   'total_count', 'unique_TCRs',
                                   '%unique_TCRs',
                                   'Frac_CoV_TCRs',
                                   'Frac_SARS-CoV-2_only_TCRs',
                                   'N_CoV_TCRs', '%CoV-common_TCRs',
                                   'N_SARS-CoV-2_only_TCRs',
                                   '%SC2-unique_TCRs',
                                   'N_SARS-CoV-2_TCRs', '%allSC2_TCRs',
                                   'N_CoV_Eps', 'scaled_N_CoV-common_Eps',
                                   '%CoV-common_Eps',
                                   'N_SARS-CoV-2_only_Eps',
                                   'scaled_N_SC2-unique_Eps',
                                   '%SC2-unique_Eps',
                                   'Frac_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_Eps',
                                   'all_SARS-CoV-2_epitopes', 'source']]
all_data_scaled.reset_index(drop=True, inplace=True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS1-c number of days sampled per patient
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
ax = sns.histplot(all_data_scaled.groupby('patient_id')['day'].nunique(),
                  stat="count", discrete=True, color='tab:green')
ax.set(xlabel='No. of data points per patient', ylabel='No. of patients')
plt.tight_layout()
plt.savefig(f'{folder_out}/FigS1c_N_dataPoints_perPerson.jpg', dpi=600)
plt.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS4extra - longitudinal dynamics per person (in case I need form the whole DS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
long_data_only = all_data_scaled[all_data_scaled.duplicated(subset='patient_id',
                                                            keep=False)].copy()
for parameter in ['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs']:
    mpl.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(6.6, 5))
    sns.lineplot(x='week', y=parameter, data=long_data_only,
                 hue_order=['Yes', 'No'], hue='critical_disease',
                 palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                 style='patient_id', markers=True, dashes=False)
    plt.ylabel(parameter)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=6)
    plt.tight_layout()
    plt.savefig(f'{folder_out}/FigS4_PersonDyn_{parameter}.png')
    plt.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TCR to Epitope ratio for Fig7abc
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_data_scaled['CoV-common_TCR-Ep_ratio'] = all_data_scaled['N_CoV_TCRs'] \
                                                      / all_data_scaled['N_CoV_Eps']
all_data_scaled['SC2-unique_TCR-Ep_ratio'] = all_data_scaled['N_SARS-CoV-2_only_TCRs'] \
                                                      / all_data_scaled['N_SARS-CoV-2_only_Eps']
# for statistics
crit = all_data_scaled[
    (all_data_scaled['critical_disease'] == 'Yes')]
noncrit = all_data_scaled[(all_data_scaled['critical_disease'] == 'No') &
                          (all_data_scaled['disease_status'] == 'recovered')]
noncrit_active = all_data_scaled[(all_data_scaled['critical_disease'] == 'No') &
                                 (all_data_scaled['disease_status'] == 'active')]

print(f"Median TCR/SC2-unique_Epitope ratio in critical patients: "
      f"{crit['SC2-unique_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{crit['SC2-unique_TCR-Ep_ratio'].describe()[3]}-"
      f"{crit['SC2-unique_TCR-Ep_ratio'].describe()[7]}]")
print(f"Median TCR/CoV-common_Epitope ratio in critical patients: "
      f"{crit['CoV-common_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{crit['CoV-common_TCR-Ep_ratio'].describe()[3]}-"
      f"{crit['CoV-common_TCR-Ep_ratio'].describe()[7]}]")

print(f"Median TCR/SC2-unique_Epitope ratio in non-critical active patients: "
      f"{noncrit_active['SC2-unique_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit_active['SC2-unique_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit_active['SC2-unique_TCR-Ep_ratio'].describe()[7]}]")
print(f"Median TCR/CoV-common_Epitope ratio in non-critical active patients: "
      f"{noncrit_active['CoV-common_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit_active['CoV-common_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit_active['CoV-common_TCR-Ep_ratio'].describe()[7]}]")

print(f"Median TCR/SC2-unique_Epitope ratio in non-critical recovered patients: "
      f"{noncrit['SC2-unique_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit['SC2-unique_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit['SC2-unique_TCR-Ep_ratio'].describe()[7]}]")
print(f"Median TCR/CoV-common_Epitope ratio in non-critical recovered patients: "
      f"{noncrit['CoV-common_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit['CoV-common_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit['CoV-common_TCR-Ep_ratio'].describe()[7]}]")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# keeping only max of a week for each person for further week Max analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_data_scaled_week_max = all_data_scaled.copy()
# assign to all days of the same week of the same patient
# the maximum values of the parameter for this week
for parameter in ['total_count', 'unique_TCRs', '%unique_TCRs',
                  'Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
                  'N_CoV_TCRs', '%CoV-common_TCRs', 'N_SARS-CoV-2_only_TCRs',
                  '%SC2-unique_TCRs',
                  'N_CoV_Eps', 'scaled_N_CoV-common_Eps', '%CoV-common_Eps',
                  'N_SARS-CoV-2_only_Eps', 'scaled_N_SC2-unique_Eps',
                  '%SC2-unique_Eps',
                  'Frac_SARS-CoV-2_TCRs', 'N_SARS-CoV-2_Eps',
                  'N_SARS-CoV-2_TCRs', '%allSC2_TCRs']:
    if parameter == 'N_SARS-CoV-2_Eps':
        # to keep the maximum of recognized epitopes from the same week
        all_data_scaled_week_max.sort_values(['patient_id', 'week', parameter],
                                             ascending=False, inplace=True)
        all_data_scaled_week_max['all_SARS-CoV-2_epitopes'] = \
            all_data_scaled_week_max.loc[
                all_data_scaled_week_max.groupby(['patient_id', 'week'])[
                    parameter].idxmax(), 'all_SARS-CoV-2_epitopes']
        all_data_scaled_week_max['all_SARS-CoV-2_epitopes'].fillna(
            method='ffill', inplace=True)
        all_data_scaled_week_max[parameter] = \
            all_data_scaled_week_max.groupby(['patient_id', 'week'])[
                parameter].transform('max')
    else:
        all_data_scaled_week_max[parameter] = \
            all_data_scaled_week_max.groupby(['patient_id', 'week'])[
                parameter].transform('max')

# sort so active is before recovered because
# if there are two data points in one week when one point is from active
# and another from recovered state, we want to keep the active record
all_data_scaled_week_max.sort_values(
    by=['critical_disease', 'patient_id', 'week', 'disease_status'],
    ascending=True, inplace=True, ignore_index=True)
# keep only 1 record with the highest value of the week
all_data_scaled_week_max.drop_duplicates(['patient_id', 'week'],
                                         keep='first', inplace=True)
all_data_scaled_week_max.to_csv(
    f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t', index=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS1-a hist of merged data distribution per disease severity as is
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
sns.countplot(x='week', hue='disease_severity', data=all_data_scaled_week_max,
              palette={'mild': 'tab:cyan', 'moderate': 'tab:blue',
                       'severe': 'tab:grey',
                       'critical': 'tab:orange', 'fatal': 'tab:red'},
              hue_order=['mild', 'moderate', 'severe', 'critical', 'fatal'],
              dodge=True)  # palette='gist_heat_r')
plt.legend(loc='upper right', title='Disease severity')
plt.xlabel('week')
plt.ylabel('No. of patients')
plt.savefig(f'{folder_out}/FigS1a_severity_distribution.jpg', dpi=600)
plt.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS1-b hist of merged data distribution per patient group - critical/non-critical
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
sns.countplot(x='week', hue='critical_disease', data=all_data_scaled_week_max,
              hue_order=['No', 'Yes'],
              dodge=True, palette={'No': 'tab:blue', 'Yes': 'tab:orange'})
plt.legend(loc='upper right', title='Patient group',
           labels=['Non-critical', 'Critical'])
plt.xlabel('week')
plt.ylabel('No. of patients')
plt.savefig(f'{folder_out}/FigS1b_criticals_distribution.jpg', dpi=600)
plt.close()
