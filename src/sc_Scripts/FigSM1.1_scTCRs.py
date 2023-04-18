import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

folder = './data/scTCRs/E_MTAB_9357'
folder_out = './results'

tcr_per_virus = pd.read_csv(
    f'{folder}/TCRexPreds_withCounts_sc9357DS_clustered.tsv', sep='\t')
meta = pd.read_excel(f'{folder}/Human subject details.xlsx',
                     sheet_name='S1.1 Patient Clinical Data',
                     engine='openpyxl')
meta2 = pd.read_excel(f'{folder}/Human subject details.xlsx',
                      sheet_name='Ordinal Scale for Grading COVID',
                      engine='openpyxl')
meta_tp = pd.read_csv(f'{folder}/patient_timePoints.tsv', sep='\t')
meta_long = meta_tp.melt(id_vars=['patient_id', 'Onset'],
                         value_vars=['TP1_weeks', 'TP2_weeks'],
                         var_name='TP', value_name='week')
meta_long.sort_values(by=['patient_id', 'week'], inplace=True)
meta_long['TP'] = meta_long['TP'].str.replace('TP', '')
meta_long['TP'] = meta_long['TP'].str.replace('_weeks', '')
meta_long['TP'] = meta_long['TP'].astype(int)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add disease severity
meta['patient_id'] = meta['Sample ID'].str.split('-', expand=True)[0]
meta['TP'] = meta['Sample ID'].str.split('-', expand=True)[1]
meta['Who Ordinal Scale'] = meta[['Who Ordinal Scale']].replace({'1 or ': ''},
                                                                regex=True)
meta = meta[meta['Who Ordinal Scale'].notna()]
meta = meta.astype({'Who Ordinal Scale': int, 'TP': int})
meta.loc[(meta['Who Ordinal Scale'] >= 6) | ((meta['Who Ordinal Scale'] < 6)
                                             & (meta[
                                                    'Patient Location'] == 'ICU')),
         'Disease_severity'] = 'critical'
meta.loc[(meta['Who Ordinal Scale'] < 6)
         & (meta['Patient Location'] != 'ICU'),
         'Disease_severity'] = 'non-critical'
# meta.loc[meta['Who Ordinal Scale'] >= 6, 'Disease_severity'] = 'critical'
# meta.loc[meta['Who Ordinal Scale'] < 6, 'Disease_severity'] = 'non-critical'

tcr_per_virus = tcr_per_virus.merge(meta[['patient_id', 'TP',
                                          'Disease_severity']],
                                    on=['patient_id', 'TP'], how='left')
# occurs only once in metadata, for TP 1 but not TP 2
tcr_per_virus.loc[tcr_per_virus['patient_id'] == 'INCOV129',
                  'Disease_severity'] = 'non-critical'
tcr_per_virus['Disease_severity'].fillna('healthy', inplace=True)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_patients = tcr_per_virus[
    tcr_per_virus['Disease_severity'] != 'healthy'].copy()
all_patients = all_patients.merge(meta_long, on=['patient_id', 'TP'])
# c+ICU 55, n 199

# unify disease severity per patient with the most severe
# c+ICU 67, n 187
all_patients.sort_values(by=['patient_id', 'Disease_severity'], inplace=True)
all_patients['Disease_severity'] = \
    all_patients.groupby('patient_id')['Disease_severity'].transform('first')

all_patients = all_patients[all_patients['week'] == 1]
all_patients.sort_values(by=['patient_id', 'total_count'],
                         inplace=True, ascending=False)
all_patients.drop_duplicates(subset=['patient_id'], inplace=True)

critical = all_patients[all_patients['Disease_severity'] == 'critical']
noncritical = all_patients[
    all_patients['Disease_severity'] == 'non-critical']

# test for differences in the original data
alpha = 0.05
time = 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# differences within groups at TP 1 between SC2-unique and CoV-common TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MU at TP 1 unique vs common
parameters1 = [('Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
                'TCR repertoire depth'),
               ('N_CoV_TCRs', 'N_SARS-CoV-2_only_TCRs', 'N_TCRs'),
               ('N_CoV_Eps', 'N_SARS-CoV-2_only_Eps', 'Diversity'),
               # ('CoV-common_response_redundancy',
               #  'SC2-unique_response_redundancy', 'Redundancy'),
               ('CoV-common_repertoire_breadth',
                'SC2-unique_repertoire_breadth', 'TCR repertoire breadth')]

for parameter in parameters1:
    filter_tp = (all_patients['week'] == time)
    for severity in all_patients['Disease_severity'].unique().tolist():
        data1 = all_patients.loc[
            (all_patients['Disease_severity'] == severity)
            & filter_tp, parameter[0]]
        data2 = all_patients.loc[
            (all_patients['Disease_severity'] == severity)
            & filter_tp, parameter[1]]
        # compare samples
        stat, p = mannwhitneyu(data1, data2)
        print('\n', time, parameter)
        print('Statistics=%.3f, p=%.3f' % (stat, p))
        print(f'Medians: CoV={data1.median()}, SC2={data2.median()}')
        # interpret
        p_adj = p * 2
        if p_adj > alpha:
            print(f'Status: {severity}: Bonferroni corrected p is {p_adj}.'
                  f'\nThe same distribution of SC2-unique and '
                  f'CoV-common TCRs within this patient group (fail to reject H0)')
        else:
            print(f'Status: {severity}: Bonferroni corrected p is {p_adj}.'
                  f'\nDifferent distribution of SC2-unique and '
                  f'CoV-common TCRs within this patient group (reject H0)')
        unique = data2.to_list()
        common = data1.to_list()
        n_pairs = 0
        for i in common:
            for j in unique:
                if i > j:
                    n_pairs += 1
        effect = n_pairs / (len(unique) * len(common))
        print(f'Effect size is {effect}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure SM1.1 like Figure 3: SC2 vs CoV in critical & non-critical (scTCRs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_plot = all_patients[all_patients['week'] == 1].copy()

# TCR metric to plot
parameter = ('Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
             'TCR repertoire depth')

long_df = data_plot.melt(
    id_vars=['patient_id', 'TP', 'Disease_severity'],
    value_vars=[parameter[0], parameter[1]],
    var_name='Specificity', value_name=f'{parameter[2]} at week1')
# rename & reorder for the plot
long_df.loc[long_df['Specificity'] == parameter[0],
            'Specificity'] = 'CoV-common TCRs'
long_df.loc[long_df['Specificity'] == parameter[1],
            'Specificity'] = 'SC2-unique TCRs'
long_df.sort_values(by=['Disease_severity', 'Specificity'],
                    inplace=True, ascending=False)

# plotting params
hue_order = ['CoV-common TCRs', 'SC2-unique TCRs']
palette = {'CoV-common TCRs': 'tab:purple', 'SC2-unique TCRs': 'tab:green'}

np.random.seed(137645)
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(9.9, 6))
ax = sns.stripplot(data=long_df, x='Disease_severity',
                   y=f'{parameter[2]} at week1',
                   hue='Specificity', hue_order=hue_order, palette=palette,
                   linewidth=1, dodge=True, jitter=True)
ax = sns.boxplot(data=long_df, x='Disease_severity',
                 y=f'{parameter[2]} at week1',
                 hue='Specificity', hue_order=hue_order, palette=palette,
                 flierprops={"marker": ""}, showmeans=True,
                 meanprops={'marker': '*', 'markerfacecolor': 'white',
                            'markeredgecolor': 'black'})
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[0:2], labels[0:2], frameon=False)
ax.set_xlabel('\nCOVID-19 severity')
ax.set_xlabel(f'\n{parameter[2]} at week1')
ax.tick_params(labelsize=10)

# stat annotation
x1, x2 = 0 - 0.2, 0 + 0.2
x3, x4 = 1 - 0.2, 1 + 0.2
y, h, col = long_df[
                f'{parameter[2]} at week1'].max() / 2 + 0.015, 0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
plt.text((x1 + x2) * .5, (y + h) * 1.008, 'p < 0.001', ha='center',
         va='bottom', color=col)
plt.plot([x3, x3, x4, x4], [y, y + h, y + h, y], lw=1.5, c=col)
plt.text((x3 + x4) * .5, (y + h) * 1.008, 'p < 0.01', ha='center',
         va='bottom', color=col)
# plt.show()
plt.savefig(
    f'{folder_out}/FigS_{parameter[2]}_week1_CoV_vs_SC2_perSeverity_1.jpg',
    dpi=600, bbox_inches='tight')
plt.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# comparison of different parameters btwn critical vs non-critical at week 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters2 = ['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
               'total_count', '%unique_TCRs']

for parameter in parameters2:
    data1 = critical.loc[critical['week'] == time, parameter]
    data2 = noncritical.loc[noncritical['week'] == time, parameter]
    # compare samples
    stat, p = mannwhitneyu(data1, data2)
    n1 = len(data1)
    val1 = data1[data1 > 0].shape[0]
    prcnt1 = 100 * val1 / n1
    n2 = len(data2)
    val2 = data2[data2 > 0].shape[0]
    prcnt2 = 100 * val2 / n2
    print('\nStatistics=%.3f, p=%.3f' % (stat, p))
    print(f'Medians: critical={data1.median()} '
          f'(n={n1}, non-zero: {val1}={prcnt1}%), '
          f'non-critical={data2.median()}'
          f' (n={n2}, non-zero: {val2}={prcnt2}%)')
    # interpret
    p_adj = p * 3
    if p_adj > alpha:
        print(f'TP={time}, parameter={parameter}: '
              f'Bonferroni corrected p is {p_adj}.'
              f'\nThe same distribution between critical and non-critical'
              f' patients (fail to reject H0)')
    else:
        print(f'TP={time}, parameter={parameter}: '
              f'Bonferroni corrected p is {p_adj}.'
              f'\nDifferent distribution between critical and non-critical'
              f' patients (reject H0)')
    noncrit = data2.to_list()
    crit = data1.to_list()
    n_pairs = 0
    for i in crit:
        for j in noncrit:
            if i > j:
                n_pairs += 1
    effect = n_pairs / (len(noncrit) * len(crit))
    print(f'Effect size is {effect}')
