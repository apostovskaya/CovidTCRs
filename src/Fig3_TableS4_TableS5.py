import sys
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
from scipy.stats import mannwhitneyu

folder = './data'
folder_out = './results'
folder_logs = './logs'

if not os.path.exists(folder_logs):
    os.makedirs(folder_logs)
    print(f'\nThe logs folder is created in the {folder_logs}')
# to have outputs in a log file
orig_stdout = sys.stdout
log_f = open(f'{folder_logs}/Fig3_TableS4_TableS5.log', 'w+')
sys.stdout = log_f

all_data_scaled_week_max = pd.read_csv(f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t')
print('Raw statistical values, without multiple testing corrections.')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# differences within groups at week 1 between SC2-unique and CoV-common TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MU at week 1 unique vs common
alpha = 0.05
for severity in ['Yes', 'No']:
    data1 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == severity) &
                                         (all_data_scaled_week_max['week'] == 1), 'Frac_CoV_TCRs']
    data2 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == severity) &
                                         (all_data_scaled_week_max['week'] == 1), 'Frac_SARS-CoV-2_only_TCRs']
    # compare samples
    stat, p = mannwhitneyu(data1, data2)
    print('\nStatistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    if severity == 'Yes':
        p_adj = p * 8
    else:
        p_adj = p * 4
    if p_adj > alpha:
        print(f'Critical={severity}: Bonferroni corrected p is {p_adj}.'
              f'\nThe same distribution of SC2-unique and '
              f'CoV-common TCRs within this patient group (fail to reject H0)')
    else:
        print(f'Critical={severity}: Bonferroni corrected p is {p_adj}.'
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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Figure 3: plot critical and not separately side-by-side
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data_long = all_data_scaled_week_max.melt(id_vars=['patient_id', 'week',
                                                       'critical_disease'],
                                              value_vars=['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs'],
                                              var_name='Specificity', value_name='Fraction')
    for parameter in [['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs', 'TCRs']]:
        data_plot = data_long[(data_long['week'] == 1) & (data_long['critical_disease'] == severity)].copy()
        data_plot.loc[data_plot['Specificity'] == 'Frac_CoV_TCRs', 'Specificity'] = 'CoV-common'
        data_plot.loc[data_plot['Specificity'] == 'Frac_SARS-CoV-2_only_TCRs', 'Specificity'] = 'SC2-unique'
        mpl.rcParams['font.family'] = 'Arial'
        fig = plt.figure(figsize=(6.6, 5))
        ax = sns.boxplot(x="Specificity", y='Fraction', data=data_plot,
                         palette={'CoV-common': 'tab:purple', 'SC2-unique': 'tab:green'}, showmeans=True,
                         meanprops={'marker': '*', 'markerfacecolor': 'white', 'markeredgecolor': 'black'})
        ax = sns.stripplot(x="Specificity", y='Fraction', data=data_plot,
                           palette={'CoV-common': 'tab:purple', 'SC2-unique': 'tab:green'}, dodge=True, linewidth=1)
        ax.set(xlabel='Recognized epitopes', ylabel=f'TCR fractions at week1')
        ax.set_ylim(top=0.9)
        plt.tight_layout()
        # statistical annotation (M-U in stats section)
        if p_adj < alpha:
            x1, x2 = 0, 1
            y, h, col = data_plot[data_plot['week'] == 1]['Fraction'].max() + 0.05, 0.1, 'k'
            plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
            plt.text((x1 + x2) * .5, (y + h) * 1.008, f"p < {alpha}", ha='center', va='bottom', color=col)
            ax.tick_params(labelsize=10)
            ax.yaxis.set_major_locator(MultipleLocator(0.2))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            plt.title(f'Critical COVID-19')
            plt.savefig(f'{folder_out}/Fig3b_box{parameter[2]}_w1_{severity}_CoV_vs_SC2.jpg',
                        bbox_inches='tight', dpi=600)
            plt.close()
        else:
            x1, x2 = 0, 1
            y, h, col = data_plot[data_plot['week'] == 1]['Fraction'].max() + 0.05, 0.05, 'k'
            plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
            plt.text((x1 + x2) * .5, (y + h) * 1.008, "ns", ha='center', va='bottom', color=col)
            ax.tick_params(labelsize=10)
            ax.yaxis.set_major_locator(MultipleLocator(0.2))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            plt.title(f'Non-ritical COVID-19')
            plt.savefig(f'{folder_out}/Fig3a_box{parameter[2]}_w1_{severity}_CoV_vs_SC2.jpg',
                        bbox_inches='tight', dpi=600)
            plt.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for supplementary table of all active weeks (Table S4)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table1 = pd.DataFrame({'week': [i for i in range(1, 9)]})
for i in range(1, 9):
    data1 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == 'Yes') &
                                         (all_data_scaled_week_max['week'] == i), 'Frac_CoV_TCRs']
    data2 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == 'Yes') &
                                         (all_data_scaled_week_max['week'] == i),
                                         'Frac_SARS-CoV-2_only_TCRs']
    data3 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == 'No') &
                                         (all_data_scaled_week_max['week'] == i), 'Frac_CoV_TCRs']
    data4 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == 'No') &
                                         (all_data_scaled_week_max['week'] == i),
                                         'Frac_SARS-CoV-2_only_TCRs']
    # critical
    stat, p = mannwhitneyu(data1, data2)
    # Bonferroni correction for multiple comparisons (8 parameters)
    p = p * 8
    p = round(p, 3)
    # probabilities cannot exceed 1
    if p > 1:
        p = 1
    table1.loc[table1['week'] == i, 'adj_p-value critical'] = p
    # AUC
    assert data1.shape[0] == data2.shape[0]
    n_pairs = 0
    for k in data1.to_list():
        for l in data2.to_list():
            if k > l:
                n_pairs += 1
    effect = n_pairs / (len(data1.to_list()) * len(data2.to_list()))
    effect = round(effect, 3)
    table1.loc[table1['week'] == i, 'AUC critical'] = effect
    table1.loc[table1['week'] == i, 'No. critical'] = int(data1.shape[0])
    # non-critical
    if i < 5:
        stat2, p2 = mannwhitneyu(data3, data4)
        # Bonferroni correction for multiple comparisons (4 parameters)
        p2 = p2 * 4
        p2 = round(p2, 3)
        # probabilities cannot exceed 1
        if p2 > 1:
            p2 = 1
        table1.loc[table1['week'] == i, 'adj_p-value non-critical'] = p2
        # AUC
        assert data3.shape[0] == data4.shape[0]
        n_pairs2 = 0
        for k2 in data3.to_list():
            for l2 in data4.to_list():
                if k2 > l2:
                    n_pairs2 += 1
        effect2 = n_pairs2 / (len(data3.to_list()) * len(data4.to_list()))
        effect2 = round(effect2, 3)
        table1.loc[table1['week'] == i, 'AUC non-critical'] = effect2
        table1.loc[table1['week'] == i, 'No. non-critical'] = int(data3.shape[0])
    else:
        table1.loc[table1['week'] == i, 'adj_p-value non-critical'] = 'NA'
        table1.loc[table1['week'] == i, 'AUC non-critical'] = 'NA'
        table1.loc[table1['week'] == i, 'No. non-critical'] = 'NA'
table1.to_excel(f'{folder_out}/Supplementary/TableS4_SC2-vs-CoV_allWeeks.xlsx', index=False)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for supplementary Table S5
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table2 = pd.DataFrame(columns=['week', 'parameter', 'adj_p-value', 'AUC'])

week1 = all_data_scaled_week_max[all_data_scaled_week_max['week'] == 1]
for parameter in ['total_count', '%unique_TCRs', 'Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs']:
    stat, p = mannwhitneyu(week1.loc[week1['critical_disease'] == 'Yes', parameter],
                           week1.loc[week1['critical_disease'] == 'No', parameter])
    print(f'\n{parameter} at week1')
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # Bonferroni correction for multiple comparisons (4 parameters)
    p = p*4
    p = round(p, 3)
    crit_tcrs = week1.loc[week1['critical_disease'] == 'Yes', parameter].to_list()
    non_crit_tcrs = week1.loc[week1['critical_disease'] == 'No', parameter].to_list()
    n_pairs = 0
    for i in non_crit_tcrs:
        for j in crit_tcrs:
            if i > j:
                n_pairs += 1
    effect = n_pairs / (len(crit_tcrs) * len(non_crit_tcrs))
    effect = round(effect, 3)
    print(f'Effect size for {parameter} is {effect}')
    new_row = pd.DataFrame(data=[{'week': 1, 'parameter': parameter,
                                  'adj_p-value': p, 'AUC': effect}])
    table2 = pd.concat([table2, new_row], ignore_index=True)

week2 = all_data_scaled_week_max[all_data_scaled_week_max['week'] == 2]
for parameter in ['N_CoV_Eps', 'scaled_N_CoV-common_Eps',
                  'N_SARS-CoV-2_only_Eps', 'scaled_N_SC2-unique_Eps']:
    stat, p = mannwhitneyu(week2.loc[week2['critical_disease'] == 'Yes', parameter],
                           week2.loc[week2['critical_disease'] == 'No', parameter])
    print(f'\n{parameter} at week2')
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # Bonferroni correction for multiple comparisons (2 parameters - CoV, SC2)
    p = p * 2
    p = round(p, 3)
    crit_tcrs = week2.loc[week2['critical_disease'] == 'Yes', parameter].to_list()
    non_crit_tcrs = week2.loc[week2['critical_disease'] == 'No', parameter].to_list()
    n_pairs = 0
    for i in non_crit_tcrs:
        for j in crit_tcrs:
            if i > j:
                n_pairs += 1
    effect = n_pairs / (len(crit_tcrs) * len(non_crit_tcrs))
    effect = round(effect, 3)
    print(f'Effect size for {parameter} is {effect}')
    new_row = pd.DataFrame([{'week': 2, 'parameter': parameter,
                             'adj_p-value': p, 'AUC': effect}])
    table2 = pd.concat([table2, new_row], ignore_index=True)
table2['parameter'] = ['total_TCRcount', '%unique_TCRs', 'Freq_CoV-common_TCRs', 'Freq_SC2-unique_TCRs',
                       'N_CoV-common_Eps', 'normN_CoV-common_Eps', 'N_SC2-unique_Eps', 'normN_SC2-unique_Eps']
# probabilities cannot exceed 1
table2.loc[table2['adj_p-value'] > 1, 'adj_p-value'] = 1
table2.to_excel(f'{folder_out}/Supplementary/TableS5_crit-vs-not_Weeks1-2.xlsx', index=False)
print('\nSignificantly different between critical and non-critical patients'
      ' after Bonferroni correction for paired comparisons (TableS5):')
print(table2.loc[table2['adj_p-value'] < alpha, ['parameter', 'week']])

sys.stdout = orig_stdout
log_f.close()

print(f'Figure 3 is added to the {folder_out}, '
      f'Table S4 and Table S5 are saved in the {folder_out}/Supplementary')
