import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
from scipy.stats import mannwhitneyu

folder = './data'
folder_out = './results'
all_data_scaled_week_max = pd.read_csv(
    f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t')
week2 = all_data_scaled_week_max[all_data_scaled_week_max['week'] == 2].copy()

alpha = 0.05
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
long_df = week2.melt(id_vars=['patient_id', 'week', 'critical_disease'],
                     value_vars=['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs'],
                     var_name='TCR specificity',
                     value_name='TCR fraction at week2')
# rename & reorder for the plot
long_df.rename(columns={'critical_disease': 'COVID-19 severity'}, inplace=True)
long_df.loc[long_df['COVID-19 severity'] == 'No',
            'COVID-19 severity'] = 'Non-critical patients'
long_df.loc[long_df['COVID-19 severity'] == 'Yes',
            'COVID-19 severity'] = 'Critical patients'
long_df.loc[long_df['TCR specificity'] == 'Frac_CoV_TCRs',
            'TCR specificity'] = 'CoV-common'
long_df.loc[long_df['TCR specificity'] == 'Frac_SARS-CoV-2_only_TCRs',
            'TCR specificity'] = 'SC2-unique'
long_df.sort_values(by=['TCR specificity', 'COVID-19 severity'],
                    inplace=True, ascending=False)

# plotting params
hue_order = ['Non-critical patients', 'Critical patients']
palette = {'Critical patients': 'tab:orange',
           'Non-critical patients': 'tab:blue'}

np.random.seed(137645)
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(9.9, 6))
ax = sns.stripplot(data=long_df, x='TCR specificity', y='TCR fraction at week2',
                   hue='COVID-19 severity', hue_order=hue_order, palette=palette,
                   linewidth=1, dodge=True, jitter=True)
ax = sns.boxplot(data=long_df, x='TCR specificity', y='TCR fraction at week2',
                 hue='COVID-19 severity', hue_order=hue_order, palette=palette,
                 flierprops={"marker": ""}, showmeans=True,
                 meanprops={'marker': '*', 'markerfacecolor': 'white',
                            'markeredgecolor': 'black'})
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[0:2], labels[0:2], frameon=False)
ax.set_xlabel('\nTCR specificity')
ax.tick_params(labelsize=10)
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
# stat annotation
x1, x2 = 0 - 0.2, 0 + 0.2
x3, x4 = 1 - 0.2, 1 + 0.2
y, h, col = long_df['TCR fraction at week2'].max() + 0.1, 0.03, 'k'
y2 = y/2 * 1.05
plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
plt.text((x1 + x2) * .5, (y + h) * 1.008, f"p < {alpha}", ha='center',
         va='bottom', color=col)
plt.plot([x3, x3, x4, x4], [y2, y2 + h, y2 + h, y2], lw=1.5, c=col)
plt.text((x3 + x4) * .5, (y2 + h) * 1.008, 'ns', ha='center',
         va='bottom', color=col)
# plt.show()
plt.savefig(f'{folder_out}/Fig4_TCRfrac_w2_crit_vs_noncrit_perSpec.jpg',
            dpi=600, bbox_inches='tight')
plt.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure S4a and S4b -> replot to have them on one plot and with A, B
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for parameter in ['total_count', '%unique_TCRs']:
    data1 = week2.loc[week2['critical_disease'] == 'Yes', parameter]
    data2 = week2.loc[week2['critical_disease'] == 'No', parameter]
    # compare samples
    stat, p = mannwhitneyu(data1, data2)
    print('\nStatistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    if p > alpha:
        print(f'The same distribution of {parameter} '
              f'in critical and non-critical patients (fail to reject H0)')
        # plot
        mpl.rcParams['font.family'] = 'Arial'
        fig = plt.figure(figsize=(6.6, 5))
        ax = sns.boxplot(x="critical_disease", y=parameter, data=week2,
                         palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                         flierprops={"marker": ""}, showmeans=True,
                         meanprops={'marker': '*',
                                    'markerfacecolor': 'white',
                                    'markeredgecolor': 'black'})
        ax = sns.stripplot(x="critical_disease", y=parameter, data=week2,
                           palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                           dodge=True, linewidth=1)
        ax.set(xlabel='COVID-19 severity', ylabel=f'{parameter} at week2')
        ax.set_xticklabels(['Non-critical', 'Critical'])
        ax.tick_params(labelsize=10)
        plt.tight_layout()
        # statistical annotation (M-U in stats section)
        x1, x2 = 0, 1
        y, h, col = week2[parameter].max() + 0.15, 0.1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        plt.text((x1 + x2) * .5, y + h, "ns", ha='center',
                 va='bottom', color=col)
        plt.savefig(f'{folder_out}/Supplementary/'
                    f'FigS4b_box_{parameter}_w2_perGroup.jpg', dpi=600)
        plt.close()
    else:
        print(f'Different distribution of {parameter} '
              f'in critical and non-critical patients (reject H0)')
        noncrit = data2.to_list()
        crit = data1.to_list()
        n_pairs = 0
        for i in crit:
            for j in noncrit:
                if i > j:
                    n_pairs += 1
        effect = n_pairs / (len(noncrit) * len(crit))
        print(f'Effect size is {effect}')
        # plot
        mpl.rcParams['font.family'] = 'Arial'
        fig = plt.figure(figsize=(6.6, 5))
        ax = sns.boxplot(x="critical_disease", y=parameter, data=week2,
                         palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                         flierprops={"marker": ""}, showmeans=True,
                         meanprops={'marker': '*',
                                    'markerfacecolor': 'white',
                                    'markeredgecolor': 'black'})
        ax = sns.stripplot(x="critical_disease", y=parameter, data=week2,
                           palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                           dodge=True, linewidth=1)
        ax.set(xlabel='COVID-19 severity', ylabel=f'{parameter} at week2')
        ax.set_xticklabels(['Non-critical', 'Critical'])
        ax.tick_params(labelsize=10)
        plt.tight_layout()
        # statistical annotation (M-U in stats section)
        x1, x2 = 0, 1
        y, h, col = week2[parameter].max() + 20000, 10000, 'k'
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        plt.text((x1 + x2) * .5, y + h, f"p < {alpha}", ha='center',
                 va='bottom', color=col)
        plt.savefig(f'{folder_out}/Supplementary/'
                    f'FigS4a_{parameter}_w2_perGroup.jpg', dpi=600)
        plt.close()
