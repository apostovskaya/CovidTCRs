import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

folder = './data'
folder_out = './results'
all_data_scaled_week_max = pd.read_csv(f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t')

all_data_scaled_week_max.rename(columns={'Frac_CoV_TCRs': 'Fraction_CoV-common_TCRs',
                                         'Frac_SARS-CoV-2_only_TCRs': 'Fraction_SC2-unique_TCRs'},
                                inplace=True)
week2 = all_data_scaled_week_max[all_data_scaled_week_max['week'] == 2]

alpha = 0.05
for parameter in ['total_count', '%unique_TCRs',
                  'Fraction_CoV-common_TCRs', 'Fraction_SC2-unique_TCRs']:
    data1 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == 'Yes') &
                                         (all_data_scaled_week_max['week'] == 2), parameter]
    data2 = all_data_scaled_week_max.loc[(all_data_scaled_week_max['critical_disease'] == 'No') &
                                         (all_data_scaled_week_max['week'] == 2), parameter]
    # compare samples
    stat, p = mannwhitneyu(data1, data2)
    print('\nStatistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    if p > alpha:
        print(f'The same distribution of {parameter} in critical and non-critical patients (fail to reject H0)')
        # plot
        mpl.rcParams['font.family'] = 'Arial'
        fig = plt.figure(figsize=(6.6, 5))
        ax = sns.boxplot(x="critical_disease", y=parameter, data=week2,
                         palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                         showmeans=True,
                         meanprops={'marker': '*',
                                    'markerfacecolor': 'white',
                                    'markeredgecolor': 'black'})
        ax = sns.stripplot(x="critical_disease", y=parameter, data=week2,
                           palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                           dodge=True, linewidth=1)
        ax.set(xlabel='COVID-19 severity', ylabel=f'{parameter}_week2')
        ax.set_xticklabels(['Non-critical', 'Critical'])
        ax.tick_params(labelsize=10)
        plt.tight_layout()
        # statistical annotation (M-U in stats section)
        x1, x2 = 0, 1
        y, h, col = week2[
                        parameter].max() + 0.15, 0.1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        plt.text((x1 + x2) * .5, y + h, "ns", ha='center',
                 va='bottom', color=col)
        if (parameter != '%unique_TCRs') & (parameter != 'total_count'):
            plt.savefig(f'{folder_out}/Fig4b_box_{parameter}_per_group_stats_w2.jpg', dpi=600)
            plt.close()
        else:
            plt.savefig(f'{folder_out}/Supplementary/FigS3b_box_{parameter}_per_group_stats_w2.jpg', dpi=600)
            plt.close()
    else:
        print(f'Different distribution of {parameter} in critical and non-critical patients (reject H0)')
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
                         showmeans=True,
                         meanprops={'marker': '*',
                                    'markerfacecolor': 'white',
                                    'markeredgecolor': 'black'})
        ax = sns.stripplot(x="critical_disease", y=parameter, data=week2,
                           palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                           dodge=True, linewidth=1)
        ax.set(xlabel='COVID-19 severity', ylabel=f'{parameter}_week2')
        ax.set_xticklabels(['Non-critical', 'Critical'])
        ax.tick_params(labelsize=10)
        plt.tight_layout()
        # statistical annotation (M-U in stats section)
        x1, x2 = 0, 1
        if parameter == 'total_count':
            y, h, col = week2[
                            parameter].max() + 20000, 10000, 'k'
        else:
            y, h, col = week2[
                            parameter].max() + 0.15, 0.1, 'k'
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        plt.text((x1 + x2) * .5, y + h, f"p < {alpha}", ha='center',
                 va='bottom', color=col)
        if (parameter != '%unique_TCRs') & (parameter != 'total_count'):
            plt.savefig(f'{folder_out}/Fig4a_box_{parameter}_per_group_stats_w2.jpg', dpi=600)
            plt.close()
        else:
            plt.savefig(f'{folder_out}/Supplementary/FigS3a_box_{parameter}_per_group_stats_w2.jpg', dpi=600)
            plt.close()
