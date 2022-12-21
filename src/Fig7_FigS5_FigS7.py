import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

folder = './data'
folder_out = './results'
folder_logs = './logs'

if not os.path.exists(folder_logs):
    os.makedirs(folder_logs)
    print(f'\nThe logs folder is created in the {folder_logs}')

all_data_scaled_week_max = pd.read_csv(
    f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t')

orig_stdout = sys.stdout
log_f = open(f'{folder_logs}/Fig7.log', 'w+')
sys.stdout = log_f

# Spearman correlation
print('Correlation between time in weeks and TCR fractions '
      'in critical and non-critical patient groups')
alpha = 0.05
crit = all_data_scaled_week_max.loc[
    (all_data_scaled_week_max['critical_disease'] == 'Yes')]
noncrit = all_data_scaled_week_max.loc[
    (all_data_scaled_week_max['critical_disease'] == 'No')]

for parameter in ['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs']:
    print(f'\n{parameter}')
    print(f'Critical: {spearmanr(crit[parameter], crit["week"])}')
    if spearmanr(crit[parameter], crit["week"])[1] < alpha:
        crit_corr = spearmanr(crit[parameter], crit["week"])
        crit_label = mpatches.Patch(color='tab:green',
                                    label=f'rho={crit_corr[0]:.2f},'
                                          f'p={crit_corr[1]:.2f}')
    else:
        crit_label = mpatches.Patch(color='tab:green',
                                    label='ns')
    print(f'Non-critical: {spearmanr(noncrit[parameter], noncrit["week"])}')
    if spearmanr(noncrit[parameter], noncrit["week"])[1] < alpha:
        noncrit_corr = spearmanr(noncrit[parameter], noncrit["week"])
        noncrit_label = mpatches.Patch(color='tab:purple',
                                       label=f'rho={noncrit_corr[0]:.2f},'
                                             f'p={noncrit_corr[1]:.2f}')
    else:
        noncrit_label = mpatches.Patch(color='tab:purple',
                                       label='ns')
sys.stdout = orig_stdout
log_f.close()
print(f'Log file "Fig7.log" is saved in the {folder_logs} directory')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 7a,b: Plots for longitudinal group dynamics of TCR Fractions
# Fraction of CoV vs SC2 TCRs together on one plot,
# disease severity on separate subplots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_long = all_data_scaled_week_max.melt(
    id_vars=['patient_id', 'week', 'critical_disease'],
    value_vars=['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs'],
    var_name='Specificity', value_name='Fraction')
# rename specificities to match paper vocabulary
data_long.loc[data_long['Specificity'] == 'Frac_CoV_TCRs',
              'Specificity'] = 'CoV-common_TCRs'
data_long.loc[data_long['Specificity'] == 'Frac_SARS-CoV-2_only_TCRs',
              'Specificity'] = 'SC2-unique_TCRs'

mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
# x_estimator=np.mean to g to have lines instead of individual dots
g = sns.lmplot(x='week', y='Fraction', data=data_long, col='critical_disease',
               scatter_kws={'s': 10},
               hue='Specificity', palette={'CoV-common_TCRs': 'tab:purple',
                                           'SC2-unique_TCRs': 'tab:green'},
               ci=95, legend=False, markers=["o", "x"])
(g.set_axis_labels("week", "TCR fraction")
 .tight_layout(w_pad=2))
fig = g.fig
a0 = fig.axes[0]
a0.set_title("Non-critical COVID-19")
a1 = fig.axes[1]
a1.set_title("Critical COVID-19")
legend1 = plt.legend(loc='upper right', title='Specificity', markerscale=2)
legend2 = plt.legend(handles=[crit_label, noncrit_label],
                     loc='lower right', handlelength=0.25)
plt.gca().add_artist(legend1)
plt.savefig(f'{folder_out}/Fig7ab_RegressDots_Fractions_crit_not_week.jpg',
            dpi=600)
plt.close()
print(f'Figure 7 (a, b) is saved in the {folder_out} directory')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS7: trend lines of depth and breadth of SC2-unique and CoV-common TCRreps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS7-a,b: trend lines for fractions (depth) of SC2-unique & CoV-common TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
g = sns.relplot(x='week', y='Fraction', data=data_long, col='critical_disease',
                kind='line',
                markers=True,
                # {'Frac_CoV_TCRs': "x", 'Frac_SARS-CoV-2_only_TCRs': "o"},
                dashes=True, style='Specificity', legend=False,
                estimator=np.median, ci=95, n_boot=1000, seed=None, sort=True,
                err_style='band',
                hue='Specificity', palette={'CoV-common_TCRs': 'tab:purple',
                                            'SC2-unique_TCRs': 'tab:green'})
(g.set_axis_labels("week", "TCR fraction")
 .tight_layout(w_pad=2))
fig = g.fig
a0 = fig.axes[0]
a0.set_title("Non-critical COVID-19")
a1 = fig.axes[1]
a1.set_title("Critical COVID-19")
plt.savefig(
    f'{folder_out}/Supplementary/FigS7ab_TrendBand_RepDepth_crit_not_week.jpg',
    dpi=600)
plt.close()

# FigS7-c, d: trend lines for % (breadth) of SC2-unique and CoV-common TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_long2 = all_data_scaled_week_max.melt(
    id_vars=['patient_id', 'week', 'critical_disease'],
    value_vars=['%CoV-common_TCRs', '%SC2-unique_TCRs'],
    var_name='Specificity', value_name='Percentage')
data_long2.loc[data_long2[
                   'Specificity'] == '%CoV-common_TCRs',
               'Specificity'] = 'CoV-common_TCRs'
data_long2.loc[data_long2[
                   'Specificity'] == '%SC2-unique_TCRs',
               'Specificity'] = 'SC2-unique_TCRs'

mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
g = sns.relplot(x='week', y='Percentage', data=data_long2,
                col='critical_disease', kind='line',
                markers=True,
                # {'Frac_CoV_TCRs': "x", 'Frac_SARS-CoV-2_only_TCRs': "o"},
                dashes=True, style='Specificity', legend=True,
                estimator=np.median, ci=95, n_boot=1000, seed=None, sort=True,
                err_style='band',
                hue='Specificity', palette={'CoV-common_TCRs': 'tab:purple',
                                            'SC2-unique_TCRs': 'tab:green'})
(g.set_axis_labels("week", "% of specific TCRs")
 .tight_layout(w_pad=2))
fig = g.fig
a0 = fig.axes[0]
a0.set_title("Non-critical COVID-19")
a1 = fig.axes[1]
a1.set_title("Critical COVID-19")
leg = g._legend
leg.set_bbox_to_anchor([0.8, 0.8])
plt.savefig(f'{folder_out}/Supplementary'
            f'/FigS7cd_TrendBand_RepBreadth_crit_not_week.jpg', dpi=600)
plt.close()
print(f'Figure S7 (a-d) is saved in the {folder_out}/Supplementary directory')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS5-a,b: personal dynamics CoV vs SC2 for each patient on their own plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# separating a dataset with only patients who have multiple data points
long_data_only_weekMax = all_data_scaled_week_max[
    all_data_scaled_week_max.duplicated(subset='patient_id',
                                        keep=False)].copy()

# log2(FC+1) (frequency + 1 -> FC -> log2) between max of a standard week
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# shifting values by 1 to calculate fold change without dealing with 0s
long_data_only_weekMax['Frac_CoV_TCRs_+1'] = long_data_only_weekMax[
                                                 'Frac_CoV_TCRs'] + 1
long_data_only_weekMax['Frac_SARS-CoV-2_only_TCRs_+1'] = \
    long_data_only_weekMax['Frac_SARS-CoV-2_only_TCRs'] + 1

# calculating fold change
for parameter in ['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs',
                  'Frac_CoV_TCRs_+1', 'Frac_SARS-CoV-2_only_TCRs_+1']:
    long_data_only_weekMax[f'weekFC_{parameter}'] = long_data_only_weekMax[
        parameter] \
        .div(
        long_data_only_weekMax.groupby('patient_id')[f'{parameter}'].shift(1))

# filling in NANs for the first time points with 1 (no change)
long_data_only_weekMax[
    ['weekFC_Frac_CoV_TCRs_+1', 'weekFC_Frac_SARS-CoV-2_only_TCRs_+1']] = \
    long_data_only_weekMax[
        ['weekFC_Frac_CoV_TCRs_+1',
         'weekFC_Frac_SARS-CoV-2_only_TCRs_+1']].fillna(1)

# transforming to log2 fold change
for fc_col in ['weekFC_Frac_CoV_TCRs_+1',
               'weekFC_Frac_SARS-CoV-2_only_TCRs_+1']:
    long_data_only_weekMax[f'log2_{fc_col}'] = np.log2(
        long_data_only_weekMax[f'{fc_col}'])

# transforming the dataset for plotting
long_data_only_weekMax_long_logFC = long_data_only_weekMax.melt(
    id_vars=['patient_id', 'week', 'critical_disease'],
    value_vars=['log2_weekFC_Frac_CoV_TCRs_+1',
                'log2_weekFC_Frac_SARS-CoV-2_only_TCRs_+1'],
    var_name='Specificity', value_name='log2FC_Fraction')
long_data_only_weekMax_long_logFC.loc[long_data_only_weekMax_long_logFC[
                                          'Specificity'] == 'log2_weekFC_Frac_SARS-CoV-2_only_TCRs_+1',
                                      'Specificity'] = 'SC2-unique'
long_data_only_weekMax_long_logFC.loc[long_data_only_weekMax_long_logFC[
                                          'Specificity'] == 'log2_weekFC_Frac_CoV_TCRs_+1',
                                      'Specificity'] = 'CoV-common'
long_data_only_weekMax_long_logFC.loc[
        long_data_only_weekMax_long_logFC['critical_disease'] == 'Yes',
        'critical_disease'] = 'Critical'
long_data_only_weekMax_long_logFC.loc[
        long_data_only_weekMax_long_logFC['critical_disease'] == 'No',
        'critical_disease'] = 'Non-critical'
# plotting
for severity in long_data_only_weekMax_long_logFC['critical_disease'].unique().tolist():
    data_plot = long_data_only_weekMax_long_logFC[
        long_data_only_weekMax_long_logFC['critical_disease'] == severity]
    for parameter in [['CoV-common', 'SC2-unique', 'TCR_log2FC']]:
        mpl.rcParams['font.family'] = 'Arial'
        fig = plt.figure(figsize=(6.6, 5))
        ax = sns.relplot(x='week', y='log2FC_Fraction', data=data_plot,
                         col='patient_id', kind='line',
                         markers=True,
                         dashes=False, style='Specificity', col_wrap=3,
                         hue='Specificity',  # linewidth=3,
                         palette={parameter[0]: 'tab:purple',
                                  parameter[1]: 'tab:green'}, legend=False)
        ax.map(plt.axhline, y=0, color=".7", dashes=(2, 1), zorder=0)
        ax.fig.suptitle(f"{severity} COVID-19")
        plt.subplots_adjust(top=0.925)
        plt.savefig(
            f'{folder_out}/Supplementary/FigS5_{parameter[2]}+1_CoV_vs_SC2_individs_{severity}.jpg',
            bbox_inches='tight', dpi=600)
        plt.close()
print(f'Figure S5 (a-b) is saved in the {folder_out}/Supplementary directory')
