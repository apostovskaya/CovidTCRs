import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

folder = './data'
folder_out = './results'
all_data_scaled_week_max = pd.read_csv(
    f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t')

parameter = '%unique_TCRs'

# Spearman correlation
crit = all_data_scaled_week_max.loc[
    (all_data_scaled_week_max['critical_disease'] == 'Yes')]
notcrit = all_data_scaled_week_max.loc[
    (all_data_scaled_week_max['critical_disease'] == 'No')]
print(f'\n{parameter}')
print(f'Critical: {spearmanr(crit[parameter], crit["week"])}')
print(f'Not critical: {spearmanr(notcrit[parameter], notcrit["week"])}')

# compare samples before and during recovery
alpha = 0.05
data1 = notcrit.loc[notcrit['week'] > 4, parameter]
data2 = notcrit.loc[notcrit['week'] <= 4, parameter]
stat, p = mannwhitneyu(data1, data2)
print('\nStatistics=%.3f, p=%.4f' % (stat, p))
# interpret
if p > alpha:
    print(f'The same distribution of {parameter} between '
          f'active and recovered non-critical patients (fail to reject H0)')
else:
    print(f'Different distribution of {parameter} between '
          f'active and recovered non-critical patients (reject H0)')
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
# Figure 6: a trend line
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_data_scaled_week_max.loc[
    all_data_scaled_week_max['critical_disease'] == 'Yes',
    'critical_disease'] = 'Critical patients'
all_data_scaled_week_max.loc[
    all_data_scaled_week_max['critical_disease'] == 'No',
    'critical_disease'] = 'Non-critical patients'
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
ax = sns.lineplot(data=all_data_scaled_week_max, x='week', y=parameter,
                  markers=True, dashes=False, style='critical_disease',
                  estimator=np.median, ci=95, n_boot=1000, seed=None,
                  sort=True, err_style='band',
                  hue_order=['Critical patients', 'Non-critical patients'],
                  hue='critical_disease',
                  palette={'Non-critical patients': 'tab:blue',
                           'Critical patients': 'tab:orange'})
# y min and max
ymin, ymax = ax.get_ylim()
plt.vlines(x=4, ymin=ymin, ymax=ymax,
           colors='tab:blue', ls='--', lw=0.75,
           label='Recovery of non-critical patients')
plt.ylabel('% unique TCRs')
plt.legend(frameon=False)
ax.tick_params(labelsize=10)
plt.tight_layout()
plt.savefig(f'{folder_out}/Fig6_TrendDots_{parameter}.jpg', dpi=600)
plt.close()
