import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.stats import spearmanr, pearsonr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Epitope metaDATA assembly (done)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PMe epitope data: common CV vs SARS-CoV-2 eps
folder_in = './data'
folder_out = './results'

# Epitope metaDATA assembly
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_data = pd.read_csv(f'{folder_in}/epitope_species.txt', sep='\t')
# epitope protein origin data
prot_data = pd.read_csv(f'{folder_in}/epitope_proteins.txt', sep='\t')
prot_data2 = pd.read_csv(f'{folder_in}/table_wproteins.txt', sep='&',
                         skipinitialspace=True,
                         names=['epitope', 'protein', 'No. training TCRs',
                                'Balanced accuracy', 'AUC ROC', 'AUC PR'])
prot_data2['AUC PR'] = prot_data2['AUC PR'].str.strip(' \\\\')
for col in prot_data2.columns.tolist():
    try:
        prot_data2[col] = prot_data2[col].str.strip()
    except AttributeError:
        continue

prot_data.rename(columns={'protein': 'SARS2protein'}, inplace=True)
prot_data[['protein_info', 'description', 'other']] = prot_data[
    'description'].str.split(pat='; ', expand=True)
prot_data = prot_data[
    ['epitope', 'protein_info', 'description', 'other', 'SARS2protein']]
prot_data['protein_info'] = prot_data['protein_info'].str.lower()
# merging
ep_data = pd.merge(species_data, prot_data, on='epitope', how='left')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data to make supp table about models (Supp_Table_3 paper-v.2022)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ep_model_data = pd.merge(ep_data, prot_data2, on='epitope', how='left')
ep_model_data = ep_model_data[['epitope', 'species', 'protein',
                               'No. training TCRs', 'Balanced accuracy',
                               'AUC ROC', 'AUC PR']]
# custom sort with regular sort
custom_sort = {'M (ORF5)': 3, 'N (ORF9)': 7, 'ORF1ab': 0, 'ORF3a': 2,
               'ORF6': 4, 'ORF7a': 5, 'ORF8': 6, 'S (ORF2)': 1}
ep_model_data['sort'] = ep_model_data['protein'].map(custom_sort)
ep_model_data.sort_values(by=['sort', 'species', 'epitope'], inplace=True)
# custom sort only
# ep_model_data.sort_values(by=['protein', 'species', 'epitope'], key=lambda x: x.map(custom_sort), inplace=True)
ep_model_data.loc[ep_model_data['species'] == 1, 'Uniqueness'] = 'SC2-unique'
ep_model_data.loc[ep_model_data['species'] > 1, 'Uniqueness'] = 'CoV-common'
ep_model_data.reset_index(drop=True, inplace=True)
ep_model_data.rename(columns={'epitope': 'Epitope', 'species': 'No. Species',
                              'protein': 'Protein'}, inplace=True)
ep_model_data = ep_model_data[
    ['Epitope', 'Uniqueness', 'No. Species', 'Protein', 'No. training TCRs',
     'Balanced accuracy', 'AUC ROC', 'AUC PR']]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# epitope data only - Do I need it for the paper? seems that not
# ep_data = pd.merge(ep_data, prot_data2[['epitope', 'protein']], on='epitope', how='left')
# ep_data = ep_data[['epitope', 'species', 'protein', 'protein_info', 'SARS2protein', 'description', 'other']]

# ep_data.to_csv(f'{folder_out}/epitope_data.txt', sep='\t', index=False)
ep_model_data.to_csv(
    f'{folder_out}/Supplementary/TableS3_epitope_model_data.tsv', sep='\t',
    index=False)

# plot this data to remake PMe preprint plot
# ep_model_data = pd.read_csv(f'{folder_out}/Supplementary/TableS3_epitope_model_data.tsv', sep='\t')
ep_model_data.rename(columns={'Uniqueness': 'Specificity'}, inplace=True)
# ep_data_toplot = ep_model_data.groupby('Protein')['No. Species']
ax = sns.boxplot(x="Protein", y='No. Species', data=ep_model_data,
                 color=".95", saturation=0.1, showmeans=False)
# ax = sns.violinplot(x="Protein", y="No. Species", data=ep_model_data,
#                     inner='quartile', cut=0, color=".95", saturation=0.1)
ax = sns.swarmplot(x="Protein", y='No. Species', data=ep_model_data,
                   hue='Specificity',
                   palette={'SC2-unique': 'tab:green',
                            'CoV-common': 'tab:purple'},
                   linewidth=0.5, size=3.5)  # color=".3", jitter=0.35,
ax.set(xlabel='Protein of origin', ylabel='Nidovirales species matches')
plt.tight_layout()
plt.savefig(f'{folder_out}/Fig1_NSpecies-vs-Protein.png')
plt.close()

# additional supplementary figure to show that there is no correlation between
# training size and AUC
# ep_model_data[cols] = ep_model_data[cols].apply(pd.to_numeric, errors='coerce')
cols = ['AUC ROC', 'Balanced accuracy', 'AUC PR']
ep_model_data_long = ep_model_data.melt(id_vars=['No. training TCRs'],
                                        value_vars=cols,
                                        var_name='Performance evaluation metric',
                                        value_name='Value')
ep_model_data_long['Value'] = \
    ep_model_data_long['Value'].str.split(' ', expand=True)[0]
ep_model_data_long['Value'] = pd.to_numeric(ep_model_data_long['Value'],
                                            errors='coerce')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure S1 ab
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.family'] = 'Arial'
fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(13, 5))
sns.histplot(x='No. training TCRs', data=ep_model_data,
             color='lightgrey', log_scale=True, ax=ax1)
# plt.axhline(y=1, ls='--', lw=0.75, color='black')
ax1.vlines(x=30, ymin=0, ymax=10.5, ls='--', lw=1.5, colors='red',
           label='Required minimum = 30 TCRs')
ax1.vlines(x=ep_model_data['No. training TCRs'].median(),
           ymin=0, ymax=10.5, ls='--', lw=1.5, colors='black',
           label=f'Median = {int(ep_model_data["No. training TCRs"].median())} '
                 f'TCRs')
ax1.legend(frameon=False, title='Positive training data size')
ax1.set_xlabel('No. positive training TCRs')
ax1.set_title('Distribution of the positive training data size\n')
# ax1.text(-0.1, 1, 'A', transform=ax.transAxes, size=20, weight='bold')
ax1.annotate('A', (-0.15, 1.1), xycoords='axes fraction',
             size=14)

sns.lineplot(x='No. training TCRs', y='Value', data=ep_model_data_long,
             hue='Performance evaluation metric', hue_order=cols,
             palette={cols[1]: 'darkgrey', cols[0]: 'grey',
                      cols[2]: 'lightgrey'}, legend=True, ax=ax2)
ax2.set(xscale='log')
ax2.set_xlabel('No. positive training TCRs')
ax2.set_ylabel('Metric value')
ax2.set_title('Performance evaluation metrics of the models\n'
              'with different size of the positive training data')
# ax2.text(2.25, 1, 'B', transform=ax.transAxes, size=20, weight='bold')
ax2.annotate('B', (1.05, 1.1), xycoords='axes fraction',
             size=14)
# remove legend title and frame
plt.gca().legend().set_title('')
plt.legend(frameon=False)
plt.tight_layout()
# plt.show()
plt.savefig(f'{folder_out}/Supplementary/FigS1ab_training_size.jpg',
            dpi=600, bbox_inches="tight")
plt.close()

# one
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
ax = sns.lineplot(x='No. training TCRs', y='Value', data=ep_model_data_long,
                  # scatter_kws={'s': 20}, x_ci='sd',
                  # markers=['v', 'x', 'o'],
                  hue='Performance evaluation metric', hue_order=cols,
                  palette={cols[1]: 'darkgrey', cols[0]: 'grey',
                           cols[2]: 'lightgrey'}, legend=True)
ax.set(xscale='log')
ax.set_ylabel('Value of the metric')
ax.set_title('Performance evaluation metrics \nof the models '
             'with different size of the training data')
# remove legend title
plt.gca().legend().set_title('')
plt.tight_layout()
# plt.show()
plt.savefig(f'{folder_out}/Supplementary/FigS8_training_size_vs_metric.jpg',
            dpi=600)
plt.close()

for metric in cols:
    stat_data = ep_model_data_long[
        ep_model_data_long['Performance evaluation metric'] == metric]
    metric_stats = pearsonr(stat_data['No. training TCRs'],
                             stat_data['Value'])
    print(f'{metric}: {metric_stats}')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# spearmanr (evaluates the monotonic relationship, based on the ranked values)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC ROC: SpearmanrResult(correlation=0.34882899750318735,
# pvalue=0.016257840610589745)
# Balanced accuracy: SpearmanrResult(correlation=0.09601255882283027,
# pvalue=0.5208795668819042)
# AUC PR: SpearmanrResult(correlation=0.22178178280193256,
# pvalue=0.1340694145395434)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pearsonr (evaluates the linear relationship between two continuous variables)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC ROC: (0.24615439752816543, 0.0953393246635227)
# Balanced accuracy: (0.07784208283297232, 0.6030088227976087)
# AUC PR: (0.14560583209476713, 0.32878661213940186)
#
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
sns.histplot(x='No. training TCRs', data=ep_model_data, color='lightgrey',
             log_scale=True)
# plt.axhline(y=1, ls='--', lw=0.75, color='black')
plt.vlines(x=30, ymin=0, ymax=10.5, ls='--', lw=1.5, colors='red',
           label='Minimum required size of the training data')
plt.vlines(x=ep_model_data['No. training TCRs'].median(),
           ymin=0, ymax=10.5, ls='--', lw=1.5, colors='black',
           label='Median size of the training data')
plt.tight_layout()
plt.legend(frameon=False)
plt.show()
plt.savefig(f'{folder_out}/Supplementary/FigS8a_severity_distribution.jpg',
            dpi=600)
plt.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ep_model_data = pd.read_csv(f'{folder_out}/'
                            f'Supplementary/TableS3_epitope_model_data.tsv',
                            sep='\t')
ep_model_data.rename(columns={'Uniqueness': 'Specificity'}, inplace=True)
# 1 column width = 85mm -> 3.346 inch ~ 3.3
# mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
# ax = plt.axes((0.1,0.1,0.5,0.8))
ax = sns.boxplot(x="Protein", y='No. Species', data=ep_model_data,
                 color=".95", saturation=0.1, showmeans=False)
ax = sns.swarmplot(x="Protein", y='No. Species', data=ep_model_data,
                   hue='Specificity',
                   palette={'SC2-unique': 'tab:green',
                            'CoV-common': 'tab:purple'},
                   linewidth=0.5, size=3.5)
ax.set(xlabel='Protein of origin', ylabel='Nidovirales species matches')
ax.tick_params(axis='x', which='major', direction='out', labelsize=8)
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.legend(title='Specificity', frameon=False)
plt.tight_layout()
plt.savefig(f'{folder_out}/Fig1_NSpecies-vs-Protein.jpg', dpi=600)
plt.close()
