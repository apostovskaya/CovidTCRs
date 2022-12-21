import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
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
prot_data2 = pd.read_csv(f'{folder_in}/table_wproteins.txt', sep='&', skipinitialspace=True,
                         names=['epitope', 'protein', 'No. training TCRs', 'Balanced accuracy', 'AUC ROC', 'AUC PR'])
prot_data2['AUC PR'] = prot_data2['AUC PR'].str.strip(' \\\\')
for col in prot_data2.columns.tolist():
    try:
        prot_data2[col] = prot_data2[col].str.strip()
    except AttributeError:
        continue

prot_data.rename(columns={'protein': 'SARS2protein'}, inplace=True)
prot_data[['protein_info', 'description', 'other']] = prot_data['description'].str.split(pat='; ', expand=True)
prot_data = prot_data[['epitope', 'protein_info', 'description', 'other', 'SARS2protein']]
prot_data['protein_info'] = prot_data['protein_info'].str.lower()
# merging
ep_data = pd.merge(species_data, prot_data, on='epitope', how='left')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data to make supp table about models (Supp_Table_3 paper-v.2022)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ep_model_data = pd.merge(ep_data, prot_data2, on='epitope', how='left')
ep_model_data = ep_model_data[['epitope', 'species', 'protein',
                               'No. training TCRs', 'Balanced accuracy', 'AUC ROC', 'AUC PR']]
# custom sort with regular sort
custom_sort = {'M (ORF5)': 3, 'N (ORF9)': 7, 'ORF1ab': 0, 'ORF3a': 2, 'ORF6': 4, 'ORF7a': 5, 'ORF8': 6, 'S (ORF2)': 1}
ep_model_data['sort'] = ep_model_data['protein'].map(custom_sort)
ep_model_data.sort_values(by=['sort', 'species', 'epitope'], inplace=True)
# custom sort only
# ep_model_data.sort_values(by=['protein', 'species', 'epitope'], key=lambda x: x.map(custom_sort), inplace=True)
ep_model_data.loc[ep_model_data['species'] == 1, 'Uniqueness'] = 'SC2-unique'
ep_model_data.loc[ep_model_data['species'] > 1, 'Uniqueness'] = 'CoV-common'
ep_model_data.reset_index(drop=True, inplace=True)
ep_model_data.rename(columns={'epitope': 'Epitope', 'species': 'No. Species', 'protein': 'Protein'}, inplace=True)
ep_model_data = ep_model_data[['Epitope', 'Uniqueness', 'No. Species', 'Protein', 'No. training TCRs',
                               'Balanced accuracy', 'AUC ROC', 'AUC PR']]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# epitope data only - Do I need it for the paper? seems that not
# ep_data = pd.merge(ep_data, prot_data2[['epitope', 'protein']], on='epitope', how='left')
# ep_data = ep_data[['epitope', 'species', 'protein', 'protein_info', 'SARS2protein', 'description', 'other']]

# ep_data.to_csv(f'{folder_out}/epitope_data.txt', sep='\t', index=False)
ep_model_data.to_csv(f'{folder_out}/Supplementary/TableS3_epitope_model_data.tsv', sep='\t', index=False)

# plot this data to remake PMe preprint plot
# ep_model_data = pd.read_csv(f'{folder_out}/Supplementary/TableS3_epitope_model_data.tsv', sep='\t')
ep_model_data.rename(columns={'Uniqueness': 'Specificity'}, inplace=True)
# ep_data_toplot = ep_model_data.groupby('Protein')['No. Species']
ax = sns.boxplot(x="Protein", y='No. Species', data=ep_model_data,
                 color=".95", saturation=0.1, showmeans=False)
# ax = sns.violinplot(x="Protein", y="No. Species", data=ep_model_data,
#                     inner='quartile', cut=0, color=".95", saturation=0.1)
ax = sns.swarmplot(x="Protein", y='No. Species', data=ep_model_data, hue='Specificity',
                   palette={'SC2-unique': 'tab:green', 'CoV-common': 'tab:purple'},
                   linewidth=0.5, size=3.5)  # color=".3", jitter=0.35,
ax.set(xlabel='Protein of origin', ylabel='Nidovirales species matches')
plt.tight_layout()
plt.savefig(f'{folder_out}/Fig1_NSpecies-vs-Protein.png')
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
