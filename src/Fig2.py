import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
from scipy.stats import mannwhitneyu

folder = './data'
folder_out = './results'

imseq_file = 'MetaForPlotting_viralTCRexPreds_withCounts_splitDS_intersected2x_meta_perSpecificity.tsv'
tcr_per_virus_inter = pd.read_csv(f'{folder}/splitDS_IMSEQ/TCRex_processed/{imseq_file}', sep='\t')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# M-U STATISTICS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of predicted SC2 TCRs to demonstrate specificity of TCRex models
alpha = 0.005
parameter = 'N_SARS-CoV-2_TCRs'
data1 = tcr_per_virus_inter.loc[tcr_per_virus_inter['T-cell_type'] == 'cd4', parameter]
data2 = tcr_per_virus_inter.loc[tcr_per_virus_inter['T-cell_type'] == 'cd8', parameter]
# compare samples
stat, p = mannwhitneyu(data1, data2)
print('\nStatistics=%.3f, p=%.3f' % (stat, p))
# Statistics=212.000, p=0.002 = Different distribution (reject H0)
# Same distribution (fail to reject H0 when a=0.001)
# interpret
if p > alpha:
    print(f'Same distribution of {parameter} between CD4 and CD8 TCRs (fail to reject H0)')
else:
    print(f'Different distribution of {parameter} between CD4 and CD8 TCRs (reject H0)')
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Figure 2:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # BOXplots with STAT significance
    # all disease groups together
    mpl.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(6.6, 5))
    plt.tight_layout()
    plt.savefig(f'{folder_out}/Fig1_NSpecies-vs-Protein.jpg', dpi=600)
    ax = sns.boxplot(x="T-cell_type", y=parameter, data=tcr_per_virus_inter,
                     palette='Set1', showmeans=True, order=['cd8', 'cd4'],
                     meanprops={'marker': '*', 'markerfacecolor': 'white',
                                'markeredgecolor': 'black'})
    ax = sns.stripplot(x="T-cell_type", y=parameter, data=tcr_per_virus_inter,
                       # tcr_per_virus_inter,
                       order=['cd8', 'cd4'], palette='Set1', dodge=True,
                       linewidth=1)
    ax.set(xlabel='T-cell population', ylabel='No. predicted TCRs')
    plt.tight_layout()
    # annotate median
    # medians are 14 vs 4 for cd8 vs cd4
    medians = tcr_per_virus_inter.groupby(['T-cell_type'], sort=False)[parameter].median()
    medians.sort_index(ascending=False, inplace=True)
    vertical_offset = tcr_per_virus_inter[
                          parameter].median() * 0.1  # offset from median for display
    for xtick in ax.get_xticks():
        ax.text(xtick, medians[xtick] + vertical_offset, medians[xtick],
                horizontalalignment='center', size='small', color='w',
                weight='semibold')
    # statistical annotation (M-U in stats section higher up)
    x1, x2 = 0, 1
    y, h, col = tcr_per_virus_inter[parameter].max() + 4, 2, 'k'
    plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
    plt.text((x1 + x2) * .5, y + h, "p < 0.01", ha='center', va='bottom',
             color=col)
    ax.tick_params(labelsize=10)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    plt.savefig(f'{folder_out}/Fig2_box_signif_{parameter}_cd4-cd8_noasympt.jpg', dpi=600)
    plt.close()
