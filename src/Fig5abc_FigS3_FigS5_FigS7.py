import sys
import os
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import matplotlib.lines as mlines

folder = './data'
folder_out = './results'  # /Supplementary
folder_logs = './logs'

if not os.path.exists(folder_logs):
    os.makedirs(folder_logs)
    print(f'\nThe logs folder is created in the {folder_logs}')

all_data_scaled_week_max = pd.read_csv(
    f'{folder}/SplitMixed_mergedCD8_scaled_weekMax.tsv', sep='\t')
# ep_data = pd.read_csv('/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/Preprint_Metadata_PMe/epitope_data.txt', sep='\t')
# ep_data = pd.read_csv(f'{folder}/epitope_data.txt', sep='\t')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS7 - longitudinal dynamics per person week Max
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for parameter in ['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs']:
    mpl.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(6.6, 5))
    plt.figure(figsize=(9, 7))
    sns.lineplot(x='week', y=parameter, data=all_data_scaled_week_max,
                 hue_order=['Yes', 'No'], hue='critical_disease',
                 palette={'No': 'tab:blue', 'Yes': 'tab:orange'},
                 style='patient_id', markers=True, dashes=True)
    plt.ylabel(parameter)
    crit_patch = mlines.Line2D([], [], color='tab:orange', marker='o',
                               linestyle='None', markersize=5,
                               markeredgewidth=1.5,
                               label='Critical')
    noncrit_patch = mlines.Line2D([], [], color='tab:blue', marker='o',
                                  linestyle='None', markersize=5,
                                  markeredgewidth=1.5,
                                  label='Non-critical')
    plt.legend(handles=[crit_patch, noncrit_patch],
               title='COVID-19 severity', loc='upper right')
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=6)
    plt.tight_layout()
    plt.savefig(
        f'{folder_out}/Supplementary/FigS7_PersonDyn_{parameter}_weekMax.jpg',
        dpi=600)
    plt.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TCR to Epitope ratio - during the first 2 weeks is what matters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_data_scaled_week_max['CoV-common_TCR-Ep_ratio'] = \
    all_data_scaled_week_max['N_CoV_TCRs'] / all_data_scaled_week_max[
        'N_CoV_Eps']
all_data_scaled_week_max['SC2-unique_TCR-Ep_ratio'] = \
    all_data_scaled_week_max['N_SARS-CoV-2_only_TCRs'] / \
    all_data_scaled_week_max['N_SARS-CoV-2_only_Eps']
all_data_scaled_week_max[['CoV-common_TCR-Ep_ratio',
                          'SC2-unique_TCR-Ep_ratio']] = \
    all_data_scaled_week_max[['CoV-common_TCR-Ep_ratio',
                              'SC2-unique_TCR-Ep_ratio']].fillna(1)
all_data_scaled_week_max.loc[
    all_data_scaled_week_max['disease_status'] == 'recovered',
    'critical_disease'] = 'recovered'

# for supplementary figures
weeks12 = all_data_scaled_week_max[all_data_scaled_week_max['week'] < 3]
weeks38 = all_data_scaled_week_max[all_data_scaled_week_max['week'] >= 3]

# for statistics
crit = all_data_scaled_week_max[(all_data_scaled_week_max['week'] == 2) &
                                (all_data_scaled_week_max[
                                     'critical_disease'] == 'Yes')]
noncrit_active = all_data_scaled_week_max[
    (all_data_scaled_week_max['week'] == 2) &
    (all_data_scaled_week_max['critical_disease'] == 'No')]
noncrit = all_data_scaled_week_max[
    (all_data_scaled_week_max['critical_disease'] == 'recovered')]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS3 - dynamics of TCR repertoire and response parameters in individuals
# between weeks 1 and 2, critical and non-critical groups separately
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_week12_dynamics(data, varbls, path):
    import numpy as np
    first = varbls[0][0]
    second = varbls[0][1]
    data.loc[data['critical_disease'] == 'Yes',
             'critical_disease'] = 'Critical COVID-19'
    data.loc[data['critical_disease'] == 'No',
             'critical_disease'] = 'Non-critical COVID-19'
    for severity_ in data['critical_disease'].unique().tolist():
        print(severity_)
        data_long = data.melt(
            id_vars=['patient_id', 'week', 'critical_disease'],
            value_vars=[first, second],
            var_name='Specificity', value_name='Count')
        for parameter_ in varbls:
            data_plot = data_long[
                ((data_long['week'] == 1) | (data_long['week'] == 2))
                & (data_long['critical_disease'] == severity_)].copy()
            data_plot['week'] = data_plot['week'].astype(str)
            data_plot.sort_values(by=['week'], inplace=True)
            # Plot dots with trend line
            mpl.rcParams['font.family'] = 'Arial'
            fig = plt.figure(figsize=(6.6, 5))
            ax = sns.relplot(x='week', y='Count', data=data_plot,
                             col='critical_disease', kind='line',
                             markers=True,
                             # {'Frac_CoV_TCRs': "x", 'Frac_SARS-CoV-2_only_TCRs': "o"},
                             dashes=[(1, 1), (1, 1)], style='Specificity',
                             estimator=np.median, ci=95, n_boot=1000,
                             seed=None, sort=True, err_style='band',
                             hue='Specificity',
                             palette={parameter_[0]: 'tab:purple',
                                      parameter_[1]: 'tab:green'}, legend=False)
            # Plot individual lines
            ax = sns.lineplot(x='week', y='Count', data=data_plot,
                              markers=False,
                              # {'Frac_CoV_TCRs': "x", 'Frac_SARS-CoV-2_only_TCRs': "o"},
                              dashes=False, style='patient_id',
                              hue='Specificity',
                              palette={parameter_[0]: 'tab:purple',
                                       parameter_[1]: 'tab:green'},
                              legend=False)
            ax = sns.stripplot(x="week", y='Count', data=data_plot,
                               hue='Specificity',
                               palette={parameter_[0]: 'tab:purple',
                                        parameter_[1]: 'tab:green'},
                               dodge=True, size=5, linewidth=1)
            ax = sns.violinplot(x="week", y='Count', data=data_plot,
                                hue='Specificity',
                                palette={parameter_[0]: 'tab:purple',
                                         parameter_[1]: 'tab:green'},
                                split=True,
                                scale='count', saturation=0.4,
                                inner='quartile', dodge=True, cut=0)  # bw=.4,
            ax.legend_.remove()
            ax.set(xlabel='week', ylabel=f'{parameter_[2]}')
            cov_patch = mlines.Line2D([], [], color='tab:purple', marker='o',
                                      linestyle='None', markersize=5,
                                      markeredgewidth=1.5,
                                      label='CoV-common')
            sc2_patch = mlines.Line2D([], [], color='tab:green', marker='o',
                                      linestyle='None', markersize=5,
                                      markeredgewidth=1.5,
                                      label='SC2-unique')
            plt.legend(handles=[cov_patch, sc2_patch], title='Specificity')
            plt.tight_layout()
            plt.title(severity_)
            plt.savefig(
                f'{path}/FigS3_{parameter_[2]}_w1-2_{severity_.split(" ")[0]}_CoV_vs_SC2.jpg',
                bbox_inches='tight', dpi=600)
            plt.close()
            print(
                f'FigS3_{parameter_[2]}_w1-2_{severity_}_CoV_vs_SC2.png is saved in {path}')


l = [[['Frac_CoV_TCRs', 'Frac_SARS-CoV-2_only_TCRs', 'TCR Fractions']],
     [['%CoV-common_TCRs', '%SC2-unique_TCRs', '%_TCRs']],
     [['N_CoV_Eps', 'N_SARS-CoV-2_only_Eps', 'Epitope counts']],
     [['CoV-common_TCR-Ep_ratio', 'SC2-unique_TCR-Ep_ratio',
       'Response redundancy']]]
for i in l:
    plot_week12_dynamics(all_data_scaled_week_max[
                             all_data_scaled_week_max['critical_disease']
                             != 'recovered'], i,
                         f'{folder_out}/Supplementary')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# N recognized epitopes vs N specific TCRs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to have outputs in a log file
orig_stdout = sys.stdout
log_f = open(f'{folder_logs}/Fig5abc.log', 'w+')
sys.stdout = log_f

print('\nCRITICAL ACTIVE patients week 2')
print(f"Median TCR/SC2-unique_Epitope ratio: "
      f"{crit['SC2-unique_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{crit['SC2-unique_TCR-Ep_ratio'].describe()[3]}-"
      f"{crit['SC2-unique_TCR-Ep_ratio'].describe()[7]}]")
print(f"Median TCR/CoV-common_Epitope ratio: "
      f"{crit['CoV-common_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{crit['CoV-common_TCR-Ep_ratio'].describe()[3]}-"
      f"{crit['CoV-common_TCR-Ep_ratio'].describe()[7]}]")

print('\nNON-CRITICAL ACTIVE patients week 2')
print(f"Median TCR/SC2-unique_Epitope ratio: "
      f"{noncrit_active['SC2-unique_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit_active['SC2-unique_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit_active['SC2-unique_TCR-Ep_ratio'].describe()[7]}]")
print(f"Median TCR/CoV-common_Epitope ratio: "
      f"{noncrit_active['CoV-common_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit_active['CoV-common_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit_active['CoV-common_TCR-Ep_ratio'].describe()[7]}]")

print('\nNON-CRITICAL RECOVERED patients')
print(f"Median TCR/SC2-unique_Epitope ratio: "
      f"{noncrit['SC2-unique_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit['SC2-unique_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit['SC2-unique_TCR-Ep_ratio'].describe()[7]}]")
print(f"Median TCR/CoV-common_Epitope ratio: "
      f"{noncrit['CoV-common_TCR-Ep_ratio'].describe()[5]}; "
      f"range=[{noncrit['CoV-common_TCR-Ep_ratio'].describe()[3]}-"
      f"{noncrit['CoV-common_TCR-Ep_ratio'].describe()[7]}]")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MU N of all epitopes in crit vs not active, week 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha = 0.05
print('\nN of epitopes recognized by critical and non-critical patients'
      ' during week 2')
for parameter in ['N_SARS-CoV-2_Eps', 'N_SARS-CoV-2_only_Eps', 'N_CoV_Eps']:
    stat, p = mannwhitneyu(crit[parameter], noncrit_active[parameter])
    # interpret
    if p > alpha:
        print(f'{parameter}: Same distribution (fail to reject H0)')
    else:
        print(f'{parameter}: Different distribution (reject H0)')
    print('    Statistics=%.3f, p=%.5f' % (stat, p))
    crit_tcrs = crit[parameter].to_list()
    non_crit_active_tcrs = noncrit_active[parameter].to_list()
    n_pairs = 0
    for i in non_crit_active_tcrs:
        for j in crit_tcrs:
            if i > j:
                n_pairs += 1
    effect = n_pairs / (len(crit_tcrs) * len(non_crit_active_tcrs))
    print(f'    Effect size for {parameter} is {effect}')

print('\nN of epitopes recognized by active non-critical patients'
      ' during week 2 and recovered non-critical patients')
for parameter in ['N_SARS-CoV-2_Eps', 'N_SARS-CoV-2_only_Eps', 'N_CoV_Eps']:
    stat, p = mannwhitneyu(noncrit_active[parameter], noncrit[parameter])
    # interpret
    if p > alpha:
        print(f'{parameter}: Same distribution (fail to reject H0)')
    else:
        print(f'{parameter}: Different distribution (reject H0)')
    print('    Statistics=%.3f, p=%.5f' % (stat, p))
    non_crit_active_tcrs = noncrit_active[parameter].to_list()
    non_crit_tcrs = noncrit[parameter].to_list()
    n_pairs = 0
    for i in non_crit_tcrs:
        for j in crit_tcrs:
            if i > j:
                n_pairs += 1
    effect = n_pairs / (len(noncrit_active) * len(noncrit))
    print(f'    Effect size for {parameter} is {effect}')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spearman & Plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_ep_redundancy_spearman(data, weeks):
    crit = data[(data['week'] >= int(weeks[0])) &
                (data['week'] <= int(weeks[1])) &
                (data['critical_disease'] == 'Yes')]
    noncrit_active = data[(data['week'] >= int(weeks[0])) &
                          (data['week'] <= int(weeks[1])) &
                          (data['critical_disease'] == 'No')]
    noncrit = data[data['critical_disease'] == 'recovered']
    data_to_plot = pd.concat([crit, noncrit_active, noncrit])
    # what time period will be plotted and for main or Supp figure
    if weeks[0] == weeks[1]:
        period = f'w{weeks[1]}'
        figure = 'Fig5'
    else:
        period = f'w{weeks[0]}-{weeks[1]}'
        figure = 'Supplementary/FigS5'

    print(f'\nN recognized epitopes vs N of specific TCRs, {period}')
    for parameter in [['N_CoV_TCRs', 'N_CoV_Eps', 'CoV-common'],
                      ['N_SARS-CoV-2_only_TCRs', 'N_SARS-CoV-2_only_Eps',
                       'SC2-unique']]:
        print(f'{parameter[2]}')
        crit_corr = spearmanr(crit[parameter[0]], crit[parameter[1]])
        print(f'Critical: {crit_corr}')
        noncrit_active_corr = spearmanr(noncrit_active[parameter[0]],
                                        noncrit_active[parameter[1]])
        print(f'Non-critical active: {noncrit_active_corr}')
        noncrit_corr = spearmanr(noncrit[parameter[0]], noncrit[parameter[1]])
        print(f'Non-critical recovered: {noncrit_corr}')

        # to plot all lines on one plot, looks messy
        # noncrit_recovered = all_data_scaled_week_max[
        #     (all_data_scaled_week_max['disease_status'] == 'recovered')
        #     & (all_data_scaled_week_max['critical_disease'] == 'No')]
        # noncrit_recov_corr = spearmanr(noncrit[parameter[0]],
        #                                noncrit[parameter[1]])
        # print(f'Non-critical recovered: {noncrit_recov_corr}')

        # simple band legend without markers
        # crit_label = mpatches.Patch(color='tab:orange', hatch='o',
        #                             label=f'rho={crit_corr[0]:.2f},'
        #                                   f'p={crit_corr[1]:.1e}')
        # noncrit_label = mpatches.Patch(color='lightskyblue',
        #                                label=f'rho={noncrit_corr[0]:.2f},'
        #                                      f'p={noncrit_corr[1]:.1e}')
        # noncrit_active_label = mpatches.Patch(color='tab:blue',
        #                                       label=f'rho={noncrit_active_corr[0]:.2f},'
        #                                             f'p={noncrit_active_corr[1]:.1e}')
        # ref_label = mpatches.Patch(color='black', label='TCR/Epitope ratio = 1')

        # legend with markers instead of bands
        crit_label = mlines.Line2D([], [], color='tab:orange',
                                   marker='v', markersize=6,
                                   label=f'Critical active {period}: rho={crit_corr[0]:.2f},'
                                         f'p={crit_corr[1]:.1e}')
        noncrit_active_label = mlines.Line2D([], [], color='tab:blue',
                                             marker='x', markersize=6,
                                             label=f'Non-critical active {period}: '
                                                   f'rho={noncrit_active_corr[0]:.2f},'
                                                   f'p={noncrit_active_corr[1]:.1e}')
        noncrit_label = mlines.Line2D([], [], color='lightskyblue',
                                      marker='o', markersize=6,
                                      label=f'Non-citical recovered: '
                                            f'rho={noncrit_corr[0]:.2f},'
                                            f'p={noncrit_corr[1]:.1e}')
        ref_label = mlines.Line2D([], [], color='black',
                                  marker='.', markersize=6,
                                  label='TCR/Epitope ratio = 1')

        mpl.rcParams['font.family'] = 'Arial'
        fig = plt.figure(figsize=(6.6, 5))
        ax = sns.lmplot(x=parameter[0], y=parameter[1],
                        data=data_to_plot,
                        scatter_kws={'s': 20}, x_ci='sd',
                        markers=['v', 'x', 'o'],
                        hue_order=['Yes', 'No', 'recovered'],
                        hue='critical_disease',
                        palette={'Yes': 'tab:orange', 'No': 'tab:blue',
                                 'recovered': 'lightskyblue'}, ci=95,
                        legend=False)
        # g = sns.lmplot(x=parameter[0], y=parameter[1],
        #                data=data_toplot,
        #                scatter_kws={'s': 20}, x_ci=None, ci=None,
        #                markers=['o', 'o', 'x', 'x', 'v'],
        #                hue_order=['critical_active_w12', 'non_critical_active_w12',
        #                           'critical_active_w3+', 'non_critical_active_w3+',
        #                           'recovered'],
        #                hue='critical_disease',
        #                palette={'non_critical_active_w12': 'tab:blue',
        #                         'critical_active_w12': 'tab:orange',
        #                         'critical_active_w3+': '#a6611a',
        #                         'non_critical_active_w3+': '#253494',
        #                         'recovered': 'lightblue'}, legend=False)
        ax.set(ylim=(0, data[parameter[1]].max() + 1))
        # plt.axis('square')
        # g.spines['bottom'].set_position('zero')
        plt.axline([0, 0], [1, 1], color='black', linestyle=':')
        plt.xlabel(f'No. {parameter[2]} TCRs')
        plt.ylabel(f'No. {parameter[2]} Eps')
        plt.legend(
            handles=[crit_label, noncrit_active_label, noncrit_label,
                     ref_label],
            loc='upper center', bbox_to_anchor=(0.5, -0.15),
            fancybox=True, shadow=True)
        # plt.legend(title='critical_disease')
        plt.tight_layout()
        # plt.show()
        plt.savefig(
            f'{folder_out}/{figure}_RegressDots_TCRs-vs-Eps_{parameter[2]}'
            f'_active_{period}_recovered_Ref.jpg', dpi=600)
        plt.close()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 5a,b (week 2 for active patients vs recovered)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ep_redundancy_spearman(all_data_scaled_week_max, '22')
# FigS5
plot_ep_redundancy_spearman(all_data_scaled_week_max, '12')
plot_ep_redundancy_spearman(all_data_scaled_week_max, '38')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EPITOPE INTERSECTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
weeks12 = all_data_scaled_week_max[all_data_scaled_week_max['week'] == 2]
dis_sev = weeks12['critical_disease'].unique().tolist()
# dict where keys - groups, values - all epitopes per group
eps_per_group = {}
# dict where key1 - group, value1 - dict where key2 - person in group,
# value2 - all epitopes per this person
eps_per_ind = {}

for severity in dis_sev:
    df1 = weeks12[weeks12['critical_disease'] == severity]
    # create a set of epitopes per group
    eps_gr = set(a.strip("''") for b in
                 df1['all_SARS-CoV-2_epitopes'].str.strip('[]').str.split(', ')
                 for a in b)
    eps_per_group[severity] = eps_gr
    eps_per_ind[severity] = {}
    for patient in df1['patient_id'].unique().tolist():
        df2 = df1.loc[df1['patient_id'] == patient, ['patient_id',
                                                     'all_SARS-CoV-2_epitopes']]
        # create a set of epitopes per person
        eps_i = set(a.strip("''") for b in
                    df2['all_SARS-CoV-2_epitopes'].str.strip('[]').str.split(
                        ', ') for a in b)
        eps_per_ind[severity][str(patient)] = eps_i
# removing empty string from epitope list
for v in eps_per_group.values():
    # print(v)
    try:
        v.remove('')
    except KeyError:
        pass
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# intersections BETWEEN groups
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a = eps_per_group['No']
b = eps_per_group['Yes']
# non-critical only
a_only = a - (a & b)
# 'ALSKGVHFV',
#  'ITEEVGHTDLMAAY',
#  'KAYNVTQAF',
#  'KLPDDFTGCV',
#  'QECVRGTTVL',
#  'RQLLFVVEV',
#  'SSNVANYQK',
#  'YFPLQSYGF',
#  'YLDAYNMMI'
# critical only
b_only = b - (a & b)
print(
    '\nEpitopes recognized only by critical patients during week 2:')
print(b_only)

print(
    '\nEpitopes recognized only by non-critical patients during week 2:')
print(a_only)

# details about 9 epitopes which occur only in non-critical patients
noncrit_only_ep = pd.DataFrame({'epitope': list(a_only)})
# noncrit_only_ep = pd.merge(noncrit_only_ep, ep_data, on='epitope', how='left')
# noncrit_only_ep.sort_values(by=['species', 'protein'], inplace=True)
# 3 unique vs 6 common: 3/6=0.5 vs 19unique/28common=0.68
# how many? 19 unique vs 28 common eps
# print(len(ep_data.loc[ep_data['species'] == 1, 'epitope'].to_list()))

# how many patients do recognize those 9 epitopes?
print('Number of non-critical patients recognizing those unique epitopes:')
epi_n_noncrit = {}
for epi in noncrit_only_ep['epitope'].to_list():
    cnt = [item for el in
           [v2 for v1 in eps_per_ind.values() for v2 in v1.values()] for item
           in el].count(epi)
    epi_n_noncrit[epi] = cnt
print(epi_n_noncrit)  # only 1 or 2 patients
# 2 circles
# eps_per_group['Critical active week 2'] = eps_per_group.pop('Yes')
# eps_per_group['Non-critical active week 2'] = eps_per_group.pop('No')
# venn2(subsets=[a, b], set_labels=('Non-critical patients', 'Critical patients'))
# venn(eps_per_group, cmap=['tab:blue', 'tab:orange'])
# plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 5c
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6.6, 5))
v = venn2(subsets=(len(a_only), len(b_only), len(a & b)),
          set_labels=('Non-critical\nactive    \nweek 2   ',
                      'Critical\n active\nweek 2'),
          set_colors=("#1f77b4",
                      "#ff7f0e"))
# set colors of circles
v.get_patch_by_id('11').set_color('#ff7f0e')
v.get_patch_by_id('11').set_edgecolor('tab:orange')
v.get_patch_by_id('10').set_color('#1f77b4')
v.get_patch_by_id('10').set_edgecolor('tab:blue')
v.get_patch_by_id('11').set_alpha(0.45)
# reposition labels
lbl = v.get_label_by_id('A')
x, y = lbl.get_position()
lbl.set_position((x - 0.2, y + 0.55))
lbl = v.get_label_by_id('B')
x, y = lbl.get_position()
lbl.set_position((x - 0.05, y + 0.35))
green_patch = mlines.Line2D([], [], color='tab:green', marker='|',
                            linestyle='None', markersize=60,
                            markeredgewidth=1.5,
                            label='5 SC2-unique: \nALSKGVHFV, \nEILDITPCSF, '
                                  '\nFADDLNQLTGY, \nNLDSKVGGNY, \nYFPLQSYGF')
purple_patch = mlines.Line2D([], [], color='tab:purple', marker='|',
                             linestyle='None', markersize=50,
                             markeredgewidth=1.5,
                             label='4 CoV-common: \nFLNGSCGSV, \nFLPRVFSAV, '
                                   '\nHLVDFQVTI, \nKLSYGIATV')

fig.legend(handles=[green_patch, purple_patch], title="SARS-CoV-2 epitopes",
           loc='center left')
plt.tight_layout()
plt.savefig(f'{folder_out}/Fig5c_Eps_intersect_crit_not_w2.jpg', dpi=600)
plt.close()
# a_subset = {key: eps_per_group[key] for key in dis_sev}
# venn(a_subset)
sys.stdout = orig_stdout
log_f.close()

print(f'Log file Fig5abc.log is created in the {folder_logs} directory.')
print(f'Figures 5 a-c are saved in the {folder_out} directory.')
print(
    f'Figures S5 a-d, S7 a-b are saved in the {folder_out}/Supplementary directory.')

# 24 'orf1ab polyprotein' vs 23 others; 4 out of 9 are ORF1ab
# print(len(ep_data.loc[ep_data['protein'] == 'ORF1ab', 'epitope'].to_list()))

# patients and time when those 9 epitopes are recognized
# for epi in noncrit_only_ep['epitope'].to_list():
#     print(f'\nEpitope: {epi}')
#     for a, b, c in zip(all_data_scaled_week_max['patient_id'],
#                        all_data_scaled_week_max['all_SARS-CoV-2_epitopes'],
#                        all_data_scaled_week_max['week']):
#         if epi in b:
#             print(f'Patient {a}, week {c}')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# intersections WITHIN groups
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Epitopes in all patients within a group
# common_noncrit = set.intersection(
#     *(set(val) for val in eps_per_ind['No'].values()))  # 0
# common_crit = set.intersection(
#     *(set(val) for val in eps_per_ind['Yes'].values()))  # 0

# modify epitope lists
# all_data_scaled_week_max = pd.read_csv(
#     f'{folder}/imseq_binder_merged_withEpitopes.tsv', sep='\t')
# all_data_scaled_week_max['all_SARS-CoV-2_epitopes'] = all_data_scaled_week_max[
#     'all_SARS-CoV-2_epitopes'] \
#     .apply(lambda x: set(x.strip('[]').split(', ')))

# are there TCRs to all available epitope models?
# all_predicted_eps = set()
# for i in all_data_scaled_week_max['all_SARS-CoV-2_epitopes']:
#     all_predicted_eps = all_predicted_eps.union(i)
#     print(len(all_predicted_eps))

# to plot 5a,b,c
# import matplotlib.gridspec as gridspec
# mpl.rcParams['font.family'] = 'Arial'
# # fig = plt.figure(figsize=(6.6, 5))
# gs = gridspec.GridSpec(2, 4)
# gs.update(wspace=0.5)
# ax1 = plt.subplot(gs[0, :2], )
# ax2 = plt.subplot(gs[0, 2:])
# ax3 = plt.subplot(gs[1, 1:3])
# plt.tight_layout()
# plt.show()
