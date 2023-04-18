import pandas as pd
from src.processPredictions.Functions import tcr_counter_sc_ds

folder = './data/scTCRs/E_MTAB_9357'

# clustering data
clust_df = pd.read_csv(f'{folder}/allCD8_clusters.tsv', sep='\t')
clust_sum = pd.read_csv(f'{folder}/CD8_clusters_summary.tsv', sep='\t')
# all TCRs per individual sample
seq_all = pd.read_csv(f'{folder}/CD8_Individual.tsv', sep='\t')
ep_data = pd.read_csv('./data/epitope_species.txt', sep='\t')
# one file with all TCRex predictions for SARS-CoV-2
tcrex_seq_all = pd.read_csv(f'{folder}/Individual.tcrex.cd8.csv', index_col=0)

tcrex_seq_all = tcrex_seq_all[tcrex_seq_all['score'] >= 0.9]
tcrex_seq_all['cloneCount'] = tcrex_seq_all['cloneCount'].fillna(0).astype('int')
# merge with epitope data from PMe
tcrex_seq_all = pd.merge(tcrex_seq_all, ep_data, on='epitope', how='left')
tcrex_seq_all['species'].fillna(0, inplace=True)
tcrex_seq_all['species'] = tcrex_seq_all['species'].astype(int)

# seq summary
seq_all.sort_values(by='ID', inplace=True)
seq_all_sum = seq_all.groupby('ID', as_index=False).agg({'cloneCount': 'sum',
                                                         'CDR3_beta': 'count'})
seq_all_sum.rename(columns={'cloneCount': 'total_count',
                            'CDR3_beta': 'unique_TCRs'}, inplace=True)

# add seq count data to clustering
clust_df.rename(columns={'junction_aa': 'CDR3_beta'}, inplace=True)
clust_df = clust_df.merge(seq_all, on='CDR3_beta', how='inner')
# 1 patient is not in clustering but is in seq_all file (which was clustered)
clust_df['cloneCount'] = clust_df['cloneCount'].fillna(0).astype('int')
# 1 patient INCOV145-AC has only 3 TCRs, won't add it to clusters
clust_df = clust_df.merge(seq_all_sum, on='ID', how='inner')

# annotate TCRs in clusters with TCRex
clusters_tcrex = clust_df.merge(tcrex_seq_all,
                                on=['CDR3_beta', 'TRBV_gene', 'TRBJ_gene',
                                    'ID', 'Disease', 'cloneCount'],
                                how='outer')

clusters_tcrex['%unique_TCRs'] = \
    100 * clusters_tcrex['unique_TCRs'] / clusters_tcrex['total_count']
clusters_tcrex['cloneFreq'] = \
    clusters_tcrex['cloneCount'] / clusters_tcrex['total_count']

clusters_tcrex.rename(columns={'ID': 'patient_id'}, inplace=True)
# clusters_tcrex.drop(columns=['score', 'bpr'], inplace=True)

for col in ['Disease', 'total_count', 'unique_TCRs', '%unique_TCRs']:
    clusters_tcrex.sort_values(by=['patient_id', col], inplace=True)
    clusters_tcrex[col] = clusters_tcrex.groupby(['patient_id'])[col].ffill()

# annotate the whole clusters with TCRex
for col in ['pathology', 'epitope', 'species']:
    clusters_tcrex.sort_values(by=['cluster', col], inplace=True)
    clusters_tcrex[col] = clusters_tcrex.groupby(['cluster'])[col].ffill()
clusters_tcrex.to_csv(f'{folder}/TCRex_allClusters.tsv', sep='\t', index=False)

clusters_tcrex = pd.read_csv(f'{folder}/TCRex_allClusters.tsv', sep='\t')
# keep only annotated clusters
clusters_tcrex.dropna(subset=['pathology', 'epitope'], inplace=True)
clusters_tcrex.drop('cluster', axis=1, inplace=True)
# clusters_tcrex[clusters_tcrex.isna().any(axis=1)]

# tmp df with TCR counts per covid/virus per TP per person
n = 0
tcr_per_virus = pd.DataFrame()
id_list = clusters_tcrex['patient_id'].unique().tolist()
id_list.sort()

for one_patient in id_list:
    n += 1
    counted = tcr_counter_sc_ds(clusters_tcrex, one_patient)
    # append to the df with tcr count for all individuals at all time points
    tcr_per_virus = tcr_per_virus.append(counted, ignore_index=True, sort=False)
    print(f"{n}/{len(id_list)}: "
          f"{one_patient} TCR specificities have been counted!")

tcr_per_virus.loc[tcr_per_virus['patient_id'].str.startswith('Healthy'),
                  'Disease'] = 'healthy'

tcr_per_virus.loc[tcr_per_virus['patient_id'].str.startswith('INCOV'),
                  'Disease'] = 'COVID-19'
row_filter = tcr_per_virus['patient_id'].str.startswith('INCOV')
tcr_per_virus.loc[row_filter, 'TP'] = tcr_per_virus.loc[
    row_filter, 'patient_id'].str.split('-', expand=True)[1]
tcr_per_virus.loc[tcr_per_virus['TP'] == 'BL', 'TP'] = 1
tcr_per_virus.loc[tcr_per_virus['TP'] == 'AC', 'TP'] = 2
tcr_per_virus['TP'].fillna(1, inplace=True)

tcr_per_virus.loc[row_filter, 'patient_id'] = tcr_per_virus.loc[
    row_filter, 'patient_id'].str.split('-', expand=True)[0]

# add response_redundancy
tcr_per_virus['response_redundancy'] = tcr_per_virus['N_SARS-CoV-2_TCRs'] \
                                       / tcr_per_virus['N_SARS-CoV-2_Eps']
tcr_per_virus['SC2-unique_response_redundancy'] = \
    tcr_per_virus['N_SARS-CoV-2_only_TCRs'] \
    / tcr_per_virus['N_SARS-CoV-2_only_Eps']
tcr_per_virus['CoV-common_response_redundancy'] = \
    tcr_per_virus['N_CoV_TCRs'] / tcr_per_virus['N_CoV_Eps']
# add repertoire breadth
tcr_per_virus['repertoire_breadth'] = tcr_per_virus['N_SARS-CoV-2_TCRs'] \
                                       / tcr_per_virus['unique_TCRs']
tcr_per_virus['SC2-unique_repertoire_breadth'] = \
    tcr_per_virus['N_SARS-CoV-2_only_TCRs']\
    / tcr_per_virus['unique_TCRs']
tcr_per_virus['CoV-common_repertoire_breadth'] = \
    tcr_per_virus['N_CoV_TCRs'] / tcr_per_virus['unique_TCRs']
tcr_per_virus.fillna(0, inplace=True)

tcr_per_virus = pd.read_csv(
    f'{folder}/TCRexPreds_withCounts_sc9357DS_clustered.tsv',
    index=False, sep='\t')
