import pandas as pd

folder = './data/scTCRs/E_MTAB_9357'

# all TCRs per individual sample
seq_all = pd.read_csv(f'{folder}/CD8_Individual.tsv', sep='\t')
ep_data = pd.read_csv('./data/epitope_species.txt', sep='\t')
# one file with all TCRex predictions for SARS-CoV-2
tcrex_seq_all = pd.read_csv(f'{folder}/Individual.tcrex.cd8.csv', index_col=0)  # 754, 10
tcrex_seq_all = tcrex_seq_all[tcrex_seq_all['score'] >= 0.9]  # 170, 10

seq_all.sort_values(by='ID', inplace=True)
seq_all_sum = seq_all.groupby('ID', as_index=False).agg({'cloneCount': 'sum',
                                                         'CDR3_beta': 'count'})
seq_all_sum.rename(columns={'cloneCount': 'total_count',
                            'CDR3_beta': 'unique_TCRs'}, inplace=True)

tcrex_seq_all = tcrex_seq_all.merge(seq_all_sum, on='ID', how='outer')
tcrex_seq_all['%unique_TCRs'] = \
    100 * tcrex_seq_all['unique_TCRs'] / tcrex_seq_all['total_count']
tcrex_seq_all['cloneFreq'] = \
    tcrex_seq_all['cloneCount'] / tcrex_seq_all['total_count']
tcrex_seq_all.rename(columns={'ID': 'patient_id'}, inplace=True)

# merge with epitope data from PMe
tcrex_seq_all = pd.merge(tcrex_seq_all, ep_data, on='epitope', how='left')
tcrex_seq_all['species'].fillna(0, inplace=True)
tcrex_seq_all['species'] = tcrex_seq_all['species'].astype(int)
tcrex_seq_all.to_csv(f'{folder}/TCRex_Cnts.tsv', sep='\t', index=False)
