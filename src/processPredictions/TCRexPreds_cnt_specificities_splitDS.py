import pandas as pd
from Functions import tcr_counter_split_ds

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            File processing: count X-specific TCRs per file, combine, save
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder = './data/splitDS_IMSEQ'
# output directory
directory_new = f'{folder}/TCRex_processed'

all_intersect2 = pd.read_csv((directory_new + '/viralTCRexPreds_withCounts_splitDS_intersected2x_meta.tsv'), sep='\t')
# 47 covid epitopes, 43 other viral epitopes

# PMe epitope data: common CV vs SARS-CoV-2 eps
ep_data = pd.read_csv('./data/epitope_species.txt', sep='\t')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TCR counts for pathology
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter out CDR3s with frequency less than 1 in 100.000 to compensate for different seq depth
all_intersect2 = all_intersect2[all_intersect2['freq'] > 1 / 100000]  # 3337 -> 3332
# 47 covid epitopes, 43 other viral epitopes

# filter out CDR3s with a score < 0.9
all_intersect2 = all_intersect2[all_intersect2['score'] >= 0.9]  # 3332 -> 817 ... ;(
# didn't explicitly filter, but all those 817 TCRs have bpr <= 1e-4, the same as in IR-Binder
# 25 covid epitopes (16 cr), 20 other viral epitopes

# merge with epitope data from PMe
all_intersect2 = pd.merge(all_intersect2, ep_data, on='epitope', how='left')
all_intersect2['species'].fillna(0, inplace=True)
all_intersect2['species'] = all_intersect2['species'].astype(int)

new_pat_id = []
for patient in all_intersect2['patient_id'].unique().tolist():
    for cell in ['cd4', 'cd8']:
        new_pat_id.append(f'{str(patient)}_{cell}')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTERSECTED CDR3s from repeats (occurs in at least 2 out of 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tcr_per_virus_inter = pd.DataFrame()

# tmp df with TCR counts per covid/virus per TP per person per cell type
n = 0
for one_patient in new_pat_id:
    n += 1
    counted = tcr_counter_split_ds(all_intersect2, one_patient)
    # append to the df with tcr count info for all individuals at all time points
    tcr_per_virus_inter = tcr_per_virus_inter.append(counted, ignore_index=True, sort=False)
    print(f"{n}/{len(new_pat_id)}: " + one_patient + " TCR counts have been added!")
    # log_file.write(f"{n}/{len(new_pat_id)}: " + one_patient + " TCR counts have been added!")
# log_file.close()

# tcr_per_virus_inter['fraq_SARS-CoV-2_TCRs'] = (tcr_per_virus_inter['SARS-CoV-2_TCRs']
#                                                / tcr_per_virus_inter['Total_TCRs'] * 100).astype(int)
tcr_per_virus_inter.sort_values(by=['disease_severity', 'day', 'new_id', 'T-cell_type'], inplace=True)
tcr_per_virus_inter.reset_index(inplace=True, drop=True)
# to solve rounding problems
# tcr_per_virus['fraq_SARS-CoV-2_TCRs'] = tcr_per_virus['fraq_SARS-CoV-2_TCRs'].apply(lambda x: round(x, 0))

# save one file with all TCR counts (with cutoff and freq filters; lists of viruses, eps and covid eps)
tcr_per_virus_inter.to_csv(f'{directory_new}/viralTCRexPreds_withCounts_splitDS_intersected2x_meta_perSpecificity.tsv',
                           index=False, sep='\t')

