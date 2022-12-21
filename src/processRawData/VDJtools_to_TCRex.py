import os
import pandas as pd

# IMSEQ data
folder = '/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/'
genes = ['A', 'B']

for gene in genes:
    path_in = folder + f'VDJtools_format_TR{gene}/'
    path_out = folder + f'TCRex_format_TR{gene}/'

    for file in os.listdir(path_in):
        data_in = pd.read_csv(path_in + file, sep='\t')
        data_in.rename(columns={'cdr3aa': 'CDR3_beta', 'v': 'TRBV_gene', 'j': 'TRBJ_gene'}, inplace=True)
        # ? remove TCRs with count of 1 (oe below 5? 10?)
        data_in = data_in[data_in['count'] >= 5]
        data_in = data_in[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene']]
        # remove non-canonical TCRs
        # check VJ-gene annotation
        data_in.to_csv(path_or_buf=(path_out + file), sep='\t', index=False)


# IR-Binder data
folder = '/Users/apost/Documents/CloudMail/PhD_2020/IMSEQ/IR-Binder_PRJEB38339/'
path_in = folder + 'MiXCR_parsed/'
if not os.path.exists(folder + 'TCRex_format'):
    os.makedirs(folder + 'TCRex_format')

# Patients, Donors
for label in ['Pt', 'HD']:
    # list of all IDs of patients (filenames match the start of TCRex prediction file names)
    id_list = [one_id for one_id in os.listdir(folder + 'MiXCR_parsed/')
               if one_id.startswith(label)]
    id_list.sort()
    # output location
    if not os.path.exists(folder + f'TCRex_format/{label}'):
        os.makedirs(folder + f'TCRex_format/{label}')
    path_out = folder + f'TCRex_format/{label}/'

    for file in id_list:
        data_in = pd.read_csv(path_in + file, sep='\t')
        # data_in = data_in[data_in['count'] >= 5]
        data_in = data_in[['TRBV_gene', 'CDR3_beta', 'TRBJ_gene']]
        # remove non-canonical TCRs
        # check VJ-gene annotation
        data_in.to_csv(path_or_buf=(path_out + file), sep='\t', index=False)