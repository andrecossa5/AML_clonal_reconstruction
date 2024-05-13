"""
Analyize mut stats.
"""

import os
import pandas as pd
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'var_selection')


##


# Read data
sample = 'sAML1'
path_results = os.path.join(path_results, sample)


jobs_df = pd.read_csv(os.path.join(path_results, f'{sample}_job.csv'), index_col=0)
jobs_df = jobs_df.reset_index().rename(columns={'index':'job'})
vars_df = pd.read_csv(os.path.join(path_results, f'{sample}_vars.csv'), index_col=0)
vars_df.reset_index(names='var').iloc[:,:2].drop_duplicates()['var'].value_counts().describe()
dataset_df = pd.read_csv(os.path.join(path_results, f'{sample}_dataset.csv'), index_col=0)
# 
# # Dataset level
# cols = ['density', 'genomes_redundancy', 'median_n_cells_per_var', 'median_n_vars_per_cell', 'n_vars', 'C>T', 'T>C' ]
# df_large = dataset_df.pivot_table(index='job', columns='metric', values='value')[cols]
# 
# df_large.describe()
# df_large['n_vars'].sort_values().reset_index().merge(jobs_df, on='job')


##



dataset_df.groupby('job_id').size().describe()
dataset_df['metric'].value_counts()


# Define transitions and transversions
patterns = [ 'A>C', 'T>G', 'A>T', 'A>G', 'G>A', 'C>G', 'C>A', 'T>A', 'G>C', 'G>T', 'N>T', 'C>T', 'T>C' ]
transitions = [pattern for pattern in patterns if pattern in ['A>G', 'G>A', 'C>T', 'T>C']]
transversions = [pattern for pattern in patterns if pattern not in transitions and 'N' not in pattern]

(
    dataset_df
    .loc[dataset_df['metric'].isin(patterns)]
    .pivot_table(values='value', columns='metric', index='job_id')
    .describe()
    .loc['50%']
    .sort_values(ascending=False)
)
(
    dataset_df
    .loc[dataset_df['metric'].isin(patterns)]
    .pivot_table(values='value', columns='metric', index='job_id')
    .assign(
        transitions=lambda x: x[transitions].sum(axis=1),
        transversions=lambda x: x[transversions].sum(axis=1),
        ratio_mut_signatures=lambda x: x['transitions']/x['transversions']
    )
    .sort_values('ratio_mut_signatures', ascending=False)
    .drop_duplicates()
    # .describe()
    # .loc['50%']
    # .sort_values(ascending=False)
)

jobs = vars_df.groupby('job_id')['prior'].median().sort_values().index[:10]
vars_ = vars_df.query('job_id in @jobs')
plt.plot(vars_['Variant_CellN'], vars_['prior'], 'o')

# plt.xlim([0,50])
# plt.ylim([0,50])
plt.show()

