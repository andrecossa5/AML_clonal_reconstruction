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
jobs_df = pd.read_csv(os.path.join(path_results, 'jobs.csv'), index_col=0)
jobs_df = jobs_df.reset_index().rename(columns={'index':'job'})
vars_df = pd.read_csv(os.path.join(path_results, 'vars_df_all.csv'), index_col=0)
dataset_df = pd.read_csv(os.path.join(path_results, 'dataset_df_all.csv'), index_col=0)

# Dataset level
cols = ['density', 'genomes_redundancy', 'median_n_cells_per_var', 'median_n_vars_per_cell', 'n_vars', 'C>T', 'T>C' ]
df_large = dataset_df.pivot_table(index='job', columns='metric', values='value')[cols]

df_large.describe()
df_large['n_vars'].sort_values().reset_index().merge(jobs_df, on='job')


##