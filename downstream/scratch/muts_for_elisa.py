"""
Muts for Elisa
"""

import os
from mito_utils.preprocessing import *


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')
path_functional = os.path.join(path_main, 'results', 'functional')


##


# Read MT_data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = meta.query('sample_id==@sample')
afm.obs = afm.obs.join(meta)


##


# Filter leukemic cells and MT variants
_, a = filter_cells_and_vars(
    afm, sample_name=sample, 
    path_priors=path_priors, 
    max_prior=.1, 
    filtering='MI_TO',
    filtering_kwargs={
        "min_n_confidently_detected" : 2,
        "af_confident_detection" : 0.05,
        "min_frac_negative" : 0
    },
    af_confident_detection=.05,
    fit_mixtures=False,
    only_positive_deltaBIC=False,
    max_AD_counts=3,
    lineage_column='malignant_class_occupancy'
)

# Save
_.to_csv(os.path.join(path_functional, 'dataset.csv'))
a.var.to_csv(os.path.join(path_functional, 'variants.csv'))


##

