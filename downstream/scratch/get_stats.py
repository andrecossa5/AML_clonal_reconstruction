#!/usr/bin/python

import os
import sys
from itertools import product
from mito_utils.preprocessing import *
import warnings
warnings.simplefilter('ignore')


##


# Args
# path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
# sample = 'sAML1'
path_main = sys.argv[1]
sample = sys.argv[2]

# Paths
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')


##


# Params
filtering_kwargs = {
    'min_site_cov' : [5, 10, 35], 
    'min_var_quality' : [30], 
    'min_frac_negative' : [.5, .75, .9],
    'min_n_positive' : [2, 5, 10],
    'low_confidence_af' : [.001, .01, .1], 
    'high_confidence_af' : [.01, .1, .5], 
    'min_prevalence_low_confidence_af' : [.001, .01, .1], 
    'min_cells_high_confidence_af' : [2, 5, 10]
}


##


def main():
    
    # Read AFM and meta
    afm = read_one_sample(path_data, sample)
    meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
    meta = meta.query('sample_id==@sample')
    meta.index = meta.index.map(lambda x: x.split('-')[0])
    afm.obs = afm.obs.join(meta[['malignant_class_occupancy']])
    jobs = list(product(*filtering_kwargs.values()))
    jobs = [ dict(zip(filtering_kwargs.keys(), j)) for j in jobs ]
    # Save jobs
    pd.DataFrame(jobs).to_csv(os.path.join(path_results, 'jobs.csv'))

    # Run filtering jobs
    n = len(jobs)
    for i, job in enumerate(jobs):
        if d['low_confidence_af'] > d['high_confidence_af']:
            print(f'Jobs: {i}/{n}')
            vars_df, dataset_df, a = filter_cells_and_vars(
                afm, filtering='weng2024', filtering_kwargs=job,
                lineage_column='malignant_class_occupancy',
                spatial_metrics=False, fit_mixtures=False, 
                path_priors=path_priors
            )
        # Save
        vars_df.to_csv(os.path.join(path_results, 'vars_df', f'job_{i}.csv'))
        dataset_df.to_csv(os.path.join(path_results, 'dataset_df', f'job_{i}.csv'))
    


##


# Run
if __name__ == '__main__': 
    main()





