#!/usr/bin/python

import sys
import os
from mito_utils.preprocessing import *
import warnings
warnings.simplefilter('ignore')


##


# Paths
path_data = sys.argv[1]
path_meta = sys.argv[2]
sample = sys.argv[3]


##


def main():

    # Read, format and retain good cells
    afm = read_one_sample(path_data, sample, with_GBC=False)
    meta = pd.read_csv(path_meta, index_col=0)
    meta =  meta.query('sample_id==@sample')
    meta.index = meta.index.map(lambda x: x.split('-')[0])
    afm.obs = afm.obs.join(meta)
    afm = afm[~afm.obs['malignant_class_occupancy'].isna(),:].copy()
    print(afm)

    # Calculate vars_df, as in in Weng et al., 2024, and Miller et al. 2022 before.
    mean_sites = afm.uns['per_position_coverage'].mean(axis=0)
    vars_df = pd.DataFrame(
        np.column_stack([
            afm.X.mean(axis=0),
            afm.var_names.map(lambda x: mean_sites[x.split('_')[0]]),
            np.ma.getdata(np.ma.mean(np.ma.masked_less_equal(afm.layers['quality'], 0), axis=0))
        ]),
        columns=['mean_af', 'mean_cov', 'quality'],
        index=afm.var_names
    )

    # Calculate the number of cells that exceed VAF thresholds 0, 1, 5, 10, 50 as in Weng et al., 2024
    vars_df = (
        vars_df
        .assign(
            n0 = lambda x: np.sum(afm.X==0, axis=0),
            n1 = lambda x: np.sum(afm.X>.01, axis=0),
            n5 = lambda x: np.sum(afm.X>.05, axis=0),
            n10 = lambda x: np.sum(afm.X>.1, axis=0),
            n50 = lambda x: np.sum(afm.X>.5, axis=0),
            Variant_CellN = lambda x: np.sum(afm.X==0, axis=0),
        )
    )

    # Filter Wang et al., 2024
    min_site_cov = 5
    min_var_quality = 30
    min_frac_negative = .5
    min_n_positive = 2
    vars_df = (
        vars_df.loc[
            (vars_df['mean_cov']>min_site_cov) & \
            (vars_df['quality']>min_var_quality) & \
            (vars_df['n0']>min_frac_negative) & \
            (vars_df['mean_cov']>min_site_cov) & \
            (vars_df['Variant_CellN']>min_n_positive) 
        ]
    )

    # Write 
    vars_df.to_csv(os.path.join(path_data, f'{sample}_vars.csv'))


    ##


# Run
if __name__ == '__main__':
    main()


##