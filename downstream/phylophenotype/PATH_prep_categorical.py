"""
PATH prep and tree viz.
"""

import os
import pickle
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.phylo import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


def prep_sample(sample):

    # Filter top
    d = (
        df.query('sample==@sample and time<=10 and time>=5')
        .groupby(['sample', 'filtering', 'bootstrap', 'solver', 'metric'])
        ['TBE'].median().sort_values(ascending=False)
        .reset_index()
        .iloc[0,:].to_dict()
    )
    print(d)

    # SET PARAMs
    filtering = d['filtering']
    solver = d['solver']
    metric = d['metric']             # default, cosine
    bootstrap = d['bootstrap']
    path_data = os.path.join(path_results, 'output', sample, filtering, 'input_folder')
    path_trees = os.path.join(
        path_results, 'output', sample, filtering, 
        solver, metric, bootstrap, 'trees.pickle'
    )

    # Read tree
    with open(path_trees, 'rb') as f:
        trees = pickle.load(f)

    # Get observed tree
    obs_tree = trees['observed']

    # Read cells, variants and build the afm
    variants = pd.read_csv(os.path.join(path_data, 'variants.csv'), index_col=0, header=None)
    m = pd.read_csv(os.path.join(path_data, 'meta.csv'), index_col=0)
    AD = load_npz(os.path.join(path_data, 'AD.npz')).astype(np.int16).A
    DP = load_npz(os.path.join(path_data, 'DP.npz')).astype(np.int16).A
    afm = AnnData(AD/DP, obs=m, var=variants)
    afm.obs = afm.obs.join(meta.query('sample_id==@sample')[["aggregated_ct", "malignant_class_occupancy"]])

    # Create annot dfs
    obs_tree.cell_meta = pd.DataFrame(
        afm.obs.loc[obs_tree.leaves].astype('str')
        ).join(
            pd.DataFrame(
            afm[obs_tree.leaves,:].X, index=obs_tree.leaves, columns=afm.var_names)
    )
    df_supports = (
        df.query('sample==@sample and filtering==@filtering and solver==@solver and bootstrap==@bootstrap and metric==@metric')
        .loc[obs_tree.internal_nodes, ['TBE', 'FBP', 'median_RF']]
    )

    return afm, obs_tree, df_supports, variants


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results/trees')
path_top_trees = os.path.join(path_results, 'top_trees')


##


# Read data
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
df = pd.read_csv(os.path.join(path_results, 'output/supports_df.csv'), index_col=0)


##


# Save individual sample top tree elements for PATH analysis
JOBS = {
    'CSF3R' : {
        'name_analysis' : 'CSF3R',
        'cov' : 'CSF3R',
        'samples' : ['sAML1'],
        'test' : lambda df,s: (df['sample_id']==s) & \
                              (df['CSF3R'].isin(['mutated', 'wild-type'])),
        'rebuild_tree' : True
    },
    'RUNX1' : {
        'name_analysis' : 'RUNX1',
        'cov' : 'RUNX1',
        'samples' : ['sAML1'],
        'test' : lambda df,s: (df['sample_id']==s) & \
                              (df['RUNX1'].isin(['mutated', 'wild-type'])),
        'rebuild_tree' : True
    },
    'FLT3' : {
        'name_analysis' : 'FLT3',
        'cov' : 'FLT3',
        'samples' : ['sAML1'],
        'test' : lambda df,s: (df['sample_id']==s) & \
                              (df['FLT3'].isin(['mutated', 'wild-type'])),
        'rebuild_tree' : True
    },
    'malignant_all_corr' : {
        'name_analysis' : 'malignant_all_corr',
        'cov' : 'malignant_class_occupancy',
        'samples' : ['sAML1', 'AML2', 'AML3', 'AML5'],
        'test' : lambda df,s: (df['sample_id']==s),               
        'rebuild_tree' : False
    },
    'malignant_aggregated_ct_t_rates' : {
        'name_analysis' : 'malignant_aggregated_ct_t_rates',
        'cov' : 'aggregated_ct',
        'samples' : ['sAML1'],
        'test' : lambda df,s: (df['sample_id']==s) & \
                              (df['malignant_class_occupancy']=='malignant') & \
                              (df['aggregated_ct'].isin(['Early_myeloid', 'Late_myeloid', 'HSC_MPP'])),
        'rebuild_tree' : True
    },
    'tme_aggregated_ct_t_rates' : {
        'name_analysis' : 'tme_aggregated_ct_t_rates',
        'cov' : 'aggregated_ct',
        'samples' : ['sAML1'],
        'test' : lambda df,s: (df['sample_id']==s) & \
                              (df['malignant_class_occupancy']=='tme'),
        'rebuild_tree' : True
    },
}

# Here we go
for key in JOBS:

    d = JOBS[key]
    for sample in d['samples']:

        # Prep
        make_folder(path_top_trees, sample, overwrite=False)
        path_sample = os.path.join(path_top_trees, sample)
        afm, obs_tree, df_supports, variants = prep_sample(sample)

        # Filter
        common = set(afm.obs_names) & set(meta.loc[d['test'](meta, sample)].index)
        to_remove = set(afm.obs_names)-common
        common = list(common)
        to_remove = list(to_remove)

        # Get tree 
        if d['rebuild_tree']:
            afm_filtered = nans_as_zeros(afm[common])
            obs_tree = build_tree(afm_filtered)
        else:
            obs_tree.remove_leaves_and_prune_lineages(to_remove)

        # Get covariate dfs
        assert all([ x in common for x in obs_tree.leaves ])
        afm = afm[obs_tree.leaves,:].copy()
        afm.obs[d['cov']] = meta.loc[obs_tree.leaves, d['cov']]
        counts = afm.obs[d['cov']].value_counts().to_frame('n_cells')
        counts['to_numeric'] = range(1,counts.shape[0]+1)
        cov_df = afm.obs[d['cov']].map(counts['to_numeric'].to_dict()).to_frame('cov')
        cov_df.to_csv(os.path.join(path_sample, f'cov_df_{d["name_analysis"]}.csv'))
        cov_df.loc[cov_df['cov'].isna()] = cov_df['cov'].value_counts().index[-1]+1
        counts.to_csv(os.path.join(path_sample, f'mapping_{d["name_analysis"]}.csv'))
        with open(os.path.join(path_sample, f'tree_{d["name_analysis"]}.newick'), 'w') as f:
            f.write(obs_tree.get_newick(record_branch_lengths=True))
            

##