"""
Prep data for unsupervised gene modules discovery.
"""

import os
import pickle
from scipy.sparse import load_npz
import scanpy as sc
from anndata import AnnData
from mito_utils.phylo import *
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
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
path_results = os.path.join(path_main, 'results', 'trees')
path_top_trees = os.path.join(path_results, 'top_trees')


##


# Read data
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
df = pd.read_csv(os.path.join(path_results, 'output/supports_df.csv'), index_col=0)
adata = sc.read(os.path.join(path_data, 'expression', 'full_adata.h5ad'))
HVGs = pd.read_csv(os.path.join(path_data, 'expression', 'HVG_final.csv'), index_col=0).iloc[:,0]

# Prep data for unsupervised PATH
key = 'malignant_unsupervised'
test = lambda df,s: (df['sample_id']==s) & \
                    (df['malignant_class_occupancy']=='malignant') & \
                    (df['aggregated_ct'].isin(['Early_myeloid', 'Late_myeloid', 'HSC_MPP']))
# Prep
samples = ['sAML1', 'AML2', 'AML3', 'AML5']
for sample in samples:

    # Prep folder    
    make_folder(path_top_trees, sample, overwrite=False)
    path_sample = os.path.join(path_top_trees, sample)
    afm, obs_tree, df_supports, variants = prep_sample(sample)

    # Filter
    common = set(afm.obs_names) & set(meta.loc[test(meta, sample)].index)
    to_remove = set(afm.obs_names)-common
    common = list(common)
    to_remove = list(to_remove)

    # Get tree 
    afm_filtered = nans_as_zeros(afm[common])
    obs_tree = build_tree(afm_filtered)

    # Save data
    assert all([ x in common for x in obs_tree.leaves ])
    M_hvgs = pd.DataFrame(adata[obs_tree.leaves, HVGs].X, index=obs_tree.leaves, columns=HVGs)
    M_hvgs.to_csv(os.path.join(path_sample, f'HVGs_{key}.csv'))
    with open(os.path.join(path_sample, f'tree_{key}.newick'), 'w') as f:
        f.write(obs_tree.get_newick(record_branch_lengths=True))


##