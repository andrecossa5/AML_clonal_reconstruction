"""
Tree visualization.
"""

import os
import pickle
from scipy.sparse import load_npz
from anndata import AnnData
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.phylo import *
from mito_utils.diagnostic_plots import *
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
        .loc[obs_tree.internal_nodes]
    )

    return afm, obs_tree, df_supports, variants


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results/trees')
path_top_trees = os.path.join(path_results, 'top_trees')

# Colors

##


# Read data
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
df = pd.read_csv(os.path.join(path_results, 'output/supports_df.csv'), index_col=0)
muts_df = pd.read_csv(os.path.join(path_data, 'muts', 'output_sAML1_new.txt'), sep='\t', index_col=0)
# muts_df = pd.read_csv(os.path.join(path_data, 'muts', 'genotypes.csv'), index_col=0)
# muts_df = pd.read_csv(os.path.join(path_data, 'muts', 'genotype_all.csv'), index_col=0).set_index('barcode')
muts_df.index = muts_df.index.map(lambda x: f'{x}-sAML1')
# muts_df['genotype'] = np.where(muts_df['genotype']=='mutated', 1,0)
# muts_df = muts_df.reset_index().pivot_table(index='barcode', values='genotype', columns='gene')


##
        

# Read colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'rb') as f:
    colors = pickle.load(f)


##
    

# Prep data
sample = 'sAML1'
test_f = lambda df,s: (df['sample_id']==s) # & \
                      # (df['malignant_class_occupancy']=='malignant') & \
                      # (df['aggregated_ct'].isin(['Early_myeloid', 'Late_myeloid', 'HSC_MPP']))
path_sample = os.path.join(path_top_trees, sample)

# Get sample objs
afm, obs_tree, df_supports, variants = prep_sample(sample)

# Filter
common = set(afm.obs_names) & set(meta.loc[test_f(meta, sample)].index) & set(muts_df.index)
to_remove = set(afm.obs_names)-common
common = list(common)
to_remove = list(to_remove)

# Get tree 
afm_filtered = nans_as_zeros(afm[common])
obs_tree = build_tree(afm_filtered)

# Add signatures and MT-SNVS to obs_tree cell_meta
obs_tree.cell_meta = (
    obs_tree.cell_meta
    .join(muts_df)
    #.join(
    #    pd.DataFrame(afm.X, columns=afm.var_names, index=afm.obs_names)
    #)
)

# sAML1 only
obs_tree.cell_meta['FLT3'].loc[lambda x: x=='-'] = 'unknown'

# obs_tree.cell_meta.groupby(['malignant_class_occupancy', 'FLT3']).size()
# obs_tree.cell_meta.query('malignant_class_occupancy=="tme" and FLT3=="WT"').index.to_list()

fig, ax = plt.subplots(figsize=(5,6.5))

plot_tree(
    obs_tree, ax=ax, 
    meta=['malignant_class_occupancy'], #'FLT3'],
    categorical_cmap_annot=colors['malignant_class'],#{'MUT':'r', 'WT':'b', 'unknown':'grey'},
    colorstrip_width=4,
    colorstrip_spacing=0.001,
    extend_branches=True,
    orient=90,
    internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
    internal_node_vmin=.5, internal_node_vmax=.9, 
    leaf_kwargs={'markersize':0},
)

fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, f'sAML1_FLT3.pdf'), dpi=500)


##


for x in muts_df.columns:
    obs_tree.cell_meta.groupby('malignant_class_occupancy')[x].value_counts()