"""
Tree visualization.
"""

import os
import pickle
from scipy.sparse import load_npz, csr_matrix
from anndata import AnnData
import cassiopeia as cs
from sklearn.preprocessing import scale
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.phylo import *
from mito_utils.embeddings_plots import *
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

# Read data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
signatures = pd.read_csv(os.path.join(path_data, 'expression', 'malignant_AUCell.csv'), index_col=0)
# HVGs = pd.read_csv(os.path.join(path_top_trees, sample, f'HVGs_malignant_unsupervised.csv'), index_col=0)
counts = pd.read_csv(os.path.join(path_data, 'expression', 'sAML1_malignant_lognorm_counts.csv'), index_col=0)
# clones = pd.read_csv(os.path.join(path_top_trees, sample, 'top_clones.csv'), index_col=0)
df = pd.read_csv(os.path.join(path_results, 'output/supports_df.csv'), index_col=0)

# Read and reformat signatures phylocorr
key = 'malignant_supervised'
sample = 'sAML1'
df_supervised = pd.read_csv(os.path.join(path_top_trees, sample, f'phylocorr_{key}.csv'), index_col=0)


##
        

# Read colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'rb') as f:
    colors = pickle.load(f)


##


# sAML1 top expanded clones
sample = 'sAML1'
test_f = lambda df,s: (df['sample_id']==s) & \
                      (df['malignant_class_occupancy']=='malignant') & \
                      (df['aggregated_ct'].isin(['Early_myeloid', 'Late_myeloid', 'HSC_MPP']))
path_sample = os.path.join(path_top_trees, sample)

# Get sample objs
afm, obs_tree, df_supports, variants = prep_sample(sample)

# Filter
common = set(afm.obs_names) & set(meta.loc[test_f(meta, sample)].index)
to_remove = set(afm.obs_names)-common
common = list(common)
to_remove = list(to_remove)

# Get tree 
afm_filtered = nans_as_zeros(afm[common])
obs_tree = build_tree(afm_filtered)

# Add signatures and MT-SNVS to obs_tree cell_meta
obs_tree.cell_meta = (
    obs_tree.cell_meta
    .join(
        pd.DataFrame(
            scale(signatures.iloc[:,1:]),
            index=signatures.index, 
            columns=signatures.columns[1:]
        )
    )
    .join(
        pd.DataFrame(afm.X, columns=afm.var_names, index=afm.obs_names)
    )
)


## 


# Plot signatures

# To plot
plastic = df_supervised.query('p>=.5 and GBC1==GBC2')['GBC1'][-3:].to_list()
inheritable = (
    df_supervised
    .query('p<=0.001 and Z>=10 and GBC1==GBC2')
    .sort_values('Z', ascending=False)['GBC1'][:3].to_list()
)
to_plot = inheritable + plastic

# Viz top inheritable signatures
fig = plt.figure(figsize=(8,6))
for i, signature in enumerate(to_plot):
    ax = fig.add_subplot(2,3,i+1)
    x = obs_tree.cell_meta[signature]
    plot_tree(
        obs_tree, ax=ax, 
        meta=[to_plot[i]],
        continuous_cmap_annot='viridis',
        colorstrip_width=7,
        colorstrip_spacing=0.001,
        extend_branches=True,
        orient=90,
        vmin_annot=np.percentile(x, 5),
        vmax_annot=np.percentile(x, 95),
        internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
        internal_node_vmin=.5, internal_node_vmax=.9, 
        leaf_kwargs={'markersize':0},
    )
    corr = df_supervised.query('GBC1==GBC2 and GBC1==@signature')['Z'].values[0]
    format_ax(ax=ax, title=f'{signature}\n z-scored corr: {corr:.2f}', title_size=8)
fig.tight_layout()
fig.savefig(os.path.join(path_top_trees, sample, f'signatures_viz.pdf'), dpi=500)


##


# Top expanded clones

# Find expanded clones
cs.tl.compute_expansion_pvalues(obs_tree, min_clade_size=(0.05 * obs_tree.n_cell), min_depth=0)
probability_threshold = 0.05
expansions = {}
for node in obs_tree.depth_first_traverse_nodes():
    if obs_tree.get_attribute(node, "expansion_pvalue") <= probability_threshold:
        expansions[node] = obs_tree.get_attribute(node, "expansion_pvalue")
L = []
for node in expansions:
    clade = get_clades(obs_tree)[node]
    L.append([node, obs_tree.cell_meta.index.isin(clade).sum()])
expanded_df = (
    pd.DataFrame(L, columns=['node', 'ncells'])
    .sort_values('ncells', ascending=False)
    .set_index('node')
    .join(pd.Series(expansions).to_frame('p'))
    .query('ncells<=1000')
)
expanded = expanded_df.index.to_list()

# Add to meta
df_ = obs_tree.cell_meta.copy()
for i, node in enumerate(expanded):
    df_[f'clone_{i+1}'] = np.where(
        df_.index.isin(get_clades(obs_tree)[node]),
        'clone', 'other'
    )
obs_tree.cell_meta = df_


##


# Clones
to_plot = obs_tree.cell_meta.columns[obs_tree.cell_meta.columns.str.contains('clone')].to_list()

fig = plt.figure(figsize=(8,8))
for i, clone in enumerate(to_plot):
    ax = fig.add_subplot(3,3,i+1)
    x = obs_tree.cell_meta[clone]
    plot_tree(
        obs_tree, ax=ax, 
        colorstrip_width=4,
        colorstrip_spacing=0.001,
        extend_branches=False,
        orient=90,
        leaf_kwargs={'markersize':1},
        internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
        cov_leaves=clone,
        cmap_leaves={'clone':'r', 'other':'grey'}
    )
    ncells = expanded_df.iloc[i]['ncells']
    p = expanded_df.iloc[i]['p']
    format_ax(ax=ax, title=f'{clone.capitalize()} \n {ncells} cells, {p:.2e} pvalue', title_size=10)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_top_trees, sample, 'expanded_clones.pdf'), dpi=500)


# Fig
fig, ax = plt.subplots(figsize=(4,4))
clone = 'clone_3'
i = 2
x = obs_tree.cell_meta[clone]
plot_tree(
    obs_tree, ax=ax, 
    colorstrip_width=4,
    colorstrip_spacing=0.001,
    extend_branches=False,
    orient=90,
    leaf_kwargs={'markersize':1},
    internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
    cov_leaves=clone,
    cmap_leaves={'clone':'r', 'other':'grey'}
)
ncells = expanded_df.iloc[i]['ncells']
p = expanded_df.iloc[i]['p']
format_ax(ax=ax, title=f'{clone.capitalize()} \n {ncells} cells, {p:.2e} pvalue', title_size=10)
fig.tight_layout()
fig.savefig(os.path.join(path_top_trees, sample, 'clone_3.pdf'), dpi=500)


##


# Save expanded clones
obs_tree.cell_meta[to_plot].to_csv(os.path.join(path_top_trees, sample, 'top_clones.csv'))


##


# UMAP
clones = pd.read_csv(os.path.join(path_top_trees, sample, 'top_clones.csv'), index_col=0)
meta = obs_tree.cell_meta.copy()
# meta = meta.join(clones)

# Build adata sAML1
adata = AnnData(
    X=csr_matrix(counts.values), 
    obs=meta, 
    var=pd.DataFrame(index=counts.columns)
)

# Umap
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
red = adata[:,adata.var['highly_variable']].copy()
sc.pp.scale(red)
sc.pp.pca(red, n_comps=50)
sc.pp.neighbors(red, n_pcs=30, n_neighbors=10)
sc.tl.umap(red)

# Umap clones
df_ = meta.join(
    pd.DataFrame(
        red.obsm['X_umap'], 
        index=red.obs_names, 
        columns=['UMAP1','UMAP2']
    )
)

# Fig
fig, ax = plt.subplots(figsize=(5,5))
colors = {'clone':'r','other':'grey'}
draw_embeddings(
    df_, cat='clone_3', ax=ax, s=30,
    legend_kwargs={
        'colors':colors,
        'bbox_to_anchor' : (1,1),
        'loc' : 'upper right', 
        'label_size' : 12,
        'ticks_size' : 10,
    },
)
ax.set(title='')
ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(path_top_trees, sample, 'clone_3_UMAP.pdf'), dpi=500)


##