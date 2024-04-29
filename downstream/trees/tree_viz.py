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

# Add covariates
meta['FLT_with_malignant'] = meta['malignant_class_occupancy'] + '_' + meta['FLT3'].astype('str')


##
        

# Read colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'rb') as f:
    colors = pickle.load(f)


##


# sAML1 only
fig, ax = plt.subplots(figsize=(5,6.5))

sample = 'sAML1'
afm, obs_tree, df_supports, variants = prep_sample(sample)
to_remove = afm.obs_names[~afm.obs_names.isin(meta.query('sample_id==@sample').index)]
obs_tree.remove_leaves_and_prune_lineages(to_remove)
df_supports = df_supports.loc[obs_tree.internal_nodes]

# df_supports.query('time<10')
# cells_ = get_clades(obs_tree)['cassiopeia_internal_node21']
# obs_tree.cell_meta.loc[cells_]['aggregated_ct'].value_counts()

plot_tree(
    obs_tree, ax=ax, 
    meta=['malignant_class_occupancy'],
    categorical_cmap_annot=colors['malignant_class'],
    colorstrip_width=4,
    colorstrip_spacing=0.001,
    extend_branches=True,
    orient=90,
    meta_internal_nodes=df_supports,
    cov_internal_nodes='TBE',
    internal_node_kwargs={'markersize':4, 'markeredgecolor':None},
    internal_node_vmin=.5, internal_node_vmax=.9, 
    leaf_kwargs={'markersize':0},
)
format_ax(
    ax=ax, 
    title=
    f'''
    {sample} \n TBE {df_supports['TBE'].median():.2f} (+-{df_supports['TBE'].std():.2f}), FBP {df_supports['FBP'].median():.2f} (+-{df_supports['FBP'].std():.2f})
    ''',
    title_size=11
)
add_legend(
    label='Cell status', 
    colors={ k:c for k,c in colors['malignant_class'].items() \
            if k in obs_tree.cell_meta['malignant_class_occupancy'].values}, 
    ax=ax,
    loc='center', artists_size=8, label_size=10, ticks_size=8, 
    bbox_to_anchor=(.25,-.1), ncols=2
)
add_cbar(
    df_supports['TBE'], palette='Spectral_r', ax=ax, label='TBE',
    ticks_size=8, label_size=10, 
    vmin=.5, vmax=.8, 
    layout=( (1-.5,-0.13,.4,.02), 'top', 'horizontal' )
    #(1-.8,-0.13,.4,.02)
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'sAML1.pdf'), dpi=500)


##


# Viz 4 of them, with malignancy class
vmin_annot = .001
vmax_annot = .1

fig, axs = plt.subplots(1,4,figsize=(14.5,5))

samples = ['sAML1', 'AML2', 'AML3', 'AML5']
for i, (ax, sample) in enumerate(zip(axs, samples)):
    afm, obs_tree, df_supports, variants = prep_sample(sample)
    ro_remove = afm.obs_names[~afm.obs_names.isin(meta.query('sample_id==@sample').index)]
    obs_tree.remove_leaves_and_prune_lineages(ro_remove)
    df_supports = df_supports.loc[obs_tree.internal_nodes]
    plot_tree(
        obs_tree, ax=ax, 
        meta=['malignant_class_occupancy'],
        categorical_cmap_annot={ k:c for k,c in colors['malignant_class'].items() \
                                if k in obs_tree.cell_meta['malignant_class_occupancy'].values},
        colorstrip_width=4,
        colorstrip_spacing=0.001,
        extend_branches=True,
        orient=90,
        meta_internal_nodes=df_supports,
        cov_internal_nodes='TBE',
        internal_node_kwargs={'markersize':4, 'markeredgecolor':None},
        internal_node_vmin=.5, internal_node_vmax=.9, 
        leaf_kwargs={'markersize':0},
    )
    format_ax(
        ax=ax, 
        title=
        f'''
        {sample} \n TBE {df_supports['TBE'].median():.2f} (+-{df_supports['TBE'].std():.2f}), FBP {df_supports['FBP'].median():.2f} (+-{df_supports['FBP'].std():.2f})
        ''',
        title_size=11
    )
    if i == 0:
        add_legend(
            label='Cell status', 
            colors={ k:c for k,c in colors['malignant_class'].items() \
                                if k in obs_tree.cell_meta['malignant_class_occupancy'].values}, 
            ax=ax,
            loc='center', artists_size=8, label_size=10, ticks_size=8, 
            bbox_to_anchor=(.25,-.1), ncols=2
        )
        add_cbar(
            df_supports['TBE'], palette='Spectral_r', ax=ax, label='TBE',
            ticks_size=8, label_size=10, 
            vmin=.5, vmax=.8, 
            layout=( (1-.4,-0.13,.4,.02), 'top', 'horizontal' )
        )

# Format and save
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'top_trees.pdf'), dpi=500)


##


# Viz
# fig, axs = plt.subplots(1,3,figsize=(13,5))
# 
# pvalues = [2.433749*10**(-27), 0.2731101, 0.1289957]
# for i,gene in enumerate(['RUNX1', 'CSF3R', 'FLT3']):
#     ax = axs[i]
#     colors = { 'wild-type':sns.color_palette('inferno')[0], 'mutated':sns.color_palette('inferno')[-2]}
#     plot_tree(
#         TREES[gene], ax=ax, meta=[gene], 
#         categorical_cmap_annot=colors, 
#         colorstrip_width=5,
#     )
#     ax.set(title=f'{gene} \n phylogenetic corr pvalue: {pvalues[i]:.2e}')
# 
# fig.tight_layout()
# fig.savefig(os.path.join(path_results, 'expressed_muts.png'), dpi=500)


##