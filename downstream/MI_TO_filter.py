"""
New, redeem-like MI_TO filter, on streoids.
"""

import os
from mito_utils.preprocessing import *
from mito_utils.kNN import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
from mito_utils.diagnostic_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')
import warnings
warnings.simplefilter('ignore')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_colors = os.path.join(path_data, 'meta', 'colors.pickle')
path_results = os.path.join(path_main, 'results', 'figure_out_sAML1')


##


# Read AFM and meta
sample = 'sAML1'

# Read, format and retain good cells
afm = read_one_sample(path_data, sample, with_GBC=False)
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
meta =  meta.query('sample_id=="sAML1"')
meta.index = meta.index.map(lambda x: x.split('-')[0])
afm.obs = afm.obs.join(meta)
afm = afm[~afm.obs['malignant_class_occupancy'].isna(),:].copy()


##


# Filter now
a = filter_baseline(afm)
a = filter_baseline_old(afm)
# a = filter_sites(a)
a = filter_density(afm, density=.07)
a = filter_weng2024(afm, make_vars_df(afm))


from mito_utils.dimred import reduce_dimensions

df = reduce_dimensions(a, method='UMAP', metric='jaccard', n_comps=2)
X = df[0][:,:2]
fig, ax = plt.subplots(figsize=(7,7))
ax.plot(X[:,0], X[:,1], 'ko')
fig.tight_layout()
plt.show()    


tree = build_tree(a, metric='jaccard')

fig, ax = plt.subplots(figsize=(7,7))
plot_tree(tree, ax=ax)
fig.tight_layout()
plt.show()    


low_af = .01
afm.uns['per_position_coverage'].mean(axis=0).describe()
X_bin = np.where(a.X>=low_af,1,0)


##


# 1. n cells per var and n vars per cell
# vars_df.loc[vars_, 'Variant_CellN'].describe()
pd.Series(np.sum(X_bin, axis=1)).describe()
pd.Series(np.sum(X_bin, axis=0)).describe()


# 2. Connectedness
D = pairwise_distances(X_bin, metric=lambda x, y: np.sum(np.logical_and(x, y)))
cell_conn = np.ma.masked_equal(D, np.diag(D)).mean(axis=1).data
pd.Series(cell_conn).describe()

# 3. Sparseness and n genotypes occurrence
np.sum(a.X==1) / np.product(a.shape)
d = AFM_to_seqs(a)
pd.Series(d).value_counts().describe()

# 4. Annot
a.var_names.map(lambda x: x.split('_')[1]).value_counts() # .sum()

# 4. Drop-out simulation
#...

# 5. TME-leukemic clone?
tme = a.obs['malignant_class_occupancy']=='tme'
tme_prevalences = np.sum(X_bin[tme,:], axis=0)/tme.sum()
leukemic = a.obs['malignant_class_occupancy']=='malignant'
leukemic_prevalences = np.sum(X_bin[leukemic,:], axis=0)/leukemic.sum()

np.sum((tme_prevalences<=.1) & (leukemic_prevalences>=.2))
corr = np.corrcoef(leukemic_prevalences, tme_prevalences)[0,1]

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot(tme_prevalences, leukemic_prevalences, 'ko', markersize=3)
sns.regplot(x=tme_prevalences, y=leukemic_prevalences, ax=ax, scatter=False)
format_ax(ax, title=f'TME-malignant correlation: r2={corr:.2f}', 
          xlabel='TME prevalence', ylabel='Malignant prevalence', reduced_spines=True)
ax.set_xlim((-0.05,1.05))
ax.set_ylim((-0.05,1.05))
fig.tight_layout()
plt.show()



# 6. Tree mutations support

# Post-process tree
# tree_collapsed = tree.copy()
# tree_collapsed.collapse_mutationless_edges(True)
# len(tree_collapsed.internal_nodes) / len(tree.internal_nodes)

# Get clade muts
clades = { 
    x : (get_internal_node_muts(tree, x), tree.leaves_in_subtree(x)) \
    for x in tree.internal_nodes if x != 'root'
}

# Quantify the prevalence ratio (clade vs rest mutation prevalence) of mutations assigned 
# to each internal node. Get the most supported value and mutation
stats = []
for c in clades:
    top_mut = assess_internal_node_muts(a, clades, c, high_af=.01)
    s = top_mut.iloc[0,:]
    stats.append(s)

df_stats = pd.concat(stats, axis=1).T.reset_index().rename(columns={'index':'mut'}).set_index('clade')
final_muts = df_stats['mut'].unique()


##


# Rebuild with filtered muts

# Build tree
_, a = filter_cells_and_vars(afm, variants=final_muts)
X_bin = np.where(a.X>=high_af,1,0)

# Tree
tree = build_tree(a, t=high_af, solver='NJ')
tree_collapsed = tree.copy()
tree_collapsed.collapse_mutationless_edges(True)
len(tree_collapsed.internal_nodes) / len(tree.internal_nodes)

# Format into df
clades = { 
    x : (get_internal_node_muts(tree_collapsed, x), tree_collapsed.leaves_in_subtree(x)) \
    for x in tree_collapsed.internal_nodes if x != 'root'
}
stats = []
for c in clades:
    top_mut = assess_internal_node_muts(a, clades, c, high_af=high_af)
    s = top_mut.iloc[0,:]
    stats.append(s)

# Final supporting muts stats
df_stats = pd.concat(stats, axis=1).T.reset_index().rename(columns={'index':'mut'}).set_index('clade')
final_muts = df_stats['mut'].unique()
final_muts.size

df_stats.sort_values('ncells', ascending=False)


##


# Viz individual muts supporting some clades
# fig, ax = plt.subplots(figsize=(4,4.5))
# sns.kdeplot(a[cells, top_mut].X.flatten(), ax=ax, fill=True)
# sns.kdeplot(a[other_cells, top_mut].X.flatten(), ax=ax, fill=True)
# median_clade = np.median(a[cells, top_mut].X.flatten())
# median_rest = np.median(a[other_cells, top_mut].X.flatten())
# format_ax(
#     ax, title=f'{top_mut}: \n Median AF clade={median_clade:.2f} \n median AF rest={median_rest:.2f}', 
#     yticks=[], xlabel='AF'
# )
# for x in [median_clade, median_rest]:
#     ax.axvline(x, color='grey', linestyle='--')
# ax.spines[['right', 'left', 'top']].set_visible(False)
# fig.tight_layout()
# plt.show()


##


# Viz tree
fig, ax = plt.subplots(figsize=(7,7))
plot_tree(
    tree_collapsed, 
    ax=ax, orient=90, extend_branches=True,
    leaf_kwargs={'markersize':5},
    internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
    cov_leaves='malignant_class_occupancy', cmap_leaves={'malignant':'r', 'tme':'b'}
)
fig.tight_layout()
plt.show()      


##


























































# # Filter cellsls
# afm = filter_cells_coverage(afm, mean_coverage=50)
# filter_baseline(afm).var_names
# 
# afm.uns['per_position_quality'].isna().values.all()
# 
# 
# 
# # Shallow filter ad hoc:
# test_sites = pd.Series(
#     np.mean(afm.uns['per_position_coverage'], axis=0) > 10,
#     index=afm.uns['per_position_coverage'].columns
# )
# sites = test_sites[test_sites].index
# test_vars_site_coverage = (
#     afm.var_names
#     .map(lambda x: x.split('_')[0] in sites)
#     .to_numpy(dtype=bool)
# )
# test_vars_site_coverage.sum()
# 
# # Test 2-4: 
# # variants seen in at least 2 cells;
# # variants with AF>0.01 in at least 2 cells;
# test_vars_quality = np.nanmean(afm.layers['quality'], axis=0) > 30
# 
# np.nanmin(afm.layers['quality'])
# 
# afm.layers['quality'][:,2000:2020]
# 
# test_vars_coverage = np.sum(afm.layers['coverage']>0, axis=0) > 2
# test_vars_AF = np.sum(afm.X > 0.01, axis=0) > 2
# 
# # Filter vars and sites
# test_vars = test_vars_site_coverage & test_vars_quality & test_vars_coverage & test_vars_AF 
# filtered = afm[:, test_vars].copy()
# filtered = remove_excluded_sites(filtered)