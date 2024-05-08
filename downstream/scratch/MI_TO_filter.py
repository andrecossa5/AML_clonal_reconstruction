#!/usr/bin/python

import os
import argv
from itertools import product
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
# path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
# path_data = os.path.join(path_main, 'data')
path_main = argv[1]


##


# Read AFM and meta
# sample = 'sAML1'
sample = argv[2]

# Read, format and retain good cells
afm = read_one_sample(path_data, sample)

# Read cells meta
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
meta = meta.query('sample_id==@sample')
meta.index = meta.index.map(lambda x: x.split('-')[0])
afm.obs = afm.obs.join(meta[['malignant_class_occupancy']])

# Filter cells
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

# Product
jobs = list(product(*filtering_kwargs.values()))
jobs = [ dict(zip(filtering_kwargs.keys(), j)) for j in jobs ]

i = 0
for j in jobs:
    if d['low_confidence_af'] > d['high_confidence_af']:
        print(i)
        i += 1
        vars_df, dataset_df, a = filter_cells_and_vars(
            afm, filtering='weng2024', filtering_kwargs=j,
            lineage_column='malignant_class_occupancy',
            spatial_metrics=False, fit_mixtures=False, 
            path_priors='/Users/IEO5505/Desktop/AML_clonal_reconstruction/data/vars_df/priors.csv'
        )

# Save




























##


# Distances, kNNs
# df_weights = pd.read_csv('/Users/IEO5505/Desktop/AML_clonal_reconstruction/data/vars_df/priors.csv', index_col=0)
# w = 1-df_weights.loc[a.var_names,'median_prevalence'].values
# D = pair_d(a, weights=w)
# 
# tree_uw = build_tree(a)
# tree_w = build_tree(a, weights=w)
# fig, axs = plt.subplots(1,2,figsize=(10,5))
# plot_tree(tree_uw, ax=axs[0])
# plot_tree(tree_w, ax=axs[1])
# fig.tight_layout()
# plt.show()



##








# from mito_utils.clustering import *
 
# D = np.corrcoef(a.X.T>=0.01)
# D[np.isnan(D)] = 0  
# D[np.diag_indices(D.shape[0])] = 1
# linkage_matrix = linkage(D, method='weighted')
# order = leaves_list(linkage_matrix)

# np.percentile(D.flatten(), 99)

# fig, axs = plt.subplots(1,2)
# axs[0].imshow(D[np.ix_(order, order)])
# # axs[1].imshow(D[np.ix_(idx, idx)])
# fig.tight_layout()
# plt.show()



# Get clade muts
# clades = { 
#     x : (get_internal_node_muts(tree, x), tree_collapsed.leaves_in_subtree(x)) \
#     for x in tree_collapsed.internal_nodes if x != 'root'
# }
# 
# # Quantify the prevalence ratio (clade vs rest mutation prevalence) of mutations assigned 
# # to each internal node. Get the most supported value and mutation
# stats = []
# for c in clades:
#     top_mut = assess_internal_node_muts(a, clades, c, high_af=.01)
#     s = top_mut.iloc[0,:]
#     stats.append(s)
# 
# df_stats = pd.concat(stats, axis=1).T.reset_index().rename(columns={'index':'mut'}).set_index('clade')
# final_muts = df_stats['mut'].unique()


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
# fig, ax = plt.subplots(figsize=(7,7))
# plot_tree(
#     tree_collapsed, 
#     ax=ax, orient=90, extend_branches=True,
#     leaf_kwargs={'markersize':5},
#     internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
#     cov_leaves='malignant_class_occupancy', cmap_leaves={'malignant':'r', 'tme':'b'}
# )
# fig.tight_layout()
# plt.show()      


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



