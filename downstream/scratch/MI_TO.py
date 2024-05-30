"""
Trial filter_cells_and_vars.
"""

import os
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from mito_utils.preprocessing import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
from mito_utils.plotting_base import *
from mito_utils.utils import ji
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')


##


# Read data
sample = 'sAML1'
path_results = os.path.join(path_results, sample)

# AFM
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
meta = meta.query('sample_id==@sample')
afm.obs = (
    afm.obs
    .join(
    meta[['nCount_RNA', 'percent_MT', 
          'malignant_class_occupancy', 'aggregated_ct']])
)






























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


