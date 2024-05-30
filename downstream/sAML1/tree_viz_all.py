"""
sAML1 trees visualization.
"""

import os
import pickle
from mito_utils.preprocessing import *
from mito_utils.phylo_plots import *
from mito_utils.embeddings_plots import *
from mito_utils.plotting_base import *
from mito_utils.dimred import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')
path_phylo = os.path.join(path_main, 'results', 'phylo_sAML1_all')
path_results = os.path.join(path_main, 'results', 'sAML1')


##


# Read MT_data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = meta.query('sample_id==@sample')
afm.obs = afm.obs.join(meta)

# Filter leukemic cells with MT-variant
_, _, lineages_df = filter_cells_and_vars(
    afm, sample_name=sample, path_priors=path_priors, 
    max_prior=.05, 
    filtering='GT_enriched',
    af_confident_detection=0.05,
    lineage_column='aggregated_ct',
    with_clones_df=True
)

from itertools import chain
variants = list(chain.from_iterable(lineages_df['vars'].map(lambda x: x.split('; ')).to_list()))
variants = [ x for x in variants if x != '' ]

_, a = filter_cells_and_vars(
    afm, sample_name=sample,
    variants=variants,
    af_confident_detection=0.05
)




##


# Colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'rb') as f:
    colors = pickle.load(f)
colors_malignant = { k:v for k,v in colors['malignant_class'].items() if k != 'hBM'}
colors_malignant[np.nan] = 'grey'


##


# Load tree and supports
# path_trees = os.path.join(path_phylo, sample, 'n5_5_optimized', 'NJ', 'jaccard', 'feature_resampling')
# with open(os.path.join(path_trees, 'trees.pickle'), 'rb') as f:
#     trees = pickle.load(f)
# 
# # Get observed tree and add meta
# obs_tree = trees['observed']
# obs_tree.cell_meta = (
#     a.obs.loc[obs_tree.leaves]
#     .join(pd.DataFrame(a[obs_tree.leaves,:].X, index=obs_tree.leaves, columns=a.var_names))
# )

# Load nodes supports
# df_supports = pd.read_csv(os.path.join(path_trees, 'extended_supports.csv'), index_col=0)


# Viz tree
D = pair_d(a.X, t=.05, weights=1-a.var['prior'].values, metric='jaccard')

D



# 
from mito_utils.clustering import *
# 
order = leaves_list(linkage(D))
# 
plt.imshow(D[np.ix_(order, order)], cmap='viridis_r')
plt.show()



obs_tree = build_tree(a, weights=1-a.var['prior'].values, metric='cosine', t=.05, solver='UPMGA')
obs_tree.cell_meta = (
    a.obs.loc[obs_tree.leaves]
    .join(pd.DataFrame(a[obs_tree.leaves,:].X, index=obs_tree.leaves, columns=a.var_names))
)
obs_tree.collapse_mutationless_edges(True)

fig, ax = plt.subplots(figsize=(5,5))
plot_tree(obs_tree, ax=ax, meta=['malignant_class_occupancy'], categorical_cmap_annot=colors_malignant, colorstrip_width=1)
fig.tight_layout()
plt.show()


a




# Here we go
results = {}
n_samples = 100
n_vars_sampling = round((a.shape[1] / 100) * 90)

L = []
for _ in range(n_samples):
    vars_ = np.random.choice(a.var_names, size=n_vars_sampling, replace=False)
    D = pair_d(a[:,vars_].X, t=.05, weights=1-a[:,vars_].var['prior'].values, metric='jaccard')
    L.append(D.flatten())

np.mean(np.corrcoef(np.array(L)))







# Filter leukemic cells and MT variants
_, a = filter_cells_and_vars(
    afm, sample_name=sample, 
    path_priors=path_priors, 
    max_prior=.1, 
    filtering='MI_TO',
    filtering_kwargs={
        "min_n_confidently_detected" : 2,
        "af_confident_detection" : 0.05,
        "min_frac_negative" : 0
    },
    af_confident_detection=.05,
    fit_mixtures=True,
    only_positive_deltaBIC=True,
    max_AD_counts=3,
    lineage_column='malignant_class_occupancy'
)
obs_tree = build_tree(a, weights=1-a.var['prior'].values, metric='jaccard', t=.05, solver='NJ')

