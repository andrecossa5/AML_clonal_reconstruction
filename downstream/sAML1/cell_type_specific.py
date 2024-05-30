"""
sAML1 MT-SNVs characterization
"""

import os
from mito_utils.clustering import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')



##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'sAML1')
path_var_selection = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')
# path_top_trees = os.path.join(path_results, 'top_trees')

# Read data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = meta.query('sample_id==@sample')
afm.obs = (
    afm.obs
    .join(
    meta[['nCount_RNA', 'percent_MT', 
          'malignant_class_occupancy', 'aggregated_ct']])
)

# Find mut
rank_clone_variants(afm, 'malignant_class_occupancy', 'malignant')


##


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


# JI matrix
X_ct = pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names).join(a.obs[['aggregated_ct']]).groupby('aggregated_ct').median()

D = pair_d(X_ct, metric='jaccard')
order = leaves_list(linkage(D, method='weighted'))

fig, ax = plt.subplots(figsize=(7,5))
ax.imshow(D[np.ix_(order, order)], cmap='viridis_r')
add_cbar(D.flatten(), palette='viridis_r', layout='outside', ax=ax, label='1-jaccard', label_size=12, ticks_size=10)
format_ax(ax, xticks=X_ct.index[order], yticks=X_ct.index[order], rotx=90, title=f'n={X_ct.shape[1]} cell-type-biased MT-SNVs')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'lineage_biased_muts.png'), dpi=500)



##


# Individual ct associations
fig, axs = plt.subplots(2,1,figsize=(4,8))

x = np.mean(a[a.obs['aggregated_ct']=='T_CD4'].X, axis=0)
y = np.mean(a[a.obs['aggregated_ct']=='T_CD8'].X, axis=0)
axs[0].plot(x, y, 'ko', alpha=.3)
corr = np.corrcoef(x,y)[1,0]
format_ax(axs[0], title=f'T-cells: CD4 vs CD8 \n Pearson\'s r: {corr:.2f}', xlabel='Mean AF in CD4 T cells', ylabel='Mean AF in CD8 T cells')

x = np.mean(a[a.obs['aggregated_ct']=='Mono'].X, axis=0)
y = np.mean(a[a.obs['aggregated_ct']=='DC'].X, axis=0)
axs[1].plot(x, y, 'ko', alpha=.3)
corr = np.corrcoef(x,y)[1,0]
format_ax(axs[1], title=f'T-cells: Monocytes vs Dendritic Cells \n Pearson\'s r: {corr:.2f}', xlabel='Mean AF in Monocytes', ylabel='Mean AF in Dendritic cells')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'indidivual_corr.png'), dpi=500)


##

