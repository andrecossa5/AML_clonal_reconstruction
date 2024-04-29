"""
Viz signatures correlation patterns
"""

import os
import pickle
from scipy.cluster import hierarchy
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
from mito_utils.phylo_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


def cluster_df(df):
    D = pairwise_distances(df, metric='hamming') * df.shape[1]
    D = D-D.min() / D.std()
    linkage_matrix = hierarchy.linkage(squareform(D), method='weighted')
    order = hierarchy.leaves_list(linkage_matrix)
    return order


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results/trees')
path_top_trees = os.path.join(path_results, 'top_trees')

# Read data
key = 'malignant_supervised'
sample = 'sAML1'

# Read and reformat signatures phylocorr
df = pd.read_csv(os.path.join(path_top_trees, sample, f'phylocorr_{key}.csv'), index_col=0)
df_ = df.pivot_table(index='GBC1', columns='GBC2', values='Z')


##


# Viz 
fig, ax = plt.subplots(figsize=(6.5,5.5))

# vmin = -10
# vmax = 10
order = cluster_df(df_)
df_sorted = df_.iloc[order, order]
im = ax.imshow(df_sorted, cmap='inferno', interpolation='nearest')#, vmin=vmin, vmax=vmax)
format_ax(ax=ax, xlabel='', ylabel='', title=sample,
          xticks=[], yticks=df_sorted.index, rotx=90, xticks_size=3, yticks_size=3)
add_cbar(
    df_sorted.values.flatten(), 
    palette='inferno', ax=ax, label='z-score phylogenetic correlation', 
    layout=( (1.15,0,.05,1), 'right', 'vertical' ), label_size=10, ticks_size=8,# vmin=vmin, vmax=vmax

)
fig.subplots_adjust(left=.4, right=.7, top=.9, bottom=.25)
fig.savefig(os.path.join(path_top_trees, sample, f'phylocorr_{key}.pdf'), dpi=500)


##



# Genes z-scores phylocorr visualization
samples = ['sAML1', 'AML2', 'AML3', 'AML5']


##


# Viz 
for sample in samples:

    order = np.diag(df_).argsort()[::-1]
    x = np.diag(df_)[order][:20]
    fig, ax = plt.subplots(figsize=(7,4))
    stem_plot(pd.DataFrame(x, columns=['z']), x='z', ax=ax)
    format_ax(ax, xlabel='Phylogenetic correlation z-score', yticks=df_.index[order][:20], yticks_size=7)
    fig.subplots_adjust(left=.6, right=.9, bottom=.2, top=.9)
    fig.savefig(os.path.join(path_top_trees, sample, f'malignant_supervised_rank.pdf'), dpi=500)


##
    