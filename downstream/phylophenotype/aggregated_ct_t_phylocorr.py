"""
Viz inferred transition rates.
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


def place_anno(df, colors, ax, orientation, pos, reverse=False):
    axins = ax.inset_axes(pos) 
    idx = df.index if not reverse else df.index [::-1]
    cmap = matplotlib.colors.ListedColormap([ colors[x] for x in idx ])
    cb = plt.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap), 
        cax=axins, orientation=orientation
    )
    cb.ax.set(xticks=[], yticks=[])


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results/trees')
path_top_trees = os.path.join(path_results, 'top_trees')

# Colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'rb') as f:
    colors = pickle.load(f)

##


# Read data
key = 'tme_aggregated_ct_t_rates' # 'malignant_aggregated_ct_t_rates' # 'tme_aggregated_ct_t_rates' #
sample = 'sAML1'

# Map and phylocorr
mapping = (
    pd.read_csv(
        os.path.join(path_top_trees, sample, f'mapping_{key}.csv'),
        index_col=0
    )['to_numeric'].to_dict()
)
df = pd.read_csv(os.path.join(path_top_trees, sample, f'phylocorr_{key}.csv'), index_col=0)
df['GBC1'] = df['GBC1'].map({v:k for k,v in mapping.items()})
df['GBC2'] = df['GBC2'].map({v:k for k,v in mapping.items()})
df_ = df.pivot_table(index='GBC1', columns='GBC2', values='Z')


##


# Viz 
fig, ax = plt.subplots(figsize=(5.5,5.5))

vmin = -10
vmax = 10
order = cluster_df(df_)
df_sorted = df_.iloc[order, order]
im = ax.imshow(df_sorted, cmap='inferno', interpolation='nearest', vmin=vmin, vmax=vmax)
orientation = 'horizontal'
pos = (0, 1, 1, 0.04)
place_anno(df_sorted, colors['cell_type'], ax=ax, orientation=orientation, pos=pos)
orientation = 'vertical'
pos = (1, .0, 0.04, 1)
place_anno(df_sorted, colors['cell_type'], ax=ax, orientation=orientation, pos=pos, reverse=True)
format_ax(ax=ax, xlabel='', ylabel='', xticks=df_sorted.columns, yticks=df_sorted.index, rotx=90)
add_cbar(
    df_sorted.values.flatten(), 
    palette='inferno', ax=ax, label='z-score phylogenetic correlation', 
    layout=( (1.15,0,.05,1), 'right', 'vertical' ), label_size=10, ticks_size=8, vmin=vmin, vmax=vmax

)
fig.subplots_adjust(left=.25, right=.7, top=.9, bottom=.25)
fig.savefig(os.path.join(path_top_trees, sample, f'phylocorr_{key}.pdf'), dpi=500)


##
