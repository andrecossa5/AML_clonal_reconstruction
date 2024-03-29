"""
Expression data visualization.
"""

import os
from plotting_utils._utils import *
from plotting_utils._plotting_base import *
from scipy.cluster import hierarchy
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
from matplotlib.colors import LinearSegmentedColormap
matplotlib.use('macOSX')


##


def cluster_df(df):
    D = pairwise_distances(df, metric='hamming') * df.shape[1]
    D = D-D.min() / D.std()
    linkage_matrix = hierarchy.linkage(squareform(D), method='weighted')
    order = hierarchy.leaves_list(linkage_matrix)
    return order


##


# Set paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/'
path_data = os.path.join(path_main, 'data', 'meta')
path_results = os.path.join(path_main, 'results', 'muts')


##


# Read data
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
muts = [
    'SRSF2', 'STAG2', 'FRYL', 'IDH1', 'EZH2', 'SLC35B2', 
    'CSF3R', 'RUNX1', 'SETD2', 'FLT3', 'HDAC7', 'MACF1',
]


##


# Genotype matrix
n_min = 8
cells = (~meta[muts].isna()).sum(axis=1).loc[lambda x:x>=n_min].index
df = pd.DataFrame(
    np.select(
        [ meta.loc[cells, muts].isna(), meta.loc[cells, muts]=='wild-type', meta.loc[cells, muts]=='mutated' ],
        [ 0, 1, 2 ]
    ),
    columns=muts, index=cells
).T

# Cluster cells
muts_order = cluster_df(df)
cells_order = cluster_df(df.T)

# Plot heatmap
fig, ax = plt.subplots(figsize=(6,5))
colors = ['lightgrey', sns.color_palette('inferno')[0], sns.color_palette('inferno')[-2]]
cmap = sns.color_palette(colors)
sns.heatmap(df.iloc[muts_order, cells_order], cmap=cmap, ax=ax, cbar=False)
ax.set_xticks([])
format_ax(ax, 
    title=f'sAML1-sample, {cells.size} cells\n SCM-seq mutation status',
    yticks_size=10
)
c = { k:v for k,v in zip(['NA', 'wt', 'mutated'], colors) }
add_legend(ax=ax, label='Mutational status', colors=c, loc='center', bbox_to_anchor=(.5, -.1),
           artists_size=10, label_size=10, ticks_size=8, ncols=3)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'SCM_seq_muts.png'), dpi=400)


##
