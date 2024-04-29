"""
Expression data visualization.
"""

import os
import scanpy as sc
from plotting_utils._utils import *
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'expression')


##


# Read data
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
adata = sc.read(os.path.join(path_data, 'expression', 'full_adata.h5ad'))
markers = pd.read_csv(os.path.join(path_data, 'expression', 'aggr_lin_markers.csv'), index_col=0)
markers = markers.rename(columns={'pct.1':'perc_1', 'pct.2':'perc_2'}) # R to python convention

# Filter markers for plotting
n = 3
logFC = 2
perc_1 = .5
perc_2 = .2
genes = {}
group = markers['cluster'].unique()
for group in group:
    genes[group] = (
        markers.query('cluster==@group and avg_log2FC>=@logFC and p_val_adj<=0.01 and perc_1>=@perc_1 and perc_2<=@perc_2')
        .sort_values('avg_log2FC', ascending=False)
        .index[:n]
        .to_list()
    )
from itertools import chain
genes_ = list(chain.from_iterable([ genes[x] for x in genes ]))

# Reformat due to anndata duplication
correction_d = {}
for x in genes_:
    if x in adata.var_names:
        correction_d[x] = x
    else:
        correction_d[x] = x[:-1]
new_index = []
for x in adata.var_names:
    if x not in correction_d:
        new_index.append(x)
    else:
        new_index.append(correction_d[x])
adata.var_names = new_index

# Calc pct e logFC stats
genes_ = list(correction_d.values())
df = pd.DataFrame(adata[:,genes_].X, index=adata.obs_names, columns=genes_).join(adata.obs[['aggregated_ct']])

pct = (
    df
    .groupby('aggregated_ct')
    .apply(lambda x: (x>0).sum(axis=0))
    .div(
        df.groupby('aggregated_ct').size(), axis=0
    )
    .T.drop_duplicates().T
    .reset_index()
    .melt(id_vars='aggregated_ct', var_name='gene', value_name='pct')
)
L = []
for group in pct['aggregated_ct'].unique():
    l = []
    df_group = df.query('aggregated_ct==@group')
    df_other = df.query('aggregated_ct!=@group')
    for gene in genes_:
        x = df_group[gene].values.mean() + 0.0001
        y = df_other[gene].values.mean() + 0.0001
        l.append(np.log2(x)-np.log2(y))
    L.append(l)
logFC = (
    pd.DataFrame(
        L, columns=genes_, index=pct['aggregated_ct'].unique()
    )
    .reset_index()
    .rename(columns={'index':'aggregated_ct'})
    .melt(id_vars='aggregated_ct', var_name='gene', value_name='log2FC')
)

# Here we go
df_stats = logFC.merge(pct, on=['aggregated_ct', 'gene'])
df_stats['aggregated_ct'] = pd.Categorical(df_stats['aggregated_ct'], categories=markers['cluster'].unique())


##


# Viz
fig, ax = plt.subplots(figsize=(7,4.5))
sns.scatterplot(
    data=df_stats, y='aggregated_ct', x='gene', size='pct', hue='log2FC', 
    palette='mako', ax=ax, sizes=(1, 70)
)
format_ax(ax, title='Markers', xlabel='Top 3 marker genes', ylabel='Cell types', rotx=90, xticks_size=7)
ax.legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'dot_plot.pdf'), dpi=300)


##