"""
Clones specific gene-expression patterns using HVGs.
"""

import os
import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
import textalloc as ta
from sklearn.preprocessing import scale
from scipy.sparse import csr_matrix
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results/trees')
path_top_trees = os.path.join(path_results, 'top_trees')

# Read data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
# signatures = pd.read_csv(os.path.join(path_data, 'expression', 'malignant_AUCell.csv'), index_col=0)
# HVGs = pd.read_csv(os.path.join(path_top_trees, sample, f'HVGs_malignant_unsupervised.csv'), index_col=0)
counts = pd.read_csv(os.path.join(path_data, 'expression', 'sAML1_malignant_lognorm_counts.csv'), index_col=0)
clones = pd.read_csv(os.path.join(path_top_trees, sample, 'top_clones.csv'), index_col=0)

# Build anndata
adata = AnnData(X=csr_matrix(counts.values), obs=clones, var=pd.DataFrame(index=counts.columns))
# adata = AnnData(X=csr_matrix(HVGs.values), obs=clones, var=pd.DataFrame(index=HVGs.columns))
# adata = AnnData(X=csr_matrix(signatures.values), obs=clones, var=pd.DataFrame(index=signatures.columns))


##


# Expression pp and viz
# sc.pp.highly_variable_genes(adata)
# red = adata[:,adata.var['highly_variable']].copy()
# sc.pp.scale(red)
# sc.pp.pca(red, n_comps=30)
# sc.pp.neighbors(red, n_neighbors=10)
# sc.tl.umap(red)
# sc.pl.umap(red, color=clones.columns)


## 


# DE
test_genes = (adata.X.A>0).sum(axis=0)>.05*adata.shape[0]
adata = adata[:,test_genes].copy()
DE_clones = {}
for clone in clones.columns:
    sc.tl.rank_genes_groups(adata, method='wilcoxon', groupby=clone)
    df_res = (
        pd.DataFrame({
            'gene' : adata.uns['rank_genes_groups']['names']['clone'],
            'logFC' : adata.uns['rank_genes_groups']['logfoldchanges']['clone'],
            'pval_adj' : adata.uns['rank_genes_groups']['pvals_adj']['clone'],
            'scores' : adata.uns['rank_genes_groups']['scores']['clone'],
        })
        .sort_values('logFC', ascending=False)
    )
    if df_res.query('pval_adj<=0.1').shape[0]>50:
        DE_clones[clone] = df_res

# Viz DE
DE_clones.keys()
clone = 'clone_5'
df_res = DE_clones[clone].set_index('gene')

# Volcano
fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot(df_res.query('logFC>=1 and pval_adj<=0.1')['logFC'],
        -np.log10(df_res.query('logFC>=1 and pval_adj<=0.1')['pval_adj']), 'r.', markersize=5)
ax.plot(df_res.query('logFC<=-1 and pval_adj<=0.1')['logFC'],
        -np.log10(df_res.query('logFC<=-1 and pval_adj<=0.1')['pval_adj']), 'b.', markersize=5)
test = ( (df_res['logFC']<1) & (df_res['logFC']>-1) ) | (df_res['pval_adj']>0.1)
ax.plot(df_res.loc[test]['logFC'], -np.log10(df_res.loc[test]['pval_adj']), 'k.', markersize=2.5)
ax.set(xlim=(-2.5,2.5))
format_ax(ax, title=f'{clone} vs rest', xlabel='logFC', ylabel='-log10(pvalue)', reduced_spines=True)
test = ( (df_res['logFC']<1.5) & (df_res['logFC']>-1.5) ) | (df_res['pval_adj']>0.01)
ta.allocate_text(
    fig, ax, 
    df_res.loc[~test]['logFC'],
    -np.log10(df_res.loc[~test]['pval_adj']),
    df_res.loc[~test].index,
    x_scatter=df_res['logFC'], y_scatter=df_res['pval_adj'], 
    linecolor='black', textsize=8, 
    max_distance=.2, linewidth=0.5, nbr_candidates=100
)

# Save
fig.tight_layout()
plt.show()

# fig.savefig(os.path.join(path_top_trees, sample, f'{clone}_DEGs.png'), dpi=500)


##