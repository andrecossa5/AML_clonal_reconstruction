"""
Prova analysis sAML1.
"""

import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from mito_utils.utils import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from matplotlib.gridspec import GridSpec
from Cellula.plotting._plotting import plot_heatmap
from Cellula.plotting._colors import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._pp_recipes import *
from Cellula.preprocessing._neighbors import *
from Cellula.preprocessing._embeddings import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/'
path_data = os.path.join(path_main, 'data', 'meta')
path_results = os.path.join(path_main, 'results')

# Args
sample = 'sAML1'
make_folder(path_results, sample)
path_results = os.path.join(path_results, sample)


##


# Read all data
adata_all = sc.read(os.path.join(path_data, 'cc_mito.h5ad'))
cell_assignment = pd.read_csv(
    os.path.join(path_results, '../vireo/', sample, 'cell_assignments.csv'),
    index_col=0
)
mito_clones_df = pd.read_csv(
    os.path.join(path_results, '../vireo/', sample, 'labels.csv'), 
    index_col=0
)
mito_umap = pd.read_csv(
    os.path.join(path_results, '../vireo/', sample, 'MT_SNVs_umap.csv'), 
    index_col=0
)


##


# Visualization and reports

# UMAP (expr) all samples
adata_all.obs['sample'] = adata_all.obs['sample'].astype(str)
adata_all.obsm['UMAP'].columns = ['UMAP1', 'UMAP2']
df_ = adata_all.obs.join(adata_all.obsm['UMAP'])
df_['SRSF2_genotype'] = df_['SRSF2_genotype'].astype(str)
df_['SRSF2_genotype'] = np.where(df_['SRSF2_genotype'] == 'SRSF2 mutated', 'mut', 'wt')

# Coarse-annotation + SRF
fig, axs = plt.subplots(1,2, figsize=(11, 5))

# SRSF2
draw_embeddings(
    df_, 
    cat='SRSF2_genotype', 
    ax=axs[0],
    legend_kwargs={
        'bbox_to_anchor':(1,1), 
        'loc':'upper left', 
        'ncols':1,
        'colors' : {'mut':'b', 'wt':'grey'}
    }
)
format_ax(axs[0], xlabel='UMAP1 (expr)', ylabel='UMAP2 (expr)', title='SRSF2 status')

# Sample 
sample_colors = create_palette(df_, 'sample', 'twilight')
draw_embeddings(
    df_, 
    cat='sample', 
    ax=axs[1],
    legend_kwargs={
        'bbox_to_anchor':(1,1), 
        'loc':'upper left', 
        'ncols':1,
        'colors' : sample_colors
    }
)
format_ax(axs[1], xlabel='UMAP1 (expr)', ylabel='UMAP2 (expr)', title='Sample')
fig.tight_layout()
plt.show()

##

# Fine-annotation
fig, ax = plt.subplots(figsize=(10.5, 5))
draw_embeddings(
    df_, 
    cat='lineage', 
    ax=ax,
    legend_kwargs={'bbox_to_anchor':(1,1), 'loc':'upper left', 'ncols':2}
)
format_ax(ax, xlabel='UMAP1 (expr)', ylabel='UMAP2 (expr)', title='Refined cell type')
ax.text(.72, .05, f'n cells: {adata_all.shape[0]}', transform=ax.transAxes)
fig.tight_layout()
plt.show()


##


# sAML1-only: short Cellula workflow
cells_ = adata_all.obs.query('sample == @sample').index
sAML1 = anndata.AnnData(
    X=adata_all[cells_,:].layers['counts'],
    obs=adata_all[cells_,:].obs.loc[:, 
        [   
            'nCount_RNA',  'nFeature_RNA', 'percent_MT',
            'SRSF2_genotype', 'aggregated_lineage',
            'STAG2', 'FRYL', 'IDH1', 'SRSF2', 'EZH2', 
            'RUNX1', 'SETD2', 'FLT3', 'HDAC7', 'CSF3R', 'LSC_Ng1',
            'pt.Myelocytes', 
        ]
    ],
    var=adata_all.var,
    layers={'raw':adata_all[cells_,:].layers['counts']}
)
sAML1.obs = sAML1.obs.rename(columns={'nCount_RNA':'nUMIs', 'percent_MT':'mito_perc'})

# Create layers, PCA, kNN and embeddings
_, red = standard_pp_recipe(sAML1)
red.obsm['scaled|original|X_pca'] = pca(red)['X_pca']
compute_kNN(red, layer='scaled', k=15)                     # auto n_npcs
embs = embeddings(red, paga_groups='aggregated_lineage')
sAML1.obsm['X_pca'] = red.obsm['scaled|original|X_pca']
sAML1.obsm['X_umap'] = embs.loc[:, ['UMAP1', 'UMAP2']].values

# Format meta

# Add mito-clones labels to adata, take only relevant fields
sAML1.obs = sAML1.obs.join(mito_clones_df.rename(columns={'0':'MT_clone'}))
sAML1 = sAML1[~sAML1.obs['MT_clone'].isna(), :].copy()
meta = sAML1.obs.loc[:,
    [
        'SRSF2_genotype', 'LSC_Ng1', 'aggregated_lineage', 'MT_clone',
        'STAG2', 'FRYL', 'IDH1', 'SRSF2', 'EZH2', 'RUNX1', 'SETD2', 
        'FLT3', 'HDAC7', 'CSF3R', 
    ]
]


##


# UMAP (expression) clones MI_TO, lineages, nuclear mutations
fig, axs = plt.subplots(1,3, figsize=(14,4))

# Lineage
draw_embeddings(
    embs.loc[meta.index],
    cat='aggregated_lineage', 
    ax=axs[0],
    legend_kwargs={
        'bbox_to_anchor':(1,1), 
        'loc':'upper left',
        'label' : 'Cell type'
    }
)
format_ax(
    axs[0], xlabel='UMAP1 (expr)', ylabel='UMAP2 (expr)', 
    title='Haematopoietic lineage'
)

# MT-clones expr
df_ = pd.crosstab(meta['MT_clone'], meta['aggregated_lineage'], normalize=0).T
clone_order = df_.loc[
    ['HSC_progenitors', 'Early_myeloid', 'Late_myeloid'],
:].sum(axis=0).sort_values(ascending=False).index

clones_colors = create_palette(meta, 'MT_clone', 'rainbow')
clones_colors['unassigned'] = 'lightgrey'

draw_embeddings(
    meta.join(embs.loc[:, ['UMAP1', 'UMAP2']]),
    cat='MT_clone', 
    ax=axs[1],
    legend_kwargs={
        'bbox_to_anchor':(1,1), 'loc':'upper left',
        'colors' : clones_colors,
        'label' : 'MT-clone'
    }
)
format_ax(
    axs[1], xlabel='UMAP1 (expr)', ylabel='UMAP2 (expr)', 
    title='MT-clones'
)

# MT-clones SNVs
draw_embeddings(
    mito_umap.loc[meta.index].join(meta),
    cat='MT_clone', 
    ax=axs[2],
    legend_kwargs={
        'bbox_to_anchor':(1,1), 'loc':'upper left',
        'colors' : clones_colors
    }
)
format_ax(
    axs[2], xlabel='UMAP1 (MT-SNVs)', ylabel='UMAP2 (MT-SNVs)', 
    title='MT-clones'
)

fig.tight_layout()
plt.show()


##


# Bars clone MI-TO and lineage

# Fig
fig = plt.figure(figsize=(10, 6))
gs = GridSpec(1, 2, width_ratios=[10, 7])

# Clone MI_TO
ax = fig.add_subplot(gs[0,0])
df_ = (
    meta.groupby('MT_clone')
    .size()
    .to_frame('n_cells')
    .sort_values('n_cells', ascending=False)
)
bar(df_, 'n_cells', c='#E9E7E7', edgecolor='k', s=0.5, ax=ax)
format_ax(ax, xticks=df_.index, ylabel='n cells', title='MITO clones', rotx=90)
ax.spines[['left', 'right', 'top']].set_visible(False)

##

# lineage
ax = fig.add_subplot(gs[0,1])
df_ = (
    meta.groupby('aggregated_lineage')
    .size()
    .to_frame('n_cells')
    .sort_values('n_cells', ascending=False)
)
bar(df_, 'n_cells', c='#E9E7E7', edgecolor='k', s=0.5, ax=ax)
format_ax(ax, xticks=df_.index, ylabel='n cells', title='Haematopoietic lineage', rotx=90)
ax.spines[['left', 'right', 'top']].set_visible(False)


# Show and save
fig.tight_layout()
plt.show()


##


# Heatmap contingency
fig, ax = plt.subplots(figsize=(5.5, 5))

df_ = pd.crosstab(meta['MT_clone'], meta['aggregated_lineage'], normalize=0).T
clone_order = df_.loc[
    ['HSC_progenitors', 'Early_myeloid', 'Late_myeloid'],
:].sum(axis=0).sort_values(ascending=False).index

plot_heatmap(df_.loc[:, clone_order], ax=ax, 
    title='Haematopoietic lineages x MT-clones', annot=True, label='MT-clone fraction'
)

fig.tight_layout()
plt.show()



##


# Nuclear mutations 
df_ = (
    pd.concat(
        [meta.loc[:, ['aggregated_lineage']], meta.iloc[:, 4:]], 
        axis=1
    )
    .melt(id_vars='aggregated_lineage', var_name='gene', value_name='mut_status')
)
mut_colors = {'NA':'#E9E7E7', 'unknown':'grey', 'mutated':'r', 'wild-type':'b'}

# Viz
fig, axs = plt.subplots(2, 5, figsize=(12, 6), sharex=True, sharey=True)

i = 0
j = 0
for x in df_['gene'].unique(): 
    
    bb_plot(
        df_.query('gene==@x'), 
        cov1='aggregated_lineage', cov2='mut_status', 
        colors=mut_colors,
        legend=True if (i==1 and j==4) else False, 
        ax=axs[i,j]
    )
    format_ax(axs[i,j], title=x, xlabel='Abundance %')
    axs[i,j].spines[['left', 'right', 'top']].set_visible(False)
    
    print(x, i, j)
    j += 1
    if j == 5:
        i += 1
        j = 0
        

# Show and save
fig.tight_layout()
plt.show()


##


# All of them, mutated fraction
fig, ax = plt.subplots(figsize=(5.5, 5))

genes = meta.iloc[:, 4:].columns[~meta.iloc[:, 4:].columns.isin(['STAG2', 'FRYL', 'EZH2'])]
L = []
for gene in genes:
    print(gene)
    d = pd.crosstab(meta['MT_clone'], meta[gene], normalize=0)
    L.append(d.assign(gene=gene))




df_mut = (
    pd.concat(L, axis=0).loc[:, ['mutated', 'gene']]
    .reset_index()
    .pivot_table(index='MT_clone', values='mutated', columns='gene')
)

clone_order = df_mut.mean(axis=1).sort_values(ascending=False).index 
    
plot_heatmap(
    df_mut.loc[clone_order], 
    palette='mako', annot_size=7,
    ax=ax, title='Nuclear SNVs mutated fraction', annot=True, label='MT-clone fraction',
    ylabel='MT-clone'
)

fig.tight_layout()
plt.show()


##
