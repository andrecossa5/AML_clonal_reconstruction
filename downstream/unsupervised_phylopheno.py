"""
Viz unsupervised PATH output transition rates.
"""

import os
import pickle
from mito_utils.phylo_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


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

# Genes z-scores phylocorr visualization
samples = ['sAML1', 'AML2', 'AML3', 'AML5']


##


# Viz 
for sample in samples:
    fig, ax = plt.subplots(figsize=(4,4))
    df = pd.read_csv(
        os.path.join(path_top_trees, sample, f'ranked_genes_malignant_unsupervised.csv'),
        index_col=0
    )
    rank_plot(df, cov='Z', ax=ax, fig=fig, n_annotated=20)
    format_ax(ax, xlabel='HVG rank', ylabel='Phylogenetic correlation z-score', title=sample)
    fig.tight_layout()
    fig.savefig(os.path.join(path_top_trees, sample, f'malignant_genes.pdf'), dpi=500)


##


# GSEA viz
fig = plt.figure(figsize=(5,5))

for i, sample in enumerate(['sAML1']):#samples):

    ax = fig.add_subplot(1,1,i+1)
    # ax = fig.add_subplot(1,4,i+1)
    gsea = pd.read_csv(
        os.path.join(path_top_trees, sample, f'H_gsea_malignant_unsupervised.csv'),
        index_col=0
    )
    gsea = gsea.sort_values('NES', ascending=False)
    gsea['-log10(padj)'] = -np.log10(gsea['padj'])

    sns.scatterplot(data=gsea.head(20), x='NES', y='pathway', size='-log10(padj)', ax=ax, 
                    legend=True if i == 3 else False)
    format_ax(ax, ylabel='', xlabel='NES', reduced_spines=True, title=sample)
    ax.set(xlim=(.85, 1.75))

fig.tight_layout()
fig.savefig(os.path.join(path_top_trees, f'H_gsea_malignant_unsupervised.pdf'), dpi=500)


##