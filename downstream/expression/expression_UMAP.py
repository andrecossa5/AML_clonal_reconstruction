"""
Expression data visualization.
"""

import os
from plotting_utils._utils import *
from plotting_utils._plotting_base import *
from mito_utils.embeddings_plots import *
from matplotlib.gridspec import GridSpec
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'expression')


##


# Read data
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)

# Explore
meta.columns
meta.describe().T
meta.groupby('sample_id').size()
lineage = [
    'lineage', 'lineage_similarity_score', 'pt.Myelocytes', 'pt.cDC',
    'pt.Megakaryocte', 'pt.Erythroid', 'pt.B.cells', 'aggregated_ct',
    'singleR_label', 'singleR_score'
]
meta['aggregated_ct'].value_counts()
meta['lineage'].value_counts()
meta['singleR_label'].value_counts()
meta['condition'].value_counts()
meta['cohort'].value_counts()
meta['malignant_class_occupancy'].value_counts()
meta.columns
meta['SRSF2'].value_counts()

# Colors
colors = {
    'cell_type' : {
        'HSC_MPP' : sns.color_palette('Reds')[1],
        'Early_myeloid' : sns.color_palette('Reds')[-2],
        'Late_myeloid' : sns.color_palette('Reds')[-1],
        'T_CD4' : sns.color_palette('Greens')[1], 
        'T_gd' : sns.color_palette('Greens')[-2],
        'T_CD8' : sns.color_palette('Greens')[-1],
        'B_early' : sns.color_palette('Blues')[-3],
        'B_mature' : sns.color_palette('Blues')[-1],
        'Mono' : sns.color_palette('Oranges')[-3],
        'Eo_baso_mast' : sns.color_palette('Oranges')[-1],
        'Stromal' : 'lightgrey',
        'Erythroid' : '#7F00FF',
        'NK_cells' : 'k',
        'DC' : 'yellow',
        'MK' : '#703A10'
    },
    'malignant_class' : {'malignant':'#cc160c', 'tme':'#240ee6', 'hBM':'#c9c5bd'},
    'cohort' : {'SRSF2mut':'#e66b0e', 'SRSF2wt':'#0e73e6', 'hBM':'#c9c5bd'},
    'samples' : create_palette(meta, 'sample_id', sc.pl.palettes.vega_10)
}

# Save colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'wb') as f:
    pickle.dump(colors, f)


##


# Cell type composition across samples
fig = plt.figure(figsize=(7, 3))
gs = GridSpec(1, 2, width_ratios=[10, 2])

ax = fig.add_subplot(gs[0,0])
var = 'aggregated_ct'
bb_plot(
    meta, 'sample_id', var, legend=False,
    ax=ax, colors=colors['cell_type'],
)
ax.set(title=None)
ax.spines[['right', 'top', 'left', 'top']].set_visible(False)

ax = fig.add_subplot(gs[0,1])
counts = meta['sample_id'].value_counts()[np.unique(meta['sample_id'])]
index = range(counts.size)
ax.barh(index, counts, color='k', alpha=.6, height=0.95)
ax.spines[['left', 'right', 'top']].set_visible(False)
ax.set_yticks([])
ax.set(xlabel='n cells')

# Save
fig.tight_layout(h_pad=.9)
fig.savefig(os.path.join(path_results, 'cell_composition.pdf'), dpi=300)

# Only legend
# fig, ax = plt.subplots(figsize=(6,4))
# add_legend(ax=ax, colors=colors['cell_type'], label='Cell type', loc='center', bbox_to_anchor=(.5,.5), ncols=4, artists_size=8, label_size=8, ticks_size=6)
# fig.tight_layout()
# ax.axis('off')
# fig.savefig(os.path.join(path_results, 'legend.pdf'), dpi=300)


##


# Load umap, join and plot
umap = pd.read_csv(os.path.join(path_data, 'expression', 'umap_coord.csv'), index_col=0)
umap.columns = ['UMAP1', 'UMAP2']
df_ = umap.join(meta)


# General annot: sample_id, cohort, malignant status
fig, axs = plt.subplots(1,3, figsize=(9.5,3))

# sample_id
draw_embeddings(
    df_, 
    cat='sample_id', 
    ax=axs[0],
    axes_kwargs={'legend':False}
)
axs[0].axis('off')

# cohort
draw_embeddings(
    df_, 
    cat='cohort', 
    ax=axs[1],
    legend_kwargs={
        'bbox_to_anchor':(1.2,1), 
        'loc':'upper right', 
        'ncols':1,
        'label':'',
        'colors' : colors['cohort']
    }
)
axs[1].axis('off')

# cohort
draw_embeddings(
    df_, 
    cat='malignant_class_occupancy', 
    ax=axs[2],
    legend_kwargs={
        'bbox_to_anchor':(1.2,1), 
        'loc':'upper right', 
        'ncols':1,
        'label':'',
        'colors': colors['malignant_class']
    }
)
axs[2].axis('off')

# Save
fig.subplots_adjust(top=.9, left=.1, right=.9, bottom=.1)
fig.savefig(os.path.join(path_results, 'general_umap.pdf'), dpi=500)


##


# Fine-annotation
fig, ax = plt.subplots(figsize=(6, 5))
draw_embeddings(
    df_, 
    cat='aggregated_ct', 
    ax=ax,
    legend_kwargs={
        'bbox_to_anchor':(1,.5), 
        'loc':'center left', 
        'ncols':1,
        'label':'',
        'colors':colors['cell_type']
    }
)
ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cell_type_umap.pdf'), dpi=500)


##


