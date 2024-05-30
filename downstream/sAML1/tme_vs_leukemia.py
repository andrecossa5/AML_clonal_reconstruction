"""
sAML1 MT-SNVs characterization
"""

import os
import pickle
from mito_utils.clustering import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
matplotlib.use('macOSX')



##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'sAML1')
path_var_selection = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')
# path_top_trees = os.path.join(path_results, 'top_trees')

# Read MT_data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = meta.query('sample_id==@sample')
afm.obs = afm.obs.join(meta)

# Filter leukemic cells with MT-variant
_, a = filter_cells_and_vars(
    afm, sample_name=sample, 
    path_priors=path_priors, 
    max_prior=.1, 
    filtering='MI_TO',
    filtering_kwargs={
        "min_n_confidently_detected" : 5,
        "af_confident_detection" : 0.05,
        "min_frac_negative" : 0
    },
    fit_mixtures=True,
    only_positive_deltaBIC=True,
    max_AD_counts=2,
    lineage_column='malignant_class_occupancy'
)

# Find mut
rank_clone_variants(a, 'malignant_class_occupancy', 'malignant', filter_vars=False)
mut = '15843_T>C'

# Add to afm.obs
a.obs[mut] = a[:,mut].X.toarray().flatten()

# Load umap, join and plot
umap = pd.read_csv(os.path.join(path_data, 'expression', 'umap_coord.csv'), index_col=0)
umap.columns = ['UMAP1', 'UMAP2']
umap = umap.loc[umap.index.str.contains('sAML1')]
umap.index = umap.index.map(lambda x: x.split('-')[0])

# Viz
df_ = a.obs.join(umap)

# Colors
with open(os.path.join(path_data, 'meta', 'colors.pickle'), 'rb') as f:
    colors = pickle.load(f)
colors_malignant = { k:v for k,v in colors['malignant_class'].items() if k != 'hBM'}
colors_ct = { k:v for k,v in colors['cell_type'].items() if k in df_['aggregated_ct'].unique() }

##


fig, axs =  plt.subplots(1,3,figsize=(13,4))

draw_embeddings(
    df_.loc[~df_['malignant_class_occupancy'].isna()], 
    cat='aggregated_ct', 
    ax=axs[0],
    title='',
    legend_kwargs={
        'bbox_to_anchor':(1.2,1), 
        'loc':'upper right', 
        'ncols':1,
        'label':'Cell type',
        'colors': colors_ct
    }
)
axs[0].axis('off')

draw_embeddings(
    df_.loc[~df_['malignant_class_occupancy'].isna()], 
    cat='malignant_class_occupancy', 
    ax=axs[1],
    title='',
    legend_kwargs={
        'bbox_to_anchor':(1.2,1), 
        'loc':'upper right', 
        'ncols':1,
        'label':'Malignat status',
        'colors': colors_malignant
    }
)
axs[1].axis('off')

draw_embeddings(
    df_.loc[~df_['malignant_class_occupancy'].isna()], 
    cont=mut, 
    ax=axs[2],
    cbar_kwargs={
        'palette' : 'viridis',
        'vmin': 0,
        'vmax': .5,
        'label_size' : 12, 
        'ticks_size' : 10,  
        'layout' : 'outside'
    }
)
axs[2].axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'malingnat_vs_tme.png'), dpi=500)


##