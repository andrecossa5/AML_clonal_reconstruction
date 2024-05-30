"""
Trial filter_cells_and_vars.
"""

import os
from mito_utils.preprocessing import *
from mito_utils.utils import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Args
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
sample = 'sAML1'

# Paths
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')


##


# Read
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
meta = meta.query('sample_id==@sample')
afm.obs = (
    afm.obs
    .join(
    meta[['nCount_RNA', 'percent_MT', 
          'malignant_class_occupancy', 'aggregated_ct']])
)


##


# MT-coverage / MT-genes (and not) expression
df_ = (
    afm.obs.join(
    afm.uns['per_position_coverage']
    .mean(axis=1).to_frame('mean_MT_coverage'))
    .assign(nCount_RNA_MT=lambda x: x['nCount_RNA']*x['percent_MT']/100)
)
fig, axs = plt.subplots(1,3, figsize=(12,4))
order = df_.groupby('aggregated_ct')['mean_MT_coverage'].mean().sort_values(ascending=False).index
box(df_, 'aggregated_ct', 'mean_MT_coverage', ax=axs[0], c='r', order=order)
format_ax(axs[0], rotx=90, ylabel='Mean MT-base coverage')
box(df_, 'aggregated_ct', 'nCount_RNA_MT', ax=axs[1], c='b', order=order)
format_ax(axs[1], rotx=90, ylabel='MT-RNA nUMI')
box(df_, 'aggregated_ct', 'nCount_RNA', ax=axs[2], c='g', order=order)
format_ax(axs[2], rotx=90, ylabel='Total-RNA nUMI')
fig.suptitle('MT-genome coverage and expression')
fig.tight_layout()
plt.show()


##


# Gene mean expression vs MAESTER mean base coverage: sAML1
expr_counts_path = os.path.join(path_data, 'expression/MT_genes_expression_sAML1.csv')
expr = pd.read_csv(expr_counts_path, index_col=0)
expr = expr.loc[afm.obs_names]
mean_expr = expr.mean(axis=0)

# Annotate MT-genome sites
all_mt_genes_positions = [
    ["MT-ND1", 3307, 4262], ["MT-ND2", 4470, 5511], ["MT-CO1", 5904, 7445],
    ["MT-CO2", 7586, 8269], ["MT-ATP8", 8366, 8572], ["MT-ATP6", 8527, 9207],
    ["MT-CO3", 9207, 9990], ["MT-ND3", 10059, 10404], ["MT-ND4L", 10470, 10766],
    ["MT-ND4", 10760, 12137], ["MT-ND5", 12337, 14148], ["MT-ND6", 14149, 14673],
    ["MT-CYB", 14747, 15887], ["MT-TF", 577, 647], ["MT-TV", 1602, 1670],
    ["MT-TL1", 3230, 3304], ["MT-TI", 4263, 4331], ["MT-TQ", 4329, 4400],
    ["MT-TM", 4402, 4469], ["MT-TW", 5512, 5579], ["MT-TA", 5587, 5655],
    ["MT-TN", 5657, 5729], ["MT-TC", 5761, 5826], ["MT-TY", 5826, 5891],
    ["MT-TS1", 7518, 7585], ["MT-TD", 7513, 7585], ["MT-TK", 8295, 8364],
    ["MT-TG", 9991, 10058], ["MT-TR", 10405, 10469], ["MT-TH", 12138, 12206],
    ["MT-TS2", 12207, 12265], ["MT-TL2", 12266, 12336], ["MT-TE", 14674, 14742],
    ["MT-TT", 15888, 15953], ["MT-TP", 15956, 16023], ["12S rRNA", 648, 1601],
    ["16S rRNA", 1671, 3229]
]
mt_genes_positions = [ x for x in all_mt_genes_positions if x[0] in expr.columns ]

sites = afm.uns['per_position_coverage'].columns
annot = {}
for x in sites:
    x = int(x)
    mapped = False
    for mt_gene, start, end in mt_genes_positions:
        if x>=start and x<=end:
            annot[str(x)] = mt_gene
            mapped = True
    if not mapped:
        annot[str(x)] = 'other'

mean_site_cov = afm.uns['per_position_coverage'].T.mean(axis=1).to_frame('cov')
mean_site_cov['gene'] = pd.Series(annot)
mean_site_cov = mean_site_cov.query('gene!="other"').groupby('gene')['cov'].mean()
mean_site_cov = mean_site_cov[mean_expr.index]


##


# Viz
fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot(mean_expr.values, mean_site_cov.values, 'ko')
sns.regplot(data=pd.DataFrame({'expr':mean_expr, 'cov':mean_site_cov}), 
            x='expr', y='cov', ax=ax, scatter=False)
format_ax(ax, title='MT-transcripts counts vs site coverage',
          xlabel='Mean expression (gene nUMI, 10x)', 
          ylabel='Mean site coverage (per-site nUMI, MAESTER)')
corr = np.corrcoef(mean_expr.values, mean_site_cov.values)[0,1]
ax.text(.05, .9, f'Correlation: {corr:.2f}', transform=ax.transAxes)
fig.tight_layout()
plt.show()


##












