"""
Create a consensus set of characters.
"""

import os
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from mito_utils.preprocessing import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
from mito_utils.plotting_base import *
from mito_utils.utils import ji
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')


##


# Read data
sample = 'sAML1'
path_results = os.path.join(path_results, sample)

# AFM
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
meta = meta.query('sample_id==@sample')
afm.obs = (
    afm.obs
    .join(
    meta[['nCount_RNA', 'percent_MT', 
          'malignant_class_occupancy', 'aggregated_ct']])
)

# Jobs and dataset stats
jobs_df = pd.read_csv(os.path.join(path_results, f'{sample}_job.csv'), index_col=0).set_index('job_id')
dataset_df = pd.read_csv(os.path.join(path_results, f'{sample}_dataset.csv'), index_col=0)
vars_df = pd.read_csv(os.path.join(path_results, f'{sample}_vars.csv'), index_col=0)

df = (
    dataset_df  # .query('metric in @cols')
    .pivot_table(values='value', columns='metric', index='job_id')
    .join(jobs_df)
    .join(vars_df.groupby('job_id')[['median_af_in_positives', 'prior']].median())
)

# Redundancy among sets
jobs = vars_df['job_id'].unique()
J = np.zeros((jobs.size, jobs.size))
for i,x in enumerate(jobs):
    var_x = set(vars_df.query('job_id==@x').index)
    for j,y in enumerate(jobs):
        var_y= set(vars_df.query('job_id==@y').index)
        J[i,j] = len(var_x & var_y) / len(var_x | var_y)


L = linkage(J, method='weighted')
order = leaves_list(L)

fig, ax = plt.subplots(figsize=(5,5))
ax.imshow(J[np.ix_(order, order)])
fig.tight_layout()
plt.show()


##


# Annot:

# # Group 1:
# # "Ultra high sensitivity" call attempted (first 18 sets): 
# # params --> 'min_n_confidently_detected<10 and af_confident_detection==.01'
# # sets cluster per min_n_confidently_detected (5,2) and then per min_site_cov (20,25,50): 
# # the 3 values  min_frac_negative (.2,.75,.9) stay togehter (i.e., very high overlaps)
# group_1 = 'min_n_confidently_detected<10 and af_confident_detection==.01'
# df.query(group_1).T
# jobs_group_1 = df.query(group_1).index
# variants_group_1 = vars_df.query('job_id in @jobs_group_1').index.unique()
# variants_group_1.size
# 
# # Group 2:
# # "Ultra high sensitivity" call attempted (sets 63-72): 
# # params --> 'min_n_confidently_detected==10 and af_confident_detection==.01'
# # sets cluster per min_site_cov (20,25,50) 
# # the 3 values  min_frac_negative (.2,.75,.9) cluster together (i.e., very high overlaps)
# group_2 = 'min_n_confidently_detected==10 and af_confident_detection==.01'
# df.query(group_2).describe().T
# jobs_group_2 = df.query(group_2).index
# variants_group_2 = vars_df.query('job_id in @jobs_group_1').index.unique()
# variants_group_2.size
#
## Group 3:
# "High sensitivity" and "very rare" call attempted. n variants 500-2k
# params --> 'min_n_confidently_detected==2 and af_confident_detection==.05'
# sets cluster per min_site_cov (20,25,50) 
# the 3 values  min_frac_negative (.2,.75,.9) cluster together (i.e., very high overlaps)


##


# Group 1, aka "Ultra high sensitivity" group. n variants ~1-9k
ultra_high_jobs = df.query('af_confident_detection==.01').index
ultra_high_vars = vars_df.query('job_id in @ultra_high_jobs').index.unique()
df.loc[ultra_high_jobs].describe().T

# Group 2, aka "High sensitivity and ultra-rare" group. Variants ~2k
high_ultra_rare_jobs = df.query('af_confident_detection==.05 and min_n_confidently_detected==2').index
high_ultra_rare_vars = vars_df.query('job_id in @high_ultra_rare_jobs').index.unique()
df.loc[high_ultra_rare_jobs].describe().T

# Group 3, aka "High sensitivity and rare" group. Variants ~2k
high_rare_job = df.query('af_confident_detection==.05 and min_n_confidently_detected>=3').index
high_rare_vars = vars_df.query('job_id in @high_rare_job').index.unique()
df.loc[high_rare_job].describe().T

# Group 4: aka "Medium sensitivity and ultra-rare" group. Variants 732
medium_ultra_rare_jobs = df.query('af_confident_detection==.1 and min_n_confidently_detected==2').index
medium_ultra_rare_vars = vars_df.query('job_id in @medium_ultra_rare_jobs').index.unique()
df.loc[medium_ultra_rare_jobs].describe().T

# VG Miller, 2022, most permissive treshold in supp: vars 227
miller_jobs = df.query('af_confident_detection>=.1 and min_n_confidently_detected>=5').index
miller_vars = vars_df.query('job_id in @miller_jobs').index.unique()
df.loc[miller_jobs].describe().T

# VG Miller, 2022, less permissive treshold in supp: vars NON CE L'ABBIAMO
# df['af_confident_detection'].value_counts()
# miller_jobs = df.query('af_confident_detection>=.5 and min_n_confidently_detected>=5').index
# miller_vars = vars_df.query('job_id in @miller_jobs').index.unique()
# df.loc[miller_jobs].describe().T


##


# Rank top jobs
# weights = [3, 1, 1, 1, 3, 3, 1]
# df_ranks = pd.concat([df.iloc[:,:-1].rank(ascending=False), df.iloc[:,[-1]].rank()], axis=1).astype(int)
# top_jobs = (df_ranks * weights).median(axis=1).rank().sort_values().index[:10]
# 
# # Assess their parameters
# dataset_df.query('job_id in @top_jobs').pivot_table(values='value', columns='metric', index='job_id')
# jobs_df.query('job_id in @top_jobs')
# 
# # Create consensus set
# consensus_df = vars_df.query('job_id in @top_jobs').iloc[:,:-2].drop_duplicates()
# consensus_vars = consensus_df.index.unique().to_list()


##


# Filter and assess last time

# Filter consensus vars heteroplasmy for kimura

# Miller, but n0>.9*ncells
min_neg_cells = round(.95*afm.shape[0])
vois = vars_df.query('job_id in @miller_jobs and n0>@min_neg_cells').index.unique()

dataset_df, a = filter_cells_and_vars(
    afm,  
    sample_name='sAML1',
    filtering='MI_TO',
    # variants='vois',
    tree_kwargs = { 'metric' : 'jaccard', 'solver' : 'UPMGA', 't' : .1},
    lineage_column='malignant_class_occupancy',
    spatial_metrics=True,
    fit_mixtures=True, 
    path_priors=path_priors
)
a.var.columns

# pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names).to_csv(os.path.join(path_results, 'prova_kimura.csv'))


# Filter cells with at least one muts > t (.1 AF)
# X_bin = np.where(a.X>=.1, 1, 0)
# a = a[a.obs_names[X_bin.sum(axis=1)>=1],:]

# VG clonal grouping and ordering
pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names).to_csv('afm.csv')
corr = np.corrcoef(a.X.T)
pd.DataFrame(corr).to_csv('corr.csv')

# See R

# Classic at multiple res here
fig, axs = plt.subplots(1,5, figsize=(10,4))

for i, t in enumerate(np.linspace(.02, .15, 5)):
    D = pair_d(a, t=t, metric='jaccard', weights=1-a.var['prior'])
    order = leaves_list(linkage(D))
    axs[i].imshow(D[np.ix_(order, order)], cmap='viridis_r')
    axs[i].set(xticks=[], yticks=[], title=f'AF: {t:.2f}')

fig.tight_layout()
plt.show()



# Tree
tree = build_tree(a, solver='spectral', t=.1, weights=1-a.var['prior'].values)
tree.collapse_mutationless_edges(True)

fig, ax = plt.subplots(figsize=(5,5))
colors = {'malignant':'r', 'tme':'b', np.nan : 'grey'}
plot_tree(
    tree, 
    categorical_cmap_annot=colors,
    meta=['malignant_class_occupancy'],
    colorstrip_width=6,
    ax=ax
)
fig.tight_layout()
plt.show()



from mito_utils.dimred import *
from mito_utils.kNN import *


X, _ = reduce_dimensions(a, 'PCA', metric='jaccard', n_comps=30)
u = UMAP(n_neighbors=5, metric='jaccard')
X = u.fit_transform(X)
plt.plot(X[:,0], X[:,1], 'ko', alpha=.2)
plt.show()

L = []
for _ in range(25):
    variants = a.var_names.to_series().sample(round(a.shape[1]*.9), replace=True)
    index, dists, conn = kNN_graph(a[:, variants].X, nn_kwargs={'metric':'jaccard'}, k=5)
    L.append(dists)

X = L[0].A
[ X.A for X in L ])

np.sum([L[0].A, L[1].A])

pd.Series(np.concatenate(l_)).value_counts(normalize=True)













##


# np.sum((a.var['deltaBIC']>0) & (a.var['enriched_tme']>0) & (a.var['enriched_malignant']>0))
t = .001
plt.plot(
    (a[a.obs['aggregated_ct']=='NK_cells'].X>=t).sum(axis=0)/(a.obs['aggregated_ct']=='NK_cells').sum(axis=0),
    (a[a.obs['aggregated_ct']=='T_CD4'].X>=t).sum(axis=0)/(a.obs['aggregated_ct']=='T_CD4').sum(axis=0),
    # (a[a.obs['malignant_class_occupancy']=='malignant'].X>=t).sum(axis=0)/(a.obs['malignant_class_occupancy']=='malignant').sum(axis=0),
    # np.median(a[a.obs['malignant_class_occupancy']=='malignant'].X, axis=0),
    'ko', alpha=.2
)
# plt.xlim((-.01,.09))
# plt.ylim((-.01,.09))
# plt.xlim((-.01,.09))
# plt.ylim((-.01,.09))
plt.show()


early_myeloid_malignant = a.obs.query('aggregated_ct in ["Early_myeloid", "HSC_MPP"] and malignant_class_occupancy == "malignant"').index
early_myeloid_tme = a.obs.query('aggregated_ct in ["Early_myeloid", "HSC_MPP"] and malignant_class_occupancy == "tme"').index

t = .05
plt.plot(
    (a[early_myeloid_malignant,:].X>=t).sum(axis=0)/early_myeloid_malignant.size,
    (a[early_myeloid_tme,:].X>=t).sum(axis=0)/early_myeloid_tme.size,
    # (a[a.obs['malignant_class_occupancy']=='malignant'].X>=t).sum(axis=0)/(a.obs['malignant_class_occupancy']=='malignant').sum(axis=0),
    # np.median(a[a.obs['malignant_class_occupancy']=='malignant'].X, axis=0),
    'ko', alpha=.2
)
plt.xlim((-.01,.3))
plt.ylim((-.01,.3))
plt.show()





# Kimura
# ...

# Read kimura results and apply last filters
sig_variants = pd.read_csv(os.path.join(path_results, 'results_kimura.csv'), index_col=0)['sig_variants'].unique()
final_vars_df = a.var.loc[a.var_names.isin(sig_variants)]

# Final filtering
final_variants = final_vars_df.loc[final_vars_df['enriched_malignant'] | final_vars_df['enriched_tme']].index
dataset_df, a = filter_cells_and_vars(
    afm,  
    variants=final_variants.to_list(),
    tree_kwargs = { 'metric' : 'jaccard', 'solver' : 'UPMGA', 't' : .1},
    lineage_column='malignant_class_occupancy',
    spatial_metrics=True,
    fit_mixtures=False, 
    path_priors=path_priors
)


##