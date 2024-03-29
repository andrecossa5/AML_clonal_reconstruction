"""
New, redeem-like MI_TO filter, on streoids.
"""

import os
from mito_utils.preprocessing import *
from mito_utils.kNN import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
from mito_utils.diagnostic_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')
import warnings
warnings.simplefilter('ignore')


##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_colors = os.path.join(path_data, 'meta', 'colors.pickle')
path_results = os.path.join(path_main, 'results', 'figure_out_sAML1')



##


# Read AFM and meta
sample = 'sAML1'

plt.plot(afm.uns['per_position_coverage'].mean(axis=0).values, 'k-')
plt.show()

afm.uns['per_position_coverage'].mean(axis=0).values.mean()


# Read, format and retain good cells
afm = read_one_sample(path_data, sample, with_GBC=False)
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
meta =  meta.query('sample_id=="sAML1"')
meta.index = meta.index.map(lambda x: x.split('-')[0])
afm.obs = afm.obs.join(meta)
afm = afm[~afm.obs['malignant_class_occupancy'].isna(),:].copy()

# Calculate vars_df, as in in Weng et al., 2024, and Miller et al. 2022 before.

# vars.tib <- tibble(var = rownames(af.dm),
#                    mean_af = rowMeans(af.dm),
#                    mean_cov = rowMeans(assays(maegtk)[["coverage"]])[as.numeric(cutf(rownames(af.dm), d = "_"))],
#                    quality = qual.num)

mean_sites = afm.uns['per_position_coverage'].mean(axis=0)
vars_df = pd.DataFrame(
    np.column_stack([
        afm.X.mean(axis=0),
        afm.var_names.map(lambda x: mean_sites[x.split('_')[0]]),
        np.ma.getdata(np.ma.mean(np.ma.masked_less_equal(afm.layers['quality'], 0), axis=0))
    ]),
    columns=['mean_af', 'mean_cov', 'quality'],
    index=afm.var_names
)

# Calculate the number of cells that exceed VAF thresholds 0, 1, 5, 10, 50 as in Weng et al., 2024

# vars.tib <- vars.tib %>%
#     mutate(n0 = apply(af.dm, 1, function(x) sum(x == 0))) %>%
#     mutate(n1 = apply(af.dm, 1, function(x) sum(x > 1))) %>%
#     mutate(n5 = apply(af.dm, 1, function(x) sum(x > 5))) %>%
#     mutate(n10 = apply(af.dm, 1, function(x) sum(x > 10))) %>%
#     mutate(n50 = apply(af.dm, 1, function(x) sum(x > 50)))
# Variant_CellN<-apply(af.dm,1,function(x){length(which(x>0))})
# vars.tib<-cbind(vars.tib,Variant_CellN)

vars_df = (
    vars_df
    .assign(
        n0 = lambda x: np.sum(afm.X==0, axis=0),
        n1 = lambda x: np.sum(afm.X>.01, axis=0),
        n5 = lambda x: np.sum(afm.X>.05, axis=0),
        n10 = lambda x: np.sum(afm.X>.1, axis=0),
        n50 = lambda x: np.sum(afm.X>.5, axis=0),
        Variant_CellN = lambda x: np.sum(afm.X==0, axis=0),
    )
)

# Filter Wang et al., 2024

# vars_filter.tib <- vars.tib %>% filter(mean_cov > 5, quality >= 30, n0 > 0.9*ncol(af.dm),Variant_CellN>=2)
     
## Apply the same filter as in MAESTER
# IsInfo<-function(x){
# total<-length(x)
# if(length(which(x<10))/total>0.1 & length(which(x>50))>10){
#     return("Variable")
# }else{
#     return("Non")
# }
# }
# Variability<-apply(af.dm,1,IsInfo) %>% data.frame(Info=.)

# vars_filter.tib<-Tomerge_v2(vars_filter.tib,Variability) 

min_site_cov = 5
min_var_quality = 30
min_frac_negative = .5
min_n_positive = 2
vars_df = (
    vars_df.loc[
        (vars_df['mean_cov']>min_site_cov) & \
        (vars_df['quality']>min_var_quality) & \
        (vars_df['n0']>min_frac_negative) & \
        (vars_df['mean_cov']>min_site_cov) & \
        (vars_df['Variant_CellN']>min_n_positive) 
    ]
)

# Apply the same filter as in MAESTER
# IsInfo<-function(x){
# total<-length(x)
# if(length(which(x<10))/total>0.1 & length(which(x>50))>10){
#     return("Variable")
# }else{
#     return("Non")
# }
# }
# Variability<-apply(af.dm,1,IsInfo) %>% data.frame(Info=.)

low_af = .01
high_af = .1
min_frac_cells_below_low_af = .1
min_n_cells_above_high_af = 2

test_low_homoplasy = (afm.X<low_af).sum(axis=0)/afm.shape[0] > min_frac_cells_below_low_af
test_detection_evidence = (afm.X>high_af).sum(axis=0) > min_n_cells_above_high_af
test = test_low_homoplasy & test_detection_evidence
vars_df['miller2022_patients'] = pd.Series(np.where(test, 'Variable', 'Non'), index=afm.var_names).loc[vars_df.index]


##


# Filter now
detectability_t = .0001
nt = 2

# Subsets
detectable_vars = vars_df.query('mean_af>@detectability_t').index
n5_nt_vars = vars_df.query('n5>=@nt').index
n1_nt_vars = vars_df.query('n1>=@nt').index
miller2022_patients = vars_df.query('miller2022_patients=="Variable"').index


## 


# Examine one subset
_, a = filter_cells_and_vars(afm, variants=miller2022_patients)
X_bin = np.where(a.X>=high_af,1,0)


##


# 1. n cells per var and n vars per cell
vars_df.loc[miller2022_patients, 'Variant_CellN'].describe()
pd.Series(np.sum(X_bin, axis=1)).describe()

# 2. Connectedness
D = pairwise_distances(X_bin, metric=lambda x, y: np.sum(np.logical_and(x, y)))
cell_conn = np.ma.masked_equal(D, np.diag(D)).mean(axis=1).data
pd.Series(cell_conn).describe()

# 3. Sparseness and n genotypes occurrence
np.sum(a.X==0) / np.product(a.shape)
d = AFM_to_seqs(a)
pd.Series(d).value_counts().values

# 4. Annot
miller2022_patients.map(lambda x: x.split('_')[1]).value_counts()# .sum()

# 4. Drop-out simulation



# 5. TME-leukemic clone?
tme = a.obs['malignant_class_occupancy']=='tme'
tme_prevalences = np.sum(X_bin[tme,:], axis=0)/tme.sum()
leukemic = a.obs['malignant_class_occupancy']=='malignant'
leukemic_prevalences = np.sum(X_bin[leukemic,:], axis=0)/leukemic.sum()

np.sum((tme_prevalences<=.1) & (leukemic_prevalences>=.2))
corr = np.corrcoef(leukemic_prevalences, tme_prevalences)[0,0]

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot(tme_prevalences, leukemic_prevalences, 'ko', markersize=3)
sns.regplot(x=tme_prevalences, y=leukemic_prevalences, ax=ax, scatter=False)
format_ax(ax, title=f'TME-malignant correlation: r2={corr:.2f}', 
          xlabel='TME prevalence', ylabel='Malignant prevalence', reduced_spines=True)
ax.set_xlim((-0.05,1.05))
ax.set_ylim((-0.05,1.05))
fig.tight_layout()
plt.show()

# 6. Tree mutations support

# Build tree
tree = build_tree(a, t=high_af, solver='NJ')

# Post-process tree
tree_collapsed = tree.copy()
tree_collapsed.collapse_mutationless_edges(True)
len(tree_collapsed.internal_nodes) / len(tree.internal_nodes)

# Get each internal node mutational status
def get_internal_node_muts(tree, internal_node):
    muts = tree.character_matrix.columns
    node_state = np.array(tree.get_character_states(internal_node))
    idx = np.where(node_state==1)[0]
    muts = muts[idx].to_list()
    return muts


##


def assess_internal_node_muts(a, clades, c, high_af=.01):
    """
    Assess the prevalence of all MT-SNVs assigned to a single internal node 
    within its clade (p1), and outside of it (p0). Return useful stats.
    """
    muts, cells = clades[c]
    other_cells = a.obs_names[~a.obs_names.isin(cells)]
    p1 = np.where(a[cells, muts].X>=high_af,1,0).sum(axis=0) / len(cells)
    p0 = np.where(a[other_cells, muts].X>=high_af,1,0).sum(axis=0) / len(other_cells)
    af1 = np.median(a[cells, muts].X, axis=0)
    af0 = np.median(a[other_cells, muts].X, axis=0)
    top_mut = (
        pd.DataFrame(
            {'p1':p1,'p0':p0, 'median_af0':af0, 'median_af1':af1, 
             'clade':[c]*len(muts), 'ncells':[len(cells)]*len(muts)}, 
            index=muts
        )
        .assign(p_ratio=lambda x: x['p1']/x['p0']+.0001)
        .assign(af_ratio=lambda x: x['median_af1']/x['median_af0']+.0001)
        .sort_values('p_ratio', ascending=False)
    )
    return top_mut


##


# Get clade muts
clades = { 
    x : (get_internal_node_muts(tree_collapsed, x), tree_collapsed.leaves_in_subtree(x)) \
    for x in tree_collapsed.internal_nodes if x != 'root'
}

# Quantify the prevalence ratio (clade vs rest mutation prevalence) of mutations assigned 
# to each internal node. Get the most supported value and mutation
stats = []
for c in clades:
    top_mut = assess_internal_node_muts(a, clades, c, high_af=high_af)
    s = top_mut.iloc[0,:]
    stats.append(s)

df_stats = pd.concat(stats, axis=1).T.reset_index().rename(columns={'index':'mut'}).set_index('clade')
final_muts = df_stats['mut'].unique()


##


# Rebuild with filtered muts

# Build tree
_, a = filter_cells_and_vars(afm, variants=final_muts)
X_bin = np.where(a.X>=high_af,1,0)

# Tree
tree = build_tree(a, t=high_af, solver='NJ')
tree_collapsed = tree.copy()
tree_collapsed.collapse_mutationless_edges(True)
len(tree_collapsed.internal_nodes) / len(tree.internal_nodes)

# Format into df
clades = { 
    x : (get_internal_node_muts(tree_collapsed, x), tree_collapsed.leaves_in_subtree(x)) \
    for x in tree_collapsed.internal_nodes if x != 'root'
}
stats = []
for c in clades:
    top_mut = assess_internal_node_muts(a, clades, c, high_af=high_af)
    s = top_mut.iloc[0,:]
    stats.append(s)

# Final supporting muts stats
df_stats = pd.concat(stats, axis=1).T.reset_index().rename(columns={'index':'mut'}).set_index('clade')
final_muts = df_stats['mut'].unique()
final_muts.size

df_stats.sort_values('ncells', ascending=False)


##


# Viz individual muts supporting some clades
# fig, ax = plt.subplots(figsize=(4,4.5))
# sns.kdeplot(a[cells, top_mut].X.flatten(), ax=ax, fill=True)
# sns.kdeplot(a[other_cells, top_mut].X.flatten(), ax=ax, fill=True)
# median_clade = np.median(a[cells, top_mut].X.flatten())
# median_rest = np.median(a[other_cells, top_mut].X.flatten())
# format_ax(
#     ax, title=f'{top_mut}: \n Median AF clade={median_clade:.2f} \n median AF rest={median_rest:.2f}', 
#     yticks=[], xlabel='AF'
# )
# for x in [median_clade, median_rest]:
#     ax.axvline(x, color='grey', linestyle='--')
# ax.spines[['right', 'left', 'top']].set_visible(False)
# fig.tight_layout()
# plt.show()


##


# Viz tree
fig, ax = plt.subplots(figsize=(7,7))
plot_tree(
    tree_collapsed, 
    ax=ax, orient=90, extend_branches=True,
    leaf_kwargs={'markersize':5},
    internal_node_kwargs={'markersize':0, 'markeredgecolor':None},
    cov_leaves='malignant_class_occupancy', cmap_leaves={'malignant':'r', 'tme':'b'}
)
fig.tight_layout()
plt.show()      


##


























































# # Filter cellsls
# afm = filter_cells_coverage(afm, mean_coverage=50)
# filter_baseline(afm).var_names
# 
# afm.uns['per_position_quality'].isna().values.all()
# 
# 
# 
# # Shallow filter ad hoc:
# test_sites = pd.Series(
#     np.mean(afm.uns['per_position_coverage'], axis=0) > 10,
#     index=afm.uns['per_position_coverage'].columns
# )
# sites = test_sites[test_sites].index
# test_vars_site_coverage = (
#     afm.var_names
#     .map(lambda x: x.split('_')[0] in sites)
#     .to_numpy(dtype=bool)
# )
# test_vars_site_coverage.sum()
# 
# # Test 2-4: 
# # variants seen in at least 2 cells;
# # variants with AF>0.01 in at least 2 cells;
# test_vars_quality = np.nanmean(afm.layers['quality'], axis=0) > 30
# 
# np.nanmin(afm.layers['quality'])
# 
# afm.layers['quality'][:,2000:2020]
# 
# test_vars_coverage = np.sum(afm.layers['coverage']>0, axis=0) > 2
# test_vars_AF = np.sum(afm.X > 0.01, axis=0) > 2
# 
# # Filter vars and sites
# test_vars = test_vars_site_coverage & test_vars_quality & test_vars_coverage & test_vars_AF 
# filtered = afm[:, test_vars].copy()
# filtered = remove_excluded_sites(filtered)