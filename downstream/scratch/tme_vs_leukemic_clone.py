"""
Clones specific gene-expression patterns using HVGs.
"""

import os
from scipy.sparse import csr_matrix
from mito_utils.preprocessing import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.clustering import *
matplotlib.use('macOSX')



##


# Paths
path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'sAML1')
path_var_selection = os.path.join(path_main, 'results', 'var_selection')
path_priors = os.path.join(path_data, 'vars_df', 'priors.csv')
# path_top_trees = os.path.join(path_results, 'top_trees')

# Read expression data
sample = 'sAML1'
meta = pd.read_csv(os.path.join(path_data, 'meta', 'cells_meta.csv'), index_col=0)
signatures = pd.read_csv(os.path.join(path_data, 'expression', 'malignant_AUCell.csv'), index_col=0)
counts = pd.read_csv(os.path.join(path_data, 'expression', 'sAML1_malignant_lognorm_counts.csv'), index_col=0)
umap = pd.read_csv(os.path.join(path_data, 'expression', 'umap_coord.csv'), index_col=0)

# CBC harmonization
counts.index = counts.index.map(lambda x: x.split('-')[0])
umap.index = umap.index.map(lambda x: x.split('-')[0])

# Expression
adata = AnnData(
    X=csr_matrix(counts.values),            # Log-norm counts 
    obs=meta.loc[counts.index],
    var=pd.DataFrame(index=counts.columns)
)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:,adata.var['highly_variable']].copy()
sc.pp.scale(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)

# MT
afm = read_one_sample(os.path.join(path_data, 'samples'), sample)
meta = meta.query('sample_id==@sample')
afm.obs = (
    afm.obs
    .join(
    meta[['nCount_RNA', 'percent_MT', 
          'malignant_class_occupancy', 'aggregated_ct']])
)

# Find mut
rank_clone_variants(afm, 'malignant_class_occupancy', 'malignant')


# Subset common cells
cells = list(set(adata.obs_names) & set(afm.obs_names) & set(umap.index))
adata = adata[cells,:].copy()
afm = afm[cells,:].copy()

# Problem with index...
# umap.loc[cells]
# umap.loc[cells].join(adata.obs).columns

#


# Filter leukemic cells with MT-vatiant
_, a = filter_cells_and_vars(afm, sample_name=sample, path_priors=path_priors, filtering='baseline')
test_summary = (a.var['n5']>2) & (a.var['prior']<=.05) & (a.var['mean_cov']>=25) & (a.var['quality']>=30)

variants = a.var.loc[test_summary].index


pd.Series(np.where(a[:,variants].X>0.05,1,0).sum(axis=1)).describe()

d = pair_d(a[:,variants], t=0.01)
order = leaves_list(linkage(d))

plt.imshow(1-d[np.ix_(order, order)], cmap='viridis')
plt.imshow(1-d[np.ix_(order[650:750], order[650:750])], cmap='viridis')
plt.imshow(1-d[np.ix_(order[600:800], order[600:800])], cmap='viridis')
plt.colorbar()
plt.show()

a.obs.iloc[order[600:800],:]['aggregated_ct'].value_counts()
a.obs.iloc[order[10:100],:]['aggregated_ct'].value_counts()

a.obs

remove_cells = a.obs_names[order][650:750]

calc_median_AF_in_positives(a[remove_cells,:], variants) - calc_median_AF_in_positives(a[~a.obs_names.isin(remove_cells),:], variants)

pd.Series(calc_median_AF_in_positives(a, n5_5)).sort_values(ascending=False)[:-200]

n5_5_somatic = a[:, n5_5].var['prior']<=.01



sns.kdeplot(afm[:,'15843_T>C'].X.flatten())
plt.show()


fig, ax = plt.subplot()



a[:,'15843_T>C'].var



tree = build_tree(a, solver='NJ', metric='jaccard')
tree.collapse_mutationless_edges(True)

fig, axs = plt.subplots(1,2,figsize=(10,5))
plot_tree(tree, ax=axs[0], extend_branches=False, orient=90)
plot_tree(tree, ax=axs[1], extend_branches=False, orient='down')
plt.show()



# Top expanded clones

# Find expanded clones
cs.tl.compute_expansion_pvalues(tree, min_clade_size=(0.01 * tree.n_cell), min_depth=0)
probability_threshold = 0.05
expansions = {}
for node in tree.depth_first_traverse_nodes():
    if tree.get_attribute(node, "expansion_pvalue") <= probability_threshold:
        expansions[node] = tree.get_attribute(node, "expansion_pvalue")
L = []
for node in expansions:
    clade = get_clades(tree)[node]
    L.append([node, tree.cell_meta.index.isin(clade).sum()])
expanded_df = (
    pd.DataFrame(L, columns=['node', 'ncells'])
    .sort_values('ncells', ascending=False)
    .set_index('node')
    .join(pd.Series(expansions).to_frame('p'))
    .query('ncells<=1000')
)
expanded = expanded_df.index.to_list()