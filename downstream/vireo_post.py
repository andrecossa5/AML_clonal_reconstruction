"""
Visualize MQuad-vireoSNP clonal analysis output.
"""

import re
import cassiopeia as cs
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.distances import *
from mito_utils.clustering import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from mito_utils.embeddings_plots import *
import matplotlib
# matplotlib.use('macOSX')


##

# Args
# path_main = '/Users/IEO5505/Desktop/AML_clonal_reconstruction'
path_main = sys.argv[1]
sample = sys.argv[2]
best = sys.argv[3]

# Paths 
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'clustering', sample)
path_variants = os.path.join(path_main, 'results', 'variants_selection', sample)


##


def main():

    # Data

    # AFM
    afm = read_one_sample(path_data, sample=sample, with_GBC=False)
    _path = os.path.join(path_variants, 'MQuad_vars.csv')
    vois = pd.read_csv(_path, index_col=0)['0']
    a_cells, a = filter_cells_and_vars(
        afm,
        sample=sample,
        min_cov_treshold=50,
        variants=vois
    )

    # Cells meta (all fab4) and embeddings (all fab4)
    meta = pd.read_csv(
        os.path.join(path_data, 'meta', 'fab4_meta.csv'), 
        index_col=0, low_memory=False
    )
    meta = meta.query('sample == @sample')
    expr_umap = pd.read_csv(os.path.join(path_data, 'meta', 'umap.csv'), index_col=0)

    # vireoSNP output
    clones = pd.read_csv(os.path.join(path_results, 'vireo_clones.csv'), index_col=0)
    mt_umap = pd.read_csv(os.path.join(path_results, 'MT_SNVs_umap.csv'), index_col=0)

    # Join
    df = meta.loc[clones.index].join(clones)
    df_filtered = df.loc[lambda x: x[best] != 'unassigned']


    ##


    ############################## N clones and their cells

    fig, axs = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)

    df_ = (
        clones.loc[lambda x: x[best] != 'unassigned']
        .groupby(best).size()
        .sort_values(ascending=False)
        .to_frame('n')
        .assign(freq=lambda x: round(x['n'] / x['n'].sum(), 2))
    )

    bar(df_, 'n', l=.8, s=.8, c='#DBD6D6', edgecolor='k', ax=axs[0])
    format_ax(
        axs[0], title='n cells', ylabel='n cells', 
        xticks=df_.index.astype('str'), reduced_spines=True
    )

    bar(df_, 'freq', l=.8, s=.8, c='#DBD6D6', edgecolor='k', ax=axs[1])
    format_ax(
        axs[1], title='Prevalence', ylabel='Prevalence', 
        xticks=df_.index.astype('str'), reduced_spines=True
    )

    fig.suptitle(f'{sample} MT-clones')
    fig.savefig(
        os.path.join(path_results, 'clone_frequencies.png'),
        dpi=300
    )


    ##


    ############################## Viz contingency tables

    fig, axs = plt.subplots(1,2,figsize=(10, 5.5), constrained_layout=True)

    # Lineage
    df_ = (
        pd.crosstab(df_filtered[best], df_filtered['aggregated_lineage'], normalize=1)
    )
    plot_heatmap(df_, ax=axs[0], annot=True, label='Lineage fraction', rank_diagonal=True)

    # MT clone
    df_ = (
        pd.crosstab(df_filtered[best], df_filtered['aggregated_lineage'], normalize=0)
    )
    plot_heatmap(df_, ax=axs[1], annot=True, label='MT-clone fraction', rank_diagonal=True)

    fig.suptitle(f'{sample} MT clone-lineage relationship')
    fig.savefig(
        os.path.join(path_results, 'clone_lineages_fractions.png'),
        dpi=300
    )

    ##

    # n cells
    fig, ax = plt.subplots(figsize=(4.5, 5), constrained_layout=True)
    df_ = (
        pd.crosstab(df_filtered[best], df_filtered['aggregated_lineage'])
    )
    plot_heatmap(df_, ax=ax, annot=True, label='n cells', rank_diagonal=True)
    fig.suptitle(f'{sample} MT clone-lineage relationship')
    fig.savefig(
        os.path.join(path_results, 'clone_lineages_ncells.png'),
        dpi=300
    )


    ##


    ############################## Viz clones enrichment in aggregated lineages

    # Lineage-bias
    fig = plt.figure(figsize=(8,7))

    for i, x in enumerate(df_filtered['aggregated_lineage'].unique()):

        ax = fig.add_subplot(3,2,i+1)

        df_ = (
            compute_clonal_fate_bias(df_filtered, 'aggregated_lineage', best, x)
            .reset_index()
            .rename(columns={'index':'MT_clone'})
        )
        df_['MT_clone'] = df_['MT_clone'].astype('str')
        scatter(df_, x='MT_clone', y='fate_bias', c='k', s=100, ax=ax)
        ax.hlines(y=-np.log10(.05), xmin=0, xmax=df_.shape[0]-1, linestyles='dashed', colors='r')
        format_ax(
            ax, title=f'{x} fate bias', xlabel='MT clone', 
            ylabel='-log10(FDR)', reduced_spines=True
        )

    fig.tight_layout()
    fig.savefig(
        os.path.join(path_results, 'clones_fate_bias.png'),
        dpi=300
    )


    ##


    # Cumulative clone percentages
    cross = pd.crosstab(df_filtered[best], df_filtered['lineage'], normalize=0)
    res = compute_clonal_fate_bias(df_filtered, 'aggregated_lineage', best, 'HSC_progenitors')
    blast_like_clones = res.query('FDR<=.05').index.astype('str')

    fig, ax = plt.subplots(figsize=(6, 5))

    for clone in cross.index:
        is_blast_like = clone in blast_like_clones
        x = np.sort(cross.loc[clone, :])[::-1].cumsum()
        f = 'r.-' if is_blast_like else 'k.-'
        ax.plot(x, f)

    format_ax(
        ax, title='MT clones cumulative fractions', xticks='', reduced_spines=True,
        xlabel='All cell types, ranked', ylabel='MT clone fraction'
    )
    add_legend(
        label='MT clone', colors={'Blast-like':'r', 'Others':'k'},
        loc='lower right', bbox_to_anchor=(1,0), ax=ax, ticks_size=8, label_size=10
    )
    fig.tight_layout()
    fig.savefig(
        os.path.join(path_results, 'cumulative_clone_percentages.png'),
        dpi=300
    )


    ##


    ############################## Viz clones pseudotime and derailment_score ordering 

    # Viz
    fig = plt.figure(figsize=(10,7))

    dpts = [
        'pt.Myelocytes', 'pt.cDC', 'pt.Megakaryocte', 
        'pt.Erythroid', 'pt.B.cells', 'derailment_score'
    ]

    for i, x in enumerate(dpts):

        ax = fig.add_subplot(3,2,i+1)
        order = (
            df_filtered.groupby(best)
            [x].agg('median')
            .sort_values(ascending=False)
            .index
        )
        box(df_filtered, x=best, y=x, c='white', ax=ax, order=order)
        strip(df_filtered, x=best, y=x, c='k', s=1.5, ax=ax, order=order)
        name = x.split(".")[1] if not bool(re.findall('score', x)) else x
        format_ax(ax, title=f'DPT {name}', reduced_spines=True)

    fig.tight_layout()
    fig.savefig(
        os.path.join(path_results, 'clones_pseudotimes.png'),
        dpi=300
    )


    ##


    ############################## UMAPs 

    # Viz
    df_ = expr_umap.join(df_filtered, how='right')
    res = compute_clonal_fate_bias(df_, 'aggregated_lineage', best, 'HSC_progenitors')
    blast_like_clones = res.query('FDR<=.05').index.astype('str')
    df_['MT_clone_status'] = np.where(df_[best].isin(blast_like_clones), 'Blast-like', 'Others')

    # Viz states
    fig, axs = plt.subplots(1,3,figsize=(15,4.5))

    draw_embeddings(df_, cat='aggregated_lineage', ax=axs[0],
                    legend_kwargs={'loc':'upper left', 'bbox_to_anchor':(1,1)})
    axs[0].axis('off')

    colors = create_palette(df_, best, ten_godisnot)
    draw_embeddings(
        df_, cat=best, ax=axs[1],
        legend_kwargs={'loc':'upper left', 'bbox_to_anchor':(1,1), 'colors':colors}
    )
    axs[1].axis('off')

    colors = {'Blast-like':'r', 'Others':'k'}
    draw_embeddings(
        df_, cat='MT_clone_status', ax=axs[2],
        legend_kwargs={'loc':'upper left', 'bbox_to_anchor':(1,1), 'colors':colors}
    )
    axs[2].axis('off')

    fig.tight_layout()
    fig.savefig(
        os.path.join(path_results, 'umap_clones.png'),
        dpi=300
    )


    ##


    ############################## Viz vireoSNP clones variants

    # Viz individual variants-clone specificity
    make_folder(path_variants, 'variants', overwrite=True)
    vois_df = summary_stats_vars(a)

    # Here we go
    for var in vois_df.index:

        pos_fr = 1-vois_df.loc[var, "fr_positives"]
        median_AF = vois_df.loc[var, "median_AF"]

        fig, axs = plt.subplots(1,2, figsize=(13,5.5), constrained_layout=True)

        axs[0] = plot_exclusive_variant(a, var, vois_df, ax=axs[0])
        inset = axs[0].inset_axes((.3, .25, .4, .4))

        df_ = pd.DataFrame({'h':a[:, var].X.flatten()})
        hist(df_, 'h', n=sturges(df_['h']), c='k', l=.85, ax=inset)
        t = f'median_AF {median_AF:.2f}'
        format_ax(inset, title=t, xlabel='AF', ylabel='n cells', reduced_spines=True)

        df_ = expr_umap.join(df, how='right')
        df_[var] = pd.Series(a[:, var].X.toarray().flatten(), index=df_.index)
        order = df_.groupby(best)[var].agg('mean').sort_values(ascending=False).index

        box(df_, best, var, ax=axs[1], c='white', order=order)
        strip(df_, best, var, ax=axs[1], c='k', s=2, order=order)
        format_ax(axs[1], xlabel='MT clones', ylabel='AF', title='MT clones', reduced_spines=True)

        fig.suptitle(var)
        fig.savefig(
            os.path.join(
                path_variants, 'variants', f'{var}.png'
            ),
            dpi=300
        )

    ##


    # NB: no exclusive variants
    a_cells.obs['clone'] = clones[best]
    rank_clone_variants(a_cells, 'clone', min_clone_perc=.75, max_perc_rest=.25)


    ##


    ############################## UPMGA trees

    # Perfect meta
    a = a[df_filtered.index,:].copy()
    res = compute_clonal_fate_bias(df_filtered, 'aggregated_lineage', best, 'HSC_progenitors')
    blast_like_clones = res.query('FDR<=.05').index.astype('str')
    df_filtered['MT_clone_status'] = np.where(
        df_filtered[best].isin(blast_like_clones), 'Blast-like', 'Others'
    )

    # Define trees
    a = nans_as_zeros(a)
    D = pair_d(a, metric='cosine')
    char = pd.DataFrame(
        np.where(a.X>.1, 1, 0),
        index=a.obs_names,
        columns=a.var_names
    )
    tree = cs.data.CassiopeiaTree(
        character_matrix=char, 
        dissimilarity_map=pd.DataFrame(D, index=a.obs_names, columns=a.obs_names),
        cell_meta=df_filtered
    )
    solver = cs.solver.UPGMASolver()
    solver.solve(tree)

    # Fig
    fig, ax = plt.subplots(figsize=(5,5))
    cs.pl.plot_matplotlib(
        tree, meta_data=['MT_clone_status'],
        add_root=True,
        ax=ax
    )
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, 'UPMGA_tree.png'), dpi=300)


    ##


############################## 


# Run
if __name__ == '__main__':
    main()