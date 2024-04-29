# Phylogeny inferece for AML primary patients (SCMseq cohort)

This is a repo hosting code for the SCMseq project with Chiara Caprioli. 
The goal is to identify bona-fide cellular clones from MT-SNVs and to assess their relationship with
transcriptional genetic (somatic drivers identified from WES) states.

The repo is organized as follows:

* __downstream__:

```bash
├── expression
│   ├── dot_plot.py
│   └── expression_UMAP.py
├── muts
│   ├── heatmap_SCM_muts.py
│   └── prova_new_muts_sAML1.py
├── old
│   └── df_vars_generation.py
├── phylophenotype
│   ├── PATH.r
│   ├── PATH_prep_categorical.py
│   ├── PATH_prep_hvgs.py
│   ├── aggregated_ct_t_phylocorr.py
│   ├── malignant_supervised_phylocorr.py
│   └── unsupervised_phylopheno.py
├── sAML1
│   ├── DE_top_clones.py
│   └── top_clones_viz.py
├── scratch
│   └── MI_TO_filter.py
└── trees
    ├── tree_reconstruction_performance.py
    └── tree_viz.py
```

The __downstream__ folder holds all the necessary code for downstream expression (expression), nuclear SNVs (muts), phylogenetic (trees)
and phylophenotypic (phylophenotype) analysis. Other folders (sAML1, old, scrach) contains code for sample specific analyses,
plus deprecated and actively developed scripts/utils. 

The data for each of these analyses is read from an input data folder with the following structure:

* __data__:

```bash
├── MT_vars
│   ├── AML2_vars.csv
│   ├── AML3_vars.csv
│   ├── AML5_vars.csv
│   └── sAML1_vars.csv
├── expression
│   ├── HVG_final.csv
│   ├── aggr_lin_markers.csv
│   ├── malignant_AUCell.csv
│   ├── malignant_signatures.tsv
│   ├── sAML1_malignant_lognorm_counts.csv
│   └── umap_coord.csv
├── meta
│   ├── cells_meta.csv
│   ├── colors.pickle
│   └── old_meta
│       ├── all_samples_meta.csv
│       ├── cc_cells_meta.csv
│       ├── fab4_meta.csv
│       ├── genotype_all_final.csv
│       ├── meta.csv
│       └── umap.csv
├── muts
│   ├── genotype_all.csv
│   ├── genotypes.csv
│   ├── output_sAML1_new.txt
│   └── output_sAML1_old.txt
├── other
│   └── hg38.full.blacklist.bed
└── sAML1
    ├── AFM.h5ad
    └── barcodes.txt
```

and results are written on a folder whose structure is in refactoring. See individual analysis code for specifics.
To reproduce the computing environment needed to reproduce the analysis, follows instruction at https://github.com/andrecossa5/mito_utils.

