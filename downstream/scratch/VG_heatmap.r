# VG heatmap

library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(Seurat)
library(Matrix)
library(ComplexHeatmap)
library(circlize) # for colorRamp2
library(readxl)
library(ggrepel)
library(gridExtra)


# Read data
path_corr <- '/Users/IEO5505/Desktop/AML_clonal_reconstruction/corr.csv'
cor.mat <- read.csv(path_corr, row.names=1) %>% as.matrix()

# Variant correlation & cluster
# cor.mat <- cor(t(af_subset.mat))
var.clust <- hclust(as.dist(1-cor.mat))

# Assess how correlated variants are and group them together
plot(var.clust$height, ylim = c(0, max(var.clust$height)))
# Make 23 groups of variants, i.e. 20 alone and 3 groups of two variants. This is determined empirically.
ngroups <- 5

hm1 <- Heatmap(cor.mat,
               col = colorRamp2(c(-1,0,1), c("blue", "#DDDDDD", "red")),
               cluster_columns = var.clust,
               cluster_rows = var.clust,
               # row_split = switch(ngroups < length(voi.ch), ngroups),
               # column_split = switch(ngroups < length(voi.ch), ngroups),
               show_row_dend = F, # without this the visualizationn does not complete
               show_column_dend = F, # without this the visualizationn does not complete
               row_gap = unit(0.5, "mm"),
               column_gap = unit(0.5, "mm"),
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10),
               row_title_gp = gpar(fontsize = 10),
               width = unit(100, "mm"),
               height = unit(100, "mm"),
               column_title = ngroups)
hm1

