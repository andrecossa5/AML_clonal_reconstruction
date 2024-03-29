# PATH (Landau pre-print) analysis

# Code
library(tidyverse)
library(parallel)
library(Matrix)
library(magrittr)
library(ape)
library(fgsea)
library(PATH)


##


# Bug fix
parM_I <- function(d, w, break.length=100) {
  
  y <- c(seq(1,ncol(d), break.length), (ncol(d)+1) )
  y2 <- y[-1] - 1
  y1 <- y[-length(y)]
  
  mfunc <- function(b, e, d1=d, w1=w) {
    M <- xcor(d1[,b:e ], w1)
    x1 <- diag(M$Z)
    x2 <- diag(M$Morans.I)
    out <- data.frame("Z"=x1, "I"=x2)
    return(out)
  }
  
  out <- do.call(rbind,
                mcmapply(function(b1,e1) mfunc(b=b1, e=e1, d, w), b1=y1, e1=y2,
                        mc.cores=detectCores(), SIMPLIFY=F))
  return(out)
  
}


##


xcor_gsea_I <- function(gene.data, weight.matrix, pathways, maxS=Inf, minS=5, nperm=10000) {
  
  z <- parM_I(gene.data, weight.matrix)
  z$gene <- rownames(z)
  
  z <- z %>% as_tibble() %>% group_by(gene) %>% 
    dplyr::select(Z) %>% arrange(desc(Z)) %>%
    drop_na() %$% set_names(Z, gene)
  out <- fgsea(pathways, z, maxSize=maxS, minSize=minS, scoreType="std", nPermSimple=nperm)
  res <- list(z=z, gsea=out)
  
  return(res)
  
}


##


# Set args
cov <- 'malignant_unsupervised'
sample <- 'AML5'
path_sample <- paste0('/Users/IEO5505/Desktop/AML_clonal_reconstruction/results/trees/top_trees/', sample)


##


# Load data

# Tree
tree <- ape::read.tree(paste0(path_sample, '/tree_', cov, '.newick'))
# Expression
hvgs <- read.csv(paste0(path_sample, '/HVGs_', cov, '.csv'), row.names = 1)[tree$tip.label,]
hvgs <- scale(hvgs) %>% as.matrix()
# Cat variable, to numeric
cov_df <- read.csv(paste0(path_sample, '/cov_df_', cov, '.csv'), row.names = 1)[tree$tip.label,]
table(cov_df)
X <- state2mat.sparse(cov_df)


##


# Phylo-correlation categorical

# Phylogenetic weight matrix
Winv <- inv.tree.dist(tree, node=TRUE, norm=FALSE)
modxcor <- xcor(X, Winv)
diag(modxcor$Z.score)
diag(modxcor$one.sided.pvalue)
Idf <- reshape2::melt(modxcor$Morans.I, value.name = "I")
Zdf <- reshape2::melt(modxcor$Z.score, value.name = "Z")
pdf <- reshape2::melt(modxcor$one.sided.pvalue, value.name = "p")
dfs <- list(Idf, Zdf, pdf)
results <- Reduce(
  function(x, y) full_join(x, y, by = c("Var1", "Var2")), 
  dfs
  ) %>%
  mutate(sample=sample) %>% 
  rename(GBC1=Var1, GBC2=Var2)

# Save
write.csv(results, paste0(path_sample, '/phylocorr_', cov, '.csv'))


##


# T inference
tree$states <- cov_df     
Pinf <- PATH.inference(tree, impute_branches=T, sample_rate_est=10^-6)
Ptinf.df <- data.frame(Pinf$Pt)
write.csv(Ptinf.df, paste0(path_sample, '/rates_', cov, 't_inference.csv'))


##


# Load HVGs expression

# Unsupervised discovery of hereditable gene modules
genesets = msigdbr::msigdbr(species='Homo sapiens', category='H')
pathwaysH = split(x=genesets$gene_symbol, f=genesets$gs_name)
pathwaysH <- lapply(pathwaysH, unique)

# Correlation with gene sets
Winv <- inv.tree.dist(tree, node=TRUE, norm=FALSE)
res <- xcor_gsea_I(hvgs, Winv, pathwaysH)
write.csv(data.frame(Z=res$z), paste0(path_sample, '/ranked_genes_', cov, '.csv'))
gsea_df <- res$gsea %>% arrange(desc(-padj)) %>% mutate(rank=rank(padj)) 
write.csv(gsea_df[,1:6], paste0(path_sample, '/H_gsea_', cov, '.csv'))


##