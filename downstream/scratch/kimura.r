# Kimura fitting

# Code
library(parallel)
library(tidyverse)
library(data.table)
library(kimura)

# Read data
path_data <- '/Users/IEO5505/Desktop/AML_clonal_reconstruction/results/var_selection/prova_kimura.csv'
df <- fread(path_data) %>% as.data.frame()
rownames(df) <- df[[1]]  # Set the first column as row names
df <- df[,-1, drop = FALSE]

# Kimura test
f <- function(x) { test_kimura(x)$p } 
results <- mclapply(df, f, mc.cores=8)
results <- results[sapply(results, function(x) !inherits(x, "try-error"))]

# Print significative
sig_variants <- names(which(results<=0.05)) 
length(sig_variants)


##

