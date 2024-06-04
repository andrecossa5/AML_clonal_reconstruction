
library(tidyverse)

SEED <- 4321
set.seed(SEED)

path_variants <- fs::path("/Users/ieo6983/Desktop/mito/data/functional/variants.csv")
path_anno_table <- fs::path("/Users/ieo6983/Desktop/mito/data/functional/formatted_table_wobble.csv")
path_results_plots <- fs::path("/Users/ieo6983/Desktop/mito/AML_clonal_reconstruction/results/")

file_name <- fs::path(path_results_plots, "Rplots.mt_variants_functional_annotation.pdf")
pdf(file_name, width = 7, height = 6)

##


# Annotated variants from https://github.com/EDePasquale/Mitochondrial_variants 
anno_table <- read.csv(path_anno_table, row.names = "X")

#' Variables of interest: 
#' 
#' @Biotype
#' - type of transcript or regulatory feature: protein coding, tRNA, or rRNA.
#' @Consequence 
#' - Predicted consequence of the variant
#' - levels, i.e.: missense_variant, stop_gained, stop_lost, ...
#' @SIFT: 
#' - SIFT prediction score for whether an amino acid change will affect protein function
#' - levels, i.e.: tolerated, deleterious, tolerated_low_confidence
#' @PolyPhen
#' - PolyPhen prediction score for whether an amino acid change will affect protein function
#' - levels: "benign", "possibly_damaging", "probably_damaging", "unknown"    
#' @syn_annotation: 
#' - levels: WCF_to_WCF, WCF_to_Wobble, Wobble_to_WCF, Wobble_to_Wobble
#' - glossary: WCF = Watson & Crick; Wobble = chnage in codon usage 
#' 
#' Note: Any cells containing "-" are unknown values

vars_of_int <- c("Biotype", "Consequence", "SIFT", "PolyPhen", "syn_annotation")

# Format variables names - removing (value)
anno_table$SIFT <- str_split(anno_table$SIFT, "\\(", simplify = T)[,1]
anno_table$PolyPhen <- str_split(anno_table$PolyPhen, "\\(", simplify = T)[,1]


##


dataset <- read.csv("/Users/ieo6983/Desktop/mito/data/functional/dataset.csv")
variants <- read.csv(path_variants, row.names = "X")

# Format mutations names to: m9990T>C
# m + position + reference + > + variant
variants <- variants %>% mutate(mutation = str_c("m", row.names(variants)) %>% str_replace("_", "")) %>%
  relocate(., mutation, 1)

# Annotated variants
variants_anno <- left_join(variants, anno_table, by = "mutation") 

# Plot frequency of variants annotations
df <- variants_anno %>% dplyr::select(c("mutation", all_of(vars_of_int)))
plots <- lapply(names(df[,-1]), FUN = function(var){
  freq_table <- table(df[, var])
  levels_ordered <- names(sort(-freq_table))
  df[, var] <- factor(df[,var], levels = levels_ordered)
  
  p <- ggplot(df, aes_string(x = var, fill = var)) +
    geom_bar() +
    scale_fill_brewer(palette = "Set2")+
    ggtitle(paste("Frequency of", var)) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  print(p)
}
)

#plots
#grid_plot <- do.call(grid.arrange, c(plots, ncol = 1)) 
#ggsave(grid_plot, )


##

dev.off()




