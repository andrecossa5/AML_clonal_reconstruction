# Clonal reconstruction for the SCMseq project

This is a repo hosting code for the SCMseq project with Chiara Caprioli et al. From existing Diagnostic and matched Normal bulk WES data, 
the goal is to:

1. (Re-)run Mutect2 to call SNVs.
2. Filter the raw variant call (after comparison with the previous one) and see which variants (how many?) are selected. Can we get more?
3. WIth some (definitive, with Mutect2) variants list at hand, prepare pyClone inputs and run it.
4. Assess variants presence in single-cell reads from 10x.


