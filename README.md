## CpG_score
(R codes for CpG score-based GSEA)

Please downlaod `function.R` to compute CpG-level score.

* The R function `gene_scores`
  It computes the gene level score, aggregating individual CpG site level scores annotated to the correponding genes.
  - `pval_m` : A vector of the p-values to test a mean differerence for all CpG sites.
  - `pval_v` : A vector of the p-values to test a variance differerence for all CpG sites.
  - `gID` : A list of genes annotated CpG sites.

* The R function `wieght_gene` 


 
