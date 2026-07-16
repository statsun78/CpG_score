# CpG_score
R codes for CpG score-based Gene Set Enrichment Analysis (GSEA).

Please download `function.R` to compute the CpG-level score.

---

## Functions

### 1. `gene_scores`
It computes the gene-level score by aggregating individual CpG site-level scores annotated to the corresponding genes.

#### Arguments:
* `pval_m`: A vector of p-values testing the mean difference for all CpG sites.
* `pval_v`: A vector of p-values testing the variance difference for all CpG sites.
* `gID`: A list of genes annotated to the CpG sites.

---

### 2. `weight_gene`
It computes the weight of individual genes when an adjacency matrix for a genetic network graph is provided.

#### Arguments:
* `score`: Gene-level scores computed by `gene_scores`.
* `adjm`: An adjacency matrix for the genetic network.
* `weight`: Weighting method to use:
  * `"degree"`: for degree-based weight.
  * `"topology"`: for topological-based weight.


 
