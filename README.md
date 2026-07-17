# CpG_score
R codes for CpG score-based Gene Set Enrichment Analysis (GSEA).

Please download `function.R` to compute the gene-level scores.

---

## Functions

### 1. `gene_scores`
It computes the gene-level score by aggregating individual CpG site-level scores annotated to the corresponding genes.

#### Arguments:
* `pval_m`: A vector of p-values testing the mean difference for all CpG sites.
* `pval_v`: A vector of p-values testing the variance difference for all CpG sites.
* `gID`: A list of genes annotated to the CpG sites.
* `eps`: The minimum p-value allowed. (Default: `1e-50`)

---

### 2. `weight_gene`
It computes the weight of individual genes when an adjacency matrix for a genetic network graph is provided.

#### Arguments:
* `score`: Gene-level scores computed by `gene_scores`.
* `adjm`: An adjacency matrix for the genetic network.
* `weight`: Weighting method to use:
  * `"degree"`: for degree-based weight.
  * `"topology"`: for topological-based weight.

---

## Example

```r
# 1. Define CpG sites and Gene mapping
cpg <- paste("c", 1:10, sep="")
gene <- c("gene1", "gene2", "gene3")

gID <- vector("list", length(gene))
gID[[1]] <- cpg[1:3]
gID[[2]] <- cpg[4:8]
gID[[3]] <- cpg[7:10]
names(gID) <- gene

# 2. Generate mock p-values for CpG sites
set.seed(123)
pval_m <- runif(10, 0, 0.001)
pval_v <- runif(10, 0, 0.001)
names(pval_m) <- names(pval_v) <- cpg

# 3. Compute gene-level scores
score <- gene_scores(pval_m, pval_v, gID)

# 4. Compute network-weighted gene scores
adjm <- matrix(0, length(gene), length(gene))
adjm[1, 3] <- adjm[3, 1] <- adjm[2, 3] <- adjm[3, 2] <- 1
wscore <- weight_gene(score, adjm)

# 5. Define pathways for GSEA
sID <- vector("list", 2)
sID[[1]] <- gene[1:2]
sID[[2]] <- gene[2:3]
names(sID) <- c("path1", "path2")

# 6. Run fgsea
library(fgsea)
fgsea::fgsea(pathways = sID, stats = score, scoreType = "pos", eps = 0)
fgsea::fgsea(pathways = sID, stats = wscore, scoreType = "pos", eps = 0)
```

---

## Reference
Lee, M., Yoon S., and H. Sun (2026) A CpG score-based gene set enrichment test for analysis of DNA methylation data with a hierarchical structure, *submitted*.



 
