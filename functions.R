library(dplyr)

gene_scores <- function(pval_m, pval_v, gID) {

    ## =========================================================
    ## Step 1: Convert p-values to z-scores
    ## =========================================================
    pm <- pmin(pmax(pval_m, eps), 0.5)
    pv <- pmin(pmax(pval_v, eps), 0.5)
    ms <- qnorm(pm, lower.tail = FALSE)
    vs <- qnorm(pv, lower.tail = FALSE)
    lambda <- mean(ms / (ms + vs), na.rm = TRUE)
    cpg_score <- lambda * ms + (1 - lambda) * vs

    ## =========================================================
    ## Step 2: Check whether CpGs overlap across genes
    ## =========================================================
    all_cpgs <- unlist(gID, use.names = FALSE)
    duplicated_cpgs <- unique(all_cpgs[duplicated(all_cpgs)])

    ## =========================================================
    ## Case 1: No duplicated CpGs
    ## =========================================================
    if (length(duplicated_cpgs) == 0) {
        med <- median(lengths(gID))
        gene_score <- vapply(gID,
            function(ids) {
                scores <- cpg_score[ids]
                scores <- scores[!is.na(scores)]
                if (length(scores) == 0) return(0)
                med_sum(scores, med = med)
            },
            numeric(1)
        )
        return(gene_score)
    }

    ## =========================================================
    ## Case 2: Duplicated CpGs exist
    ## =========================================================
    mapping_df <- data.frame(
        cpg  = unlist(gID, use.names = FALSE),
        gene = rep(names(gID), lengths(gID))
    )

    ## ---------------------------------------------------------
    ## Create CpG sharing patterns
    ## ---------------------------------------------------------
    cpg_patterns <- mapping_df |>
        dplyr::group_by(cpg) |>
        dplyr::summarise(
            gene_pattern = paste(sort(unique(gene)),
                                 collapse = "|"),
            .groups = "drop"
        )

    ## ---------------------------------------------------------
    ## Group CpGs by gene and sharing pattern
    ## ---------------------------------------------------------
    gene_cpg_groups <- mapping_df |>
        dplyr::left_join(cpg_patterns, by = "cpg") |>
        dplyr::group_by(gene, gene_pattern) |>
        dplyr::summarise(
            cpgs = list(cpg),
            group_size = dplyr::n(),
            .groups = "drop"
        )

    ## ---------------------------------------------------------
    ## Compute group-level scores
    ## ---------------------------------------------------------
    med <- median(gene_cpg_groups$group_size)
    group_scores <- vapply(
        gene_cpg_groups$cpgs,
        function(ids) {
            scores <- cpg_score[ids]
            scores <- scores[!is.na(scores)]
            if (length(scores) == 0) return(0)
            med_sum(scores, med = med)
        },
        numeric(1)
    )

    ## ---------------------------------------------------------
    ## Aggregate to gene-level scores
    ## ---------------------------------------------------------
    gene_scores_df <- data.frame(
        gene       = gene_cpg_groups$gene,
        score      = group_scores,
        group_size = gene_cpg_groups$group_size
    ) |>
        dplyr::group_by(gene) |>
        dplyr::summarise(
            gene_score = weighted.mean(score,
                                       w = group_size),
            .groups = "drop"
        )

    stats::setNames(
        gene_scores_df$gene_score,
        gene_scores_df$gene
    )
}

med_sum <- function(x, med) {
    n <- length(x)
    if (n == 0) {
        return(0)
    }
    sum(x) + mean(x) * (med - n)
}


weight_gene <- function(score, adjm, weight = c("degree", "topology")) {
    weight <- match.arg(weight)
    diag(adjm) <- 0
    centrality <- switch(
        weight,
        degree = {
            rowSums(adjm)
        },
        topology = {
            gg <- igraph::graph_from_adjacency_matrix(
                adjm,
                mode = "undirected",
                diag = FALSE
            )
            igraph::betweenness(
                gg,
                directed = FALSE,
                normalized = FALSE
            )
        }
    )
    ## Avoid division by zero
    max_centrality <- max(centrality, na.rm = TRUE)
    if (max_centrality == 0) {
        weighted_score <- score
    } else {
        weighted_score <- score * (1 + centrality / max_centrality)
    }
    stats::setNames(weighted_score, names(score))
}
