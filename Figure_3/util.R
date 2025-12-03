volcano_plot = function(results_df, title) {
  library(ggrepel)
  library(cowplot)
  
  results_df_sub = results_df %>%
    mutate(logP = -log10(padj))
  
  results_df_sub$gene_name = rownames(results_df_sub)
  
  x_lim = max(abs(results_df_sub$log2FoldChange))
  
  lab_df = results_df_sub %>%
    filter(DE!=0)
  
  col_vec = c("#377EB8",
              "#d3d3d3",
              "#E41A1C")
  
  names(col_vec) = c(-1, 0, 1)
  
  p = ggplot(results_df_sub %>%
               arrange(desc(padj)),
             aes(x=log2FoldChange, y=logP, color=factor(DE))) +
    geom_point(size=0.5) +
    coord_cartesian(xlim=c(-x_lim, x_lim)) +
    geom_text_repel(data = lab_df,
                    aes(x = log2FoldChange, y = logP, label = gene_name),
                    min.segment.length = .5,
                    seed = 3,
                    box.padding = .2,
                    show.legend = FALSE,
                    max.overlaps =10,
                    size=4
    ) +
    scale_color_manual(values = col_vec[names(col_vec) %in% unique(results_df_sub$DE)]) +
    theme_cowplot(18) +
    labs(x="log2(FC)",
         y="-log10(FDR)") + ggtitle(title) +
    theme(legend.position = "none"
    )
  
  return(p)
}

call_DE_genes = function(logfc, fdr, abslogfc_cutoff, fdr_cutoff) {
  ifelse(logfc > abslogfc_cutoff & fdr < fdr_cutoff, 1,
         ifelse(logfc < (-abslogfc_cutoff) & fdr < fdr_cutoff, -1, 0))
}

# Perform GO GSEA (BP, CC, MF) using any OrgDb
GOgsea_fun <- function(res, OrgDb, keyType = "SYMBOL") {
  # res: DESeq2 results or data.frame with log2FoldChange and rownames as gene IDs
  # OrgDb: AnnotationDbi OrgDb object (e.g. org.Hs.eg.db, org.Mm.eg.db)
  # keyType: type of gene IDs in res (default "SYMBOL")
  
  # ensure required packages are available
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Please install the 'clusterProfiler' and 'AnnotationDbi' packages.")
  }
  
  # 1) Build a named, sorted vector of log2 fold-changes
  geneList <- res$log2FoldChange
  names(geneList) <- AnnotationDbi::mapIds(
    OrgDb,
    keys   = rownames(res),
    column = "ENTREZID",
    keytype = keyType,
    multiVals = "first"
  )
  geneList <- geneList[!is.na(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  
  # 2) Run GSEA for each GO ontology and make results human-readable
  ontologies <- c("BP", "CC", "MF")
  gsea_results      <- list()
  readable_results  <- list()
  for (ont in ontologies) {
    res_gsea <- clusterProfiler::gseGO(
      geneList     = geneList,
      OrgDb        = OrgDb,
      ont          = ont,
      minGSSize    = 50,
      maxGSSize    = 500,
      pvalueCutoff = 1,
      eps          = 0,
      verbose      = TRUE
    )
    # store original & readable versions
    gsea_results[[ont]]     <- res_gsea
    readable_results[[ont]]  <- clusterProfiler::setReadable(res_gsea, OrgDb = OrgDb)
  }
  
  # 3) Combine all readable results into one data.frame with an ONTOLOGY column
  simple_dfs <- lapply(names(readable_results), function(ont) {
    df <- as.data.frame(readable_results[[ont]])
    df$ONTOLOGY <- ont
    df
  })
  combined_readable <- do.call(rbind, simple_dfs)
  
  # 4) Return a structured list
  list(
    BP            = readable_results[["BP"]],
    BP_original   = gsea_results[["BP"]],
    CC            = readable_results[["CC"]],
    MF            = readable_results[["MF"]],
    ALL_readable  = combined_readable
  )
}


GOgsea_hs_fun <- function(res, orgdb, keytype = "SYMBOL") {
  # res is the output of results(dds)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  geneList <- res$log2FoldChange
  names(geneList) <- mapIds(org.Hs.eg.db,
                            keys=rownames(res),
                            column="ENTREZID",
                            keytype=keytype,
                            multiVals = "first")
  geneList <- geneList[!is.na(names(geneList))]
  
  geneList <- sort(geneList, decreasing = T)
  
  goGSEA_CC <- gseGO(geneList     = geneList,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "CC",
                     #nPerm        = 1000,
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     eps = 0,
                     verbose      = T)
  
  goGSEA_CC <- setReadable(goGSEA_CC, OrgDb = org.Hs.eg.db)
  
  goGSEA_BP <- gseGO(geneList     = geneList,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "BP",
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     eps = 0,
                     verbose      = T)
  
  goGSEA_BP_r <- setReadable(goGSEA_BP, OrgDb = org.Hs.eg.db)
  
  goGSEA_MF <- gseGO(geneList     = geneList,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "MF",
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     eps = 0,
                     verbose      = T)
  
  goGSEA_MF <- setReadable(goGSEA_MF, OrgDb = org.Hs.eg.db)
  
  goGSEA_BP_simple <- as.data.frame(goGSEA_BP_r)
  goGSEA_BP_simple$ONTOLOGY <- "BP" 
  
  goGSEA_CC_simple <- as.data.frame(goGSEA_CC)
  goGSEA_CC_simple$ONTOLOGY <- "CC"
  
  goGSEA_MF_simple <- as.data.frame(goGSEA_MF)
  goGSEA_MF_simple$ONTOLOGY <- "MF"
  
  goGSEA <- rbind(goGSEA_BP_simple,goGSEA_CC_simple,goGSEA_MF_simple)
  
  return(list(BP = goGSEA_BP_r, BP_original = goGSEA_BP, CC = goGSEA_CC, MF = goGSEA_MF, ALL_readable = goGSEA))
}


GOgsea_mmus_fun <- function(res, orgdb) {
  # res is the output of results(dds)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  
  geneList <- res$log2FoldChange
  names(geneList) <- mapIds(org.Mm.eg.db,
                            keys=rownames(res),
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals = "first")
  geneList <- geneList[!is.na(names(geneList))]
  
  geneList <- sort(geneList, decreasing = T)
  
  goGSEA_CC <- gseGO(geneList     = geneList,
                     OrgDb        = org.Mm.eg.db,
                     ont          = "CC",
                     #nPerm        = 1000,
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     eps = 0,
                     verbose      = T)
  
  goGSEA_CC <- setReadable(goGSEA_CC, OrgDb = org.Mm.eg.db)
  
  goGSEA_BP <- gseGO(geneList     = geneList,
                     OrgDb        = org.Mm.eg.db,
                     ont          = "BP",
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     eps = 0,
                     verbose      = T)
  
  goGSEA_BP <- setReadable(goGSEA_BP, OrgDb = org.Mm.eg.db)
  
  goGSEA_MF <- gseGO(geneList     = geneList,
                     OrgDb        = org.Mm.eg.db,
                     ont          = "MF",
                     minGSSize    = 50,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     eps = 0,
                     verbose      = T)
  
  goGSEA_MF <- setReadable(goGSEA_MF, OrgDb = org.Mm.eg.db)
  
  goGSEA_BP_simple <- as.data.frame(goGSEA_BP)
  goGSEA_BP_simple$ONTOLOGY <- "BP" 
  
  goGSEA_CC_simple <- as.data.frame(goGSEA_CC)
  goGSEA_CC_simple$ONTOLOGY <- "CC"
  
  goGSEA_MF_simple <- as.data.frame(goGSEA_MF)
  goGSEA_MF_simple$ONTOLOGY <- "MF"
  
  goGSEA <- rbind(goGSEA_BP_simple,goGSEA_CC_simple,goGSEA_MF_simple)
  
  return(list(BP = goGSEA_BP, CC = goGSEA_CC, MF = goGSEA_MF, ALL_readable = goGSEA))
}

#' Deconvolute bulk expression profiles using non‐negative least squares (NNLS)
#'
#' This function estimates cell‐type (or signature) proportions in bulk RNA‐seq samples
#' by solving a non‐negative least squares problem against a reference signature matrix.
#'
#' @param bulk_expression A matrix or data.frame of bulk expression values 
#'   (rows = samples, columns = genes).
#' @param reference_matrix A matrix of reference signatures 
#'   (rows = genes, columns = cell types or conditions).
#' @return A data.frame where each row corresponds to a sample and each column 
#'   corresponds to a reference signature; values are the estimated proportions.
#' @examples
#' # Assuming `bulk_expr` and `ref_mat` are already loaded:
#' props <- deconvolve_nnls(bulk_expr, ref_mat)
deconvolve_nnls <- function(bulk_expression, reference_matrix) {
  # Ensure the genes are aligned
  common_genes <- intersect(colnames(bulk_expression), rownames(reference_matrix))
  bulk_filtered <- as.matrix(bulk_expression[, common_genes])
  reference_filtered <- as.matrix(reference_matrix[common_genes, ])
  
  # Solve NNLS for each sample
  proportions <- t(apply(bulk_filtered, 1, function(sample) {
    fit <- nnls(reference_filtered, sample)
    return(fit$x)
  }))
  
  # Convert to DataFrame
  rownames(proportions) <- rownames(bulk_expression)
  colnames(proportions) <- colnames(reference_matrix)
  return(as.data.frame(proportions))
}

plot_volcano <- function(res, col_vec, ylim = c(0, 50), genes_to_label, title) {
  # Prepare data
  results_df <- res %>%
    rownames_to_column("gene_name") %>%
    mutate(
      label = as.integer(gene_name %in% genes_to_label),
      logP              = -log10(padj)
    ) %>%
    arrange(padj)
  
  # Top 10 hypoxia genes to label
  label_genes <- results_df %>%
    filter(label == 1) %>%
    slice_head(n = 10) %>%
    pull(gene_name)
  
  lab_df <- results_df %>% filter(gene_name %in% label_genes)
  
  # Define your color map for DE levels
  col_vec <- setNames(
    col_vec,
    c("-1", "0", "1")
  )[ as.character(unique(results_df$DE)) ]
  
  # Plot
  p <- ggplot(results_df, aes(log2FoldChange, logP, color = factor(DE))) +
    geom_point(size = 0.1) +
    geom_text_repel(
      data            = lab_df,
      aes(label = gene_name),
      min.segment.length = 0,
      force           = 100,
      seed            = 24,
      nudge_x         = 5,
      box.padding     = 0.2,
      max.overlaps    = 5,
      size            = 5
    ) +
    coord_cartesian(xlim = c(-3, 8), ylim = ylim) +
    scale_color_manual(values = col_vec) +
    labs(
      x     = "LFC",
      y     = "-log10(FDR)",
      title = title
    ) +
    theme_cowplot(10 * scale_factor) +
    theme(legend.position = "none")
  
  return(p)
}

#' Calculate Relative Quantification (RQ) From qPCR ΔΔCq Data
#'
#' This function takes sample metadata, raw EqCq values, and gene annotations
#' to compute ΔCq, ΔΔCq, and RQ (2^-ΔΔCq) values for each target gene.
#' It returns a long-form table of RQ values by sample, gene, and group.
#'
#' @param sample_data  Data frame with sample metadata. Must include:
#'   - Sample: sample ID
#'   - Biogroup: biological conditions/groups
#'   - Refgroup: indicator (1 = reference baseline group)
#'
#' @param rq_results   Data frame containing qPCR results with columns:
#'   - Sample
#'   - Target
#'   - EqCq.Mean (mean Cq normalized using efficiency calibration)
#'
#' @param gene_data    Data frame mapping targets to roles:
#'   - Target
#'   - EndoControl (1 = endogenous control gene, 0 = target gene)
#'
#' @return A long-form tibble with:
#'   - Sample
#'   - Biogroup
#'   - Target
#'   - RQ (relative expression, 2^-ΔΔCq)
#'   - logRQ (log2-transformed RQ)
#'
#' @details
#' Steps implemented:
#'   1. Reshape EqCq results to wide format (samples × targets)
#'   2. Compute ΔCq relative to the endogenous control gene
#'   3. Compute ΔΔCq relative to the reference (Refgroup = 1)
#'   4. Convert ΔΔCq → RQ using 2^-ΔΔCq
#'   5. Return tidy long-form dataset including log2-transformed RQ
#'
#' @examples
#' \dontrun{
#' rq_table <- calc_RQ(sample_data, rq_results, gene_data)
#' }
#'
calc_RQ <- function(sample_data, rq_results, gene_data) {
  delta_eqcq = rq_results %>%
    select(Sample, Target, EqCq.Mean) %>%
    pivot_wider(values_from = EqCq.Mean, names_from = Target) %>%
    as.data.frame()
  rownames(delta_eqcq) = delta_eqcq$Sample
  delta_eqcq$Sample = NULL
  endo_eqcq = delta_eqcq[[gene_data$Target[gene_data$EndoControl==1]]]
  delta_eqcq[[gene_data$Target[gene_data$EndoControl==1]]] = NULL
  delta_eqcq = apply(delta_eqcq, 2, function(x) {x-endo_eqcq})
  delta_delta_eqcq = delta_eqcq %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    left_join(sample_data, "Sample") %>%
    drop_na(Biogroup)
  for (i in gene_data$Target[gene_data$EndoControl==0]) {
    delta_delta_eqcq[[i]] = delta_delta_eqcq[[i]] - mean(delta_delta_eqcq[[i]][delta_delta_eqcq$Refgroup==1], na.rm = T)
  }
  ddCt = delta_delta_eqcq
  for (i in gene_data$Target[gene_data$EndoControl==0]) {
    delta_delta_eqcq[[i]] = 2**(-delta_delta_eqcq[[i]])
  }
  rq = delta_delta_eqcq %>%
    select(-c(Refgroup)) %>%
    pivot_longer(!c(Sample, Biogroup), names_to = "Target", values_to = "RQ") %>%
    drop_na(RQ)
  ddCt = ddCt %>%
    select(-c(Refgroup)) %>%
    pivot_longer(!c(Sample, Biogroup), names_to = "Target", values_to = "RQ")
  rq$Biogroup = factor(rq$Biogroup, levels = c(sample_data$Biogroup[sample_data$Refgroup==1][1], unique(sample_data$Biogroup[sample_data$Refgroup==0])))
  rq$logRQ = log2(rq$RQ)
  rq_n_samples = rq %>%
    group_by(Target, Biogroup) %>%
    summarise(n = n())
  rq_summarised = rq %>%
    group_by(Target, Biogroup) %>%
    left_join(rq_n_samples, by = c("Target", "Biogroup")) %>%
    summarise(RQ_geomean = exp(mean(log(RQ))), RQ_SE = sd(RQ)/n) %>%
    distinct()
  stat.test_ddct = ddCt %>%
    group_by(Target) %>%
    t_test(RQ ~ Biogroup) %>%
    adjust_pvalue(method="fdr") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Biogroup") %>%
    mutate(y.position = y.position) %>%
    mutate(p.adj = round(p.adj, digits = 3))
  return(rq)
}