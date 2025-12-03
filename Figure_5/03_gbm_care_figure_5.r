library(tidyverse)        # dplyr, tibble, stringr, readr, etc.
library(Seurat)           # CreateSeuratObject, NormalizeData, ScaleData, AddModuleScore, GetAssayData
library(ComplexHeatmap)   # Heatmap
library(RColorBrewer)     # brewer.pal
library(circlize)         # colorRamp2 (from circlize, auto-loaded by ComplexHeatmap)
library(grid)             # for gpar(), unit()

library(AnnotationDbi)    # select(), mapIds()
library(org.Hs.eg.db)     # human annotations
library(org.Mm.eg.db)     # mouse annotations (optional)
library(GO.db)            # GO terms database

source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_analysis_utils.R")
source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_CNA_utils.R")
source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_State_utils.R")

ANALYSIS_ROOT <- "/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/"
DATA_ROOT <- paste0(ANALYSIS_ROOT, "/data/")
TABLES_ROOT <- paste0(ANALYSIS_ROOT, "/tables/")
FIGURES_ROOT  <- paste0(ANALYSIS_ROOT, "/figures/")
SDS_SAMPLE_PATH <- paste0(DATA_ROOT, "/sample_sds/")
NMSDS_SAMPLE_PATH <- paste0(DATA_ROOT, "/sample_nmsds/")

# ===== Load data =====

meta_data <- readRDS(paste0(DATA_ROOT, "meta_data_qc.RDS"))

dirs_list <- list.files(DATA_ROOT, full.names = T)
dirs_list <- dirs_list[grepl("umi_counts.RDS", dirs_list)]

ct_classification <- readRDS("/users/mzarodn2/afs/Private/GSE274546/data/celltype_meta_data_2025_01_08.RDS")
meta_data <- meta_data %>% left_join(ct_classification %>% select(CellID, CellType), "CellID")


## ===== Load 1st sample =====

s <- dirs_list[1]
message(basename(s))
umi_data <- readRDS(s)

sample_name <- str_split(str_split(s, "/", Inf, T)[,12], "_", Inf, T)[,2]
meta_data_s <- meta_data[meta_data$Sample==sample_name,]

sds <- CreateSeuratObject(counts = umi_data[,meta_data_s$CellID],
                            project = "GBM-CARE",
                            min.cells = 0, min.features = 0, 
                            meta.data = meta_data_s)
gc()

sds <- NormalizeData(sds)
gc()

# ===== Define GO sets =====

#' Annotate gene symbols and identify protein-coding genes
#'
#' This function takes a vector of gene symbols and returns basic annotations,
#' including gene name, gene type, and whether the gene is protein-coding.
#' It supports both human and mouse annotations via `org.Hs.eg.db` and
#' `org.Mm.eg.db`.
#'
#' @param symbols A character vector of gene symbols.
#' @param species Either "human" or "mouse". Determines which annotation
#'   database to use. Defaults to "human".
#'
#' @return A tibble with columns:
#'   - `input_symbol`: input gene symbol
#'   - `GENENAME`: gene full name
#'   - `GENETYPE`: gene type returned by AnnotationDbi
#'   - `is_protein_coding`: TRUE/FALSE flag for protein-coding genes
#'
#' @examples
#' annotate_pc_by_symbol(c("TP53", "PDGFRB", "CDKN2A"), species = "human")
annotate_pc_by_symbol <- function(symbols, species = c("human","mouse")) {
  species <- match.arg(species)
  db <- if (species == "human") org.Hs.eg.db else org.Mm.eg.db

  res <- AnnotationDbi::select(
    db,
    keys     = unique(symbols),
    keytype  = "SYMBOL",
    columns  = c("SYMBOL","GENENAME","GENETYPE")   # GENETYPE: e.g. "protein-coding"
  )

  res <- res %>% distinct(SYMBOL, .keep_all = TRUE)

  tibble(input_symbol = symbols) %>%
    left_join(res, by = c("input_symbol" = "SYMBOL")) %>%
    mutate(is_protein_coding = tolower(GENETYPE) == "protein-coding") %>%
    dplyr::select(input_symbol, GENENAME, GENETYPE, is_protein_coding)
}

#' Retrieve gene symbols for a given GO term
#'
#' This function returns all unique human gene symbols annotated to a specified
#' Gene Ontology (GO) term, using the `org.Hs.eg.db` annotation database.
#'
#' @param goid A GO term ID (e.g., "GO:0006915") or a vector of GO IDs.
#'
#' @return A character vector of unique gene symbols associated with the GO term.
#'
#' @examples
#' genes_for_go("GO:0006915")
genes_for_go <- function(goid) {
  df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = goid,
    keytype = "GO",
    columns = c("ENTREZID","SYMBOL")
  )
  unique(df$SYMBOL)
}

#' Compute pairwise Jaccard similarity matrix for gene sets
#'
#' This function calculates the Jaccard index for all pairwise combinations
#' of gene sets provided in a named list. Each list element should be a vector
#' of gene symbols. The output is an n × n matrix where each entry represents
#' the Jaccard similarity between two gene sets.
#'
#' @param go_list A named list of gene vectors, where each element is a gene set.
#'
#' @return A numeric matrix with row and column names matching `go_list`,
#'   containing pairwise Jaccard similarity values.
#'
#' @examples
#' compute_jaccard_matrix(list(setA = c("A","B"), setB = c("B","C")))
compute_jaccard_matrix <- function(go_list) {
  n <- length(go_list)
  mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(names(go_list), names(go_list)))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      intersect_size <- length(intersect(go_list[[i]], go_list[[j]]))
      union_size     <- length(union(go_list[[i]], go_list[[j]]))
      mat[i, j] <- ifelse(union_size == 0, 0, intersect_size / union_size)
    }
  }
  mat
}

iN1_comp_d <- read.csv("data/iN1_neuron.csv", row.names = 1)
iN2_comp_d <- read.csv("data/iN2_neuron.csv", row.names = 1)

gene_annot <- annotate_pc_by_symbol(unique(c(rownames(iN1_comp_d), rownames(iN2_comp_d))), species = "human")
pc_genes <- gene_annot[gene_annot$is_protein_coding,]$input_symbol

iN1_comp_d <- iN1_comp_d[rownames(iN1_comp_d) %in% pc_genes,]
iN2_comp_d <- iN2_comp_d[rownames(iN2_comp_d) %in% pc_genes,]

iN1_comp_d_sig <- iN1_comp_d %>% rownames_to_column("Gene") %>% filter(padj < 0.05 & baseMean > 25) %>% arrange(desc(log2FoldChange)) %>% head(50) %>% pull("Gene")
iN2_comp_d_sig <- iN2_comp_d %>% rownames_to_column("Gene") %>% filter(padj < 0.1 & baseMean > 25) %>% arrange(desc(log2FoldChange)) %>% head(50) %>% pull("Gene")
iN_comp_d_sig <- c(iN1_comp_d_sig, iN2_comp_d_sig)

# ===== Figure 5L =====

# manually curated go gene sets related to cellular response to stress
stress_go_terms <- c(
  "cellular response to stress",
  "defense response",
  "multicellular organismal response to stress",
  "muscle hypertrophy in response to stress",
  "negative regulation of translational initiation in response to stress",
  "phage shock positive regulation of translational initiation in response to stress",
  "regulation of response to stress",
  "response to anoxia",
  "response to caloric restriction",
  "response to cold",
  "response to flooding",
  "response to fluid shear stress",
  "response to heat",
  "response to herbicide",
  "response to high population density",
  "response to hydrostatic pressure",
  "response to hyperoxia",
  "response to hypoxia",
  "response to immobilization stress",
  "response to ischemia",
  "response to isolation stress",
  "response to nitrosative stress",
  "response to osmotic stress",
  "response to oxidative stress",
  "response to psychosocial stress",
  "response to starvation",
  "response to sterol depletion",
  "response to topologically incorrect protein",
  "response to water deprivation",
  "response to wounding",
  "stress response to acid chemical",
  "stress response to metal ion"
)

go_df <- AnnotationDbi::select(
  GO.db,
  keys    = stress_go_terms,
  keytype = "TERM",
  columns = c("GOID","TERM")
)

go_keys <- keys(org.Hs.eg.db, keytype="GO")
go_df <- go_df[go_df$GOID %in% go_keys,]
go_list <- lapply(go_df$GOID, genes_for_go)
names(go_list) <- go_df$TERM
go_list[["response to solid stress"]] <- iN_comp_d_sig
go_list <- lapply(go_list, function(x) {x[x %in% rownames(sds)]})
go_list <- go_list[lengths(go_list)>=5]

jaccard_mat <- compute_jaccard_matrix(go_list)

pdf("figures/jaccard_mat.pdf", width = 7.5, height = 4.75)
ht = Heatmap(jaccard_mat,
        name = "Jaccard",
        col = colorRamp2(c(0, 0.2), c("white", "darkgreen")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        width  = unit(6, "cm"),       # body width
        height = unit(6, "cm"),
        border = T,
        show_column_dend = F,
        show_column_names = FALSE,
        show_row_dend = TRUE,)
draw(ht, padding = unit(c(6, 6, 6, 6), "cm"), heatmap_legend_side = "left")  # top, right, bottom, left
dev.off()

## ===== Overlap counts ======

target_name <- "response to solid stress"
genes_orig <- unique(go_list[[target_name]])
other_sets <- go_list[setdiff(names(go_list), target_name)]
present_in <- lapply(genes_orig, function(g)
  names(Filter(function(v) g %in% v, other_sets))
)

res <- data.frame(
  gene = genes_orig,
  n_sets = vapply(present_in, length, integer(1)),
  sets = vapply(present_in, function(x) paste(x, collapse = "; "), character(1)),
  total_other_sets = length(other_sets),
  fraction = vapply(present_in, length, integer(1)) / max(1L, length(other_sets)),
  stringsAsFactors = FALSE
)

res <- res[order(-res$n_sets, res$gene), ]
message(paste0("Genes with overlap across sets: ", paste0(res$gene[res$n_sets>0], collapse = ", ")))

iN_comp_d_sig_rev <- go_list[["response to solid stress"]][!(go_list[["response to solid stress"]] %in% res$gene[res$n_sets>0])]

# ===== Compute signature scores =====

## ===== Define gene sets ======

# manually curated gene sets realated to neuron physiology
go_terms <- c(
  "cholinergic synapse",
  "dopaminergic synapse",
  "excitatory synapse",
  "GABA-ergic synapse",
  "glutamatergic synapse",
  "glycinergic synapse",
  "inhibitory synapse",
  "neuron to neuron synapse",
  "noradrenergic synapse",
  "polyadic synapse",
  "postsynapse",
  "presynapse",
  "serotonergic synapse",
  "synaptic membrane",
  "synapse",
  "postsynaptic signal transduction",
  "signal release from synapse",
  "synaptic signaling by nitric oxide",
  "synaptic signaling via neuropeptide",
  "synaptic vesicle transport",
  "trans-synaptic signaling",
  "cholesterol metabolic process",
  "phospholipid biosynthetic process",
  "lipid oxidation",
  "antioxidant activity",
  "response to hypoxia"
)

go_df <- AnnotationDbi::select(
  GO.db,
  keys    = go_terms,
  keytype = "TERM",
  columns = c("GOID","TERM")
)

go_df <- go_df[go_df$GOID %in% go_keys,]
go_list <- lapply(go_df$GOID, genes_for_go)
names(go_list) <- paste0(gsub("\\ ", "_", go_df$TERM), "_", go_df$GOID)
go_list[["response_to_solid_stress"]] <- iN_comp_d_sig_rev
go_list <- lapply(go_list, function(x) {x[x %in% rownames(sds)]})
go_list <- go_list[lengths(go_list)>=10]

## ===== Compute scores ======

#sds_list <- list()
#scores_df <- data.frame()
#
#for (i in 1:length(dirs_list)) {
#    s <- dirs_list[i]
#
#    message(paste0("[", i, "/", length(dirs_list), "]: ", basename(s)))
#    umi_data <- readRDS(s)
#
#    sample_name <- str_split(str_split(s, "/", Inf, T)[,12], "_", Inf, T)[,2]
#    meta_data_s <- meta_data[meta_data$Sample==sample_name,]
#
#    sds <- CreateSeuratObject(counts = umi_data[,meta_data_s$CellID],
#                                project = "GBM-CARE",
#                                min.cells = 0, min.features = 0, 
#                                meta.data = meta_data_s)
#    gc()
#
#    sds <- NormalizeData(sds)
#    gc()
#
#    sds <- AddModuleScore(object = sds, features = go_list, name = "NormalSigs")
#
#    d <- sds@meta.data %>% dplyr::select(CellID, Sample, CellType, starts_with("NormalSigs")) %>% filter(CellType %in% c("Excitatory neuron", "Inhibitory neuron"))
#
#    if (nrow(d)<10) next
#
#    colnames(d)[grep("NormalSigs", colnames(d))] <- names(go_list)
#
#    scores_df <- rbind(scores_df, d)
#
#    sds_list[[s]] <- sds[unlist(go_list), sds$CellType %in% c("Excitatory neuron", "Inhibitory neuron")]
#
#    gc()
#}
#
# Cache the results
# saveRDS(scores_df, "data/nloss_module_scores_filt.RDS")
# sds_merged <- Reduce(function(x, y) merge(x, y, add.cell.ids = NULL, project = "MergedProject"), sds_list)
# saveRDS(sds_merged, "data/nloss_sds.RDS")

# ===== Figure S6E =====

scores_df <- readRDS("./data/nloss_module_scores_filt.RDS")
sds_merged <- readRDS("data/nloss_sds.RDS")
sds_merged <- ScaleData(sds_merged)

X <- GetAssayData(sds_merged[,sds_merged$CellType=="Excitatory neuron"], layer = "scale.data")
X_sub <- X[go_list$response_to_solid_stress,]

pal <- brewer.pal(11, "Spectral")
stress_col_fun <- colorRamp2(
  seq(from = -0.1, to = 0.1, length.out = 11),
  rev(pal)
)

top_anno <- HeatmapAnnotation(
  Signature = scores_df[colnames(X_sub),]$response_to_solid_stress,
  col = list(
    Signature = stress_col_fun
  ),
  annotation_name_side = "left"
)

pdf("figures/excitatory_solidstress_heatmap.pdf", width = 6, height = 6)
Heatmap(X_sub, show_row_dend = F, show_column_dend = F, show_column_names = F, name="z-score", use_raster=FALSE, top_annotation = top_anno)
dev.off()

# ===== Figure 5M ======

scores_df_sub <- scores_df %>% filter(CellType == "Excitatory neuron")

sample_names <- scores_df_sub %>% group_by(Sample) %>% summarise(count = n()) %>% filter(count > 10) %>% pull(Sample)

cor_list <- lapply(sample_names, function(x) {
  y <- scores_df %>% filter(Sample == x)
  return(cor(y[setdiff(all_cols, c("CellID","Sample","CellType"))], method = "spearman"))
})

# Compute element-wise average
cor_mean <- Reduce("+", cor_list) / length(cor_list)

colnames(cor_mean)  <- gsub("_", "\\ ", str_split(colnames(cor_mean), "_GO", 2, T)[,1])
rownames(cor_mean)  <- gsub("_", "\\ ", str_split(rownames(cor_mean), "_GO", 2, T)[,1])

col_fun <- colorRamp2(
  seq(-0.5, 0.5, length.out = 11),
  rev(brewer.pal(11, "RdBu"))
)

od <- hclust(dist(cor_mean))$order
cor_mean <- cor_mean[od, od]

ht <- Heatmap(
  cor_mean,
  col = col_fun,
  show_column_names = FALSE,
  show_row_dend = TRUE,
  width  = unit(6, "cm"),       # body width
  height = unit(6, "cm"),
  border = T,
  show_column_dend = FALSE,
  name = "correlation"
)

pdf("figures/excitatory_cor_heatmap.pdf", width = 7.5, height = 4.75)
draw(ht, padding = unit(c(2, 2, 7, 1), "cm"), heatmap_legend_side = "left")
dev.off()