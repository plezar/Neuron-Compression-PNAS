source("Figure_3/util.R")
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)

# ==== iN =====

## ==== GO GSEA =====

res_Wbo2 <- read.csv("Figure_3/results/iN1_neuron.csv", row.names = 1)
res_I27 <- read.csv("Figure_3/results/iN2_neuron.csv", row.names = 1)
res_1019 <- read.csv("Figure_3/results/iN3_neuron.csv", row.names = 1)

goGSEA_Wbo2 <- GOgsea_fun(res_Wbo2, org.Hs.eg.db)
goGSEA_I27 <- GOgsea_fun(res_I27, org.Hs.eg.db)
goGSEA_1019 <- GOgsea_fun(res_1019, org.Hs.eg.db)

goGSEA_merged <- rbind(
  as.data.frame(goGSEA_Wbo2$ALL_readable) %>% mutate(CellLine = "iN #1"),
  as.data.frame(goGSEA_1019$ALL_readable) %>% mutate(CellLine = "iN #3"),
  as.data.frame(goGSEA_I27$ALL_readable) %>% mutate(CellLine = "iN #2")
)

write.csv(goGSEA_merged, "Figure_3/results/goGSEA_merged_iN.csv")

## ==== HALLMARK =====

pathways <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

colors <- brewer.pal(length(pathways), "Set1")
hallmark_colors <- setNames(colors, pathways)

geneList_Wbo2 <- res_Wbo2$log2FoldChange
names(geneList_Wbo2) <- mapIds(org.Hs.eg.db,
                               keys=rownames(res_Wbo2),
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals = "first")

geneList_I27 <- res_I27$log2FoldChange
names(geneList_I27) <- mapIds(org.Hs.eg.db,
                              keys=rownames(res_I27),
                              column="ENTREZID",
                              keytype="SYMBOL",
                              multiVals = "first")

overlapping_genes <- intersect(names(geneList_Wbo2), names(geneList_I27))

geneList_Wbo2 <- sort(geneList_Wbo2[overlapping_genes], decreasing = T)
geneList_Wbo2 <- geneList_Wbo2[!is.na(names(geneList_Wbo2))]

geneList_I27 <- sort(geneList_I27[overlapping_genes], decreasing = T)
geneList_I27 <- geneList_I27[!is.na(names(geneList_I27))]

m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

hallmarkGSEA_Wbo2 <- GSEA(geneList_Wbo2, TERM2GENE = m_t2g, pvalueCutoff=1)
hallmarkGSEA_I27 <- GSEA(geneList_I27, TERM2GENE = m_t2g, pvalueCutoff=1)

hallmarkGSEA_merged <- rbind(
  as.data.frame(hallmarkGSEA_Wbo2) %>% mutate(CellLine = "WBO2"),
  as.data.frame(hallmarkGSEA_I27) %>% mutate(CellLine = "I27")
)

write.csv(hallmarkGSEA_merged, "Figure_3/results/HALLMARK_merged_iN.csv")

save(hallmarkGSEA_Wbo2, hallmarkGSEA_I27, file = "Figure_3/results/HALLMARK_merged_iN.RData")


# ==== mGlia =====

## ==== GO GSEA =====

res_Wbo2 <- read.csv("Figure_3/results/iN1_glia.csv", row.names = 1)
res_I27 <- read.csv("Figure_3/results/iN2_glia.csv", row.names = 1)
res_1019 <- read.csv("Figure_3/results/iN3_glia.csv", row.names = 1)

goGSEA_Wbo2 <- GOgsea_fun(res_Wbo2, org.Mm.eg.db)
goGSEA_I27 <- GOgsea_fun(res_I27, org.Mm.eg.db)
goGSEA_1019 <- GOgsea_fun(res_1019, org.Mm.eg.db)

goGSEA_merged <- rbind(
  as.data.frame(goGSEA_Wbo2$ALL_readable) %>% mutate(CellLine = "mG #1"),
  as.data.frame(goGSEA_1019$ALL_readable) %>% mutate(CellLine = "mG #3"),
  as.data.frame(goGSEA_I27$ALL_readable) %>% mutate(CellLine = "mG #2")
)

write.csv(goGSEA_merged, "Figure_3/results/goGSEA_merged_mG.csv")


## ==== HALLMARK =====

geneList_Wbo2 <- res_Wbo2$log2FoldChange
names(geneList_Wbo2) <- mapIds(org.Mm.eg.db,
                               keys=rownames(res_Wbo2),
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals = "first")

geneList_I27 <- res_I27$log2FoldChange
names(geneList_I27) <- mapIds(org.Mm.eg.db,
                              keys=rownames(res_I27),
                              column="ENTREZID",
                              keytype="SYMBOL",
                              multiVals = "first")

overlapping_genes <- intersect(names(geneList_Wbo2), names(geneList_I27))

geneList_Wbo2 <- sort(geneList_Wbo2[overlapping_genes], decreasing = T)
geneList_Wbo2 <- geneList_Wbo2[!is.na(names(geneList_Wbo2))]

geneList_I27 <- sort(geneList_I27[overlapping_genes], decreasing = T)
geneList_I27 <- geneList_I27[!is.na(names(geneList_I27))]

m_df <- msigdbr(species = "Mus musculus")
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, entrez_gene)

hallmarkGSEA_Wbo2 <- GSEA(geneList_Wbo2, TERM2GENE = m_t2g, pvalueCutoff=1)
hallmarkGSEA_I27 <- GSEA(geneList_I27, TERM2GENE = m_t2g, pvalueCutoff=1)

hallmarkGSEA_merged <- rbind(
  as.data.frame(hallmarkGSEA_Wbo2) %>% mutate(CellLine = "WBO2"),
  as.data.frame(hallmarkGSEA_I27) %>% mutate(CellLine = "I27")
)

write.csv(hallmarkGSEA_merged, "Figure_3/results/HALLMARK_merged_mG.csv")
save(hallmarkGSEA_Wbo2, hallmarkGSEA_I27, file = "Figure_3/results/HALLMARK_merged_mG.RData")

