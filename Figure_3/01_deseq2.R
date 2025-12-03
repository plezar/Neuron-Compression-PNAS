library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(AnnotationDbi)
library(metaseqR2)
source("util.R")

## ===== Human iN =====
load("Figure_3/data/hs/mtx_data.RData")

meta_data$CellLine <- factor(meta_data$CellLine, levels = c("Wbo2", "I27", "1019"))
levels(meta_data$CellLine) <- c("iN #1", "iN #2", "iN #3")
meta_data$Condition <- factor(meta_data$Condition, levels = c("control", "compress"))
levels(meta_data$Condition) <- c("Ctrl", "Comp")
meta_data$group <- paste(meta_data$CellLine, meta_data$Condition, sep=" ")
meta_data$group <- factor(meta_data$group, levels = c("iN #1 Ctrl", "iN #1 Comp", "iN #2 Ctrl", "iN #2 Comp", "iN #3 Ctrl", "iN #3 Comp"))

### ===== Removing due to very low seq. depth =====
samples_to_rm <- !(rownames(meta_data) %in% c("I27_compress_05-20_Seq_07-25_S22", "I27_control_05-20_Seq_07-25_S21"))

### ===== Downsampling to remove artefacts due to seq depth =====
set.seed(44)
count_mtx_downsampled <- downsampleCounts(count_mtx[,samples_to_rm])

### ===== All samples =====
dds_downsampled <- DESeqDataSetFromMatrix(countData = count_mtx_downsampled,
                                          colData = meta_data[samples_to_rm,],
                                          design =  ~ Condition + CellLine + Batch)
smallestGroupSize <- 4
keep <- rowSums(counts(dds_downsampled) >= 10) >= smallestGroupSize
dds_downsampled <- dds_downsampled[keep,]
dds_downsampled <- estimateSizeFactors(dds_downsampled)
dds_downsampled <- DESeq(dds_downsampled)
saveRDS(dds_downsampled, "Figure_3/results/dds_downsampled_iN.rds")

### ===== iN #1 =====
dds <- DESeqDataSetFromMatrix(countData = count_mtx[,meta_data$CellLine=="iN #1"],
                              colData = meta_data[meta_data$CellLine=="iN #1",],
                              design =  ~ Condition + Batch)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)

res_Wbo2 <- results(dds, alpha = 0.05, contrast=c("Condition", "Comp", "Ctrl"), cooksCutoff=FALSE, independentFiltering=FALSE) %>% as.data.frame() %>% drop_na() %>% arrange(padj) %>% mutate(DE = call_DE_genes(log2FoldChange, padj, 0, 0.05))

### ===== iN #2 =====
dds <- DESeqDataSetFromMatrix(countData = count_mtx[,meta_data$CellLine=="iN #2"],
                              colData = meta_data[meta_data$CellLine=="iN #2",],
                              design =  ~ Condition + Batch)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)

res_I27 <- results(dds, alpha = 0.05, contrast=c("Condition", "Comp", "Ctrl"), cooksCutoff=FALSE, independentFiltering=FALSE) %>% as.data.frame() %>% drop_na() %>% arrange(padj) %>% mutate(DE = call_DE_genes(log2FoldChange, padj, 0, 0.05))

### ===== iN #3 =====
dds <- DESeqDataSetFromMatrix(countData = count_mtx[,meta_data$CellLine=="iN #3"],
                              colData = meta_data[meta_data$CellLine=="iN #3",],
                              design =  ~ Condition + Batch)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)

res_1019 <- results(dds, alpha = 0.05, contrast=c("Condition", "Comp", "Ctrl"), cooksCutoff=FALSE, independentFiltering=FALSE) %>% as.data.frame() %>% drop_na() %>% arrange(padj) %>% mutate(DE = call_DE_genes(log2FoldChange, padj, 0, 0.05))

write.csv(res_Wbo2, "Figure_3/results/iN1_neuron.csv")
write.csv(res_I27, "Figure_3/results/iN2_neuron.csv")
write.csv(res_1019, "Figure_3/results/iN3_neuron.csv")


## ===== mGlia =====
load("Figure_3/data/mmus/mtx_data.RData")

meta_data$CellLine <- factor(meta_data$CellLine, levels = c("Wbo2", "I27", "1019"))
levels(meta_data$CellLine) <- c("mG #1", "mG #2", "mG #3")
meta_data$Condition <- factor(meta_data$Condition, levels = c("control", "compress"))
levels(meta_data$Condition) <- c("Ctrl", "Comp")
meta_data$group <- paste(meta_data$CellLine, meta_data$Condition, sep=" ")
meta_data$group <- factor(meta_data$group, levels = c("mG #1 Ctrl", "mG #1 Comp", "mG #2 Ctrl", "mG #2 Comp", "mG #3 Ctrl", "mG #3 Comp"))

### ===== Downsampling to remove artefacts due to seq depth =====
set.seed(44)
count_mtx_downsampled <- downsampleCounts(count_mtx)

### ===== All samples =====
dds_downsampled <- DESeqDataSetFromMatrix(countData = count_mtx_downsampled,
                                          colData = meta_data,
                                          design =  ~ Condition + CellLine + Batch)
smallestGroupSize <- 4
keep <- rowSums(counts(dds_downsampled) >= 10) >= smallestGroupSize
dds_downsampled <- dds_downsampled[keep,]
dds_downsampled <- estimateSizeFactors(dds_downsampled)
dds_downsampled <- DESeq(dds_downsampled)
saveRDS(dds_downsampled, "Figure_3/results/dds_downsampled_mG.rds")

### ===== mG #1 =====
dds <- DESeqDataSetFromMatrix(countData = count_mtx[,meta_data$CellLine=="mG #1"],
                              colData = meta_data[meta_data$CellLine=="mG #1",],
                              design =  ~ Condition + Batch)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)

res_Wbo2 <- results(dds, alpha = 0.05, contrast=c("Condition", "Comp", "Ctrl"), cooksCutoff=FALSE, independentFiltering=FALSE) %>% as.data.frame() %>% drop_na() %>% arrange(padj) %>% mutate(DE = call_DE_genes(log2FoldChange, padj, 0, 0.05))

### ===== mG #2 =====
dds <- DESeqDataSetFromMatrix(countData = count_mtx[,meta_data$CellLine=="mG #2"],
                              colData = meta_data[meta_data$CellLine=="mG #2",],
                              design =  ~ Condition + Batch)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)

res_I27 <- results(dds, alpha = 0.05, contrast=c("Condition", "Comp", "Ctrl"), cooksCutoff=FALSE, independentFiltering=FALSE) %>% as.data.frame() %>% drop_na() %>% arrange(padj) %>% mutate(DE = call_DE_genes(log2FoldChange, padj, 0, 0.05))

### ===== mG #3 =====
dds <- DESeqDataSetFromMatrix(countData = count_mtx[,meta_data$CellLine=="mG #3"],
                              colData = meta_data[meta_data$CellLine=="mG #3",],
                              design =  ~ Condition + Batch)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
vsd <- vst(dds)

res_1019 <- results(dds, alpha = 0.05, contrast=c("Condition", "Comp", "Ctrl"), cooksCutoff=FALSE, independentFiltering=FALSE) %>% as.data.frame() %>% drop_na() %>% arrange(padj) %>% mutate(DE = call_DE_genes(log2FoldChange, padj, 0, 0.05))

write.csv(res_Wbo2, "Figure_3/results/iN1_glia.csv")
write.csv(res_I27, "Figure_3/results/iN2_glia.csv")
write.csv(res_1019, "Figure_3/results/iN3_glia.csv")