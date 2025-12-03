library(Seurat)
library(SeuratDisk)
library(patchwork)

source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_analysis_utils.R")
source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_CNA_utils.R")
#source("GBM-CARE-WT/R/GBM-CARE-WT_NMF.R")
source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_State_utils.R")

ANALYSIS_ROOT <- "/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/"
DATA_ROOT <- paste0(ANALYSIS_ROOT, "/data/")
TABLES_ROOT <- paste0(ANALYSIS_ROOT, "/tables/")
FIGURES_ROOT  <- paste0(ANALYSIS_ROOT, "/figures/")
SDS_SAMPLE_PATH <- paste0(DATA_ROOT, "/sample_sds/")
NMSDS_SAMPLE_PATH <- paste0(DATA_ROOT, "/sample_nmsds/")

meta_data <- readRDS(paste0(DATA_ROOT, "meta_data_qc.RDS"))
umi_data_all <- readRDS(paste0(DATA_ROOT, "umi_data_list.RDS"))

# Generate UMI matrix for seurat dataset ------------------------------------------------------------------------------------------

umi_data <- lapply(unique(meta_data$Sample), function(sname) {
  m <- umi_data_all[[sname]]
  m <- m[, meta_data$CellID[meta_data$Sample == sname]]
  m
})
umi_data <- do.call(cbind, umi_data)

message("Dimensions of the merged dataset:")
print(dim(umi_data))
rm(umi_data_all)
gc()

saveRDS(object = umi_data, file = paste0(DATA_ROOT, "umi_data.RDS"))


# Create a top-level Seurat object ------------------------------------------------------------------------------------------

sds <- CreateSeuratObject(counts = umi_data[, meta_data$CellID],
                            project = "GBM-CARE_integrated_analysis",
                            min.cells = 3, min.features = 0, 
                            meta.data = meta_data)
gc()

sds <- NormalizeData(sds)
gc()

sds <- FindVariableFeatures(sds, selection.method = "vst", nfeatures = 5000)
gc()

sds <- ScaleData(sds, features = rownames(sds))
gc()

sds <- RunPCA(sds, features = VariableFeatures(object = sds))
gc()

ep <- ElbowPlot(sds)
ggsave(ep, path=FIGURES_ROOT, filename="ElbowPlot.pdf", width=5, height=5, dpi=700)

sds <- FindNeighbors(sds, dims = 1:15)
gc()

sds <- FindClusters(sds, resolution = .5)
gc()

sds <- RunUMAP(sds, dims = 1:25)
gc()

dim(sds)

p1 <- DimPlot(sds, reduction = "umap", group.by = "Patient")
p2 <- DimPlot(sds, reduction = "umap", group.by = "Timepoint")
p3 <- DimPlot(sds, reduction = "umap", group.by = "Lab_processed")
p4 <- DimPlot(sds, reduction = "umap", group.by = "Tissue_Source")
# TO-DO add FeaturePlots

p <- p1 + p2 + p3 + p4
ggsave(p, path=FIGURES_ROOT, filename="UMAP.pdf", width=15, height=15, dpi=700)

saveRDS(sds, paste0(DATA_ROOT, "sds.RDS"))

umap_res <- Reductions(sds, "umap")
meta_data$UMAP1 <- umap_res@cell.embeddings[meta_data$CellID, 1]
meta_data$UMAP2 <- umap_res@cell.embeddings[meta_data$CellID, 2]
gc()

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data_qc_umap.RDS"))


# Create sample-level Seurat objects ------------------------------------------------------------------------------------------

samples <- unique(meta_data$Sample)
length(samples)

start_time <- Sys.time()
for(i in 1:length(samples)) {
  
  gc()
  
  sname <- samples[i]
  
  message(paste0("******************************* ", sname, " - start, i = ", i, " *******************************"))
  
  cellids <- rownames(meta_data)[meta_data$Sample == sname]
  
  tmp <- new_seurat_object(m = umi_data[, cellids], md = meta_data[cellids, ],
                           project = paste0("GBM-Longitudinal_", sname), nfeatures = 3000, verbose = T)
  dim(tmp)
  
  
  message("Saving")
  st <- Sys.time()
  saveRDS(object = tmp, file = paste0(SDS_SAMPLE_PATH, "sds_", sname, ".RDS"))
  et <- Sys.time()
  print(et - st)
  
  message(paste0("******************************* ", sname, " - end, i = ", i, " *******************************"))
}
end_time <- Sys.time()
end_time - start_time