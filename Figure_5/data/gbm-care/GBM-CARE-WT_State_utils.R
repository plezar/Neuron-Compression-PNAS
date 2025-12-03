
#
# This function generates a null distribution for state/cell-type classification by sampling cells and shuffling the
# expression values, thus generating random cells. These random cells are then scored for the supplied set of signatures.
# In this approach the cells are scored across all samples.
#
generate_null_dist <- function(umi_data_list, md, sigs, genes_subset = NULL, n_iter = 20, n_cells = 10000, verbose = F) {
  
  if(verbose == T)
    print("Generating NULL distribution for state/cell-type classification")
  
  start_time <- Sys.time()
  permuted_data <- lapply(1:n_iter, function(i) {
    
    if(verbose == T)
      print(paste0("Iteration ", i))
    
    cellids <- md %>% sample_n(n_cells) %>% pull(CellID)
    
    if(!is.null(genes_subset))
      genes <- genes_subset
    else
      genes <- rownames(umi_data_list[[1]])
    
    if(verbose == T)
      print("Aggregating UMI data")
    
    umi_data <- lapply(names(umi_data_list), function(sname) {
      d <- md %>% filter(Sample == sname, CellID %in% cellids)
      if(nrow(d) < 2)
        return(NULL)
      m <- umi_data_list[[sname]]
      m <- m[genes, d$CellID]
      m
    })
    umi_data <- do.call(cbind, umi_data)
    gc()
    
    cellids <- cellids[cellids %in% colnames(umi_data)]
    
    m <- umi_data[, cellids]
    dim(m)
    
    if(verbose == T)
      print("Permuting matrix")
    
    m <- t(apply(m, 1, gtools::permute))
    colnames(m) <- cellids; rownames(m) <- genes
    
    if(verbose == T)
      print("Creating object")
    
    tmp <- CreateSeuratObject(counts = m,
                              project = paste0("Perm", i),
                              min.cells = 3, min.features = 200)
    
    if(verbose == T)
      print("Normalizing")
    
    tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
    
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
    
    if(verbose == T)
      print("Scaling")
    
    all.genes <- rownames(tmp)
    tmp <- ScaleData(tmp, features = all.genes)
    
    if(verbose == T)
      print("Scoring")
    
    tmp <- AddModuleScore(object = tmp, features = sigs)
    
    d <- as_tibble(tmp@meta.data)
    colnames(d)[grep("Cluster", colnames(d))] <- names(sigs)
    
    d[, -c(1:3)]
  })
  end_time <- Sys.time()
  end_time - start_time
  
  permuted_data <- do.call(rbind, permuted_data)
  
  return(permuted_data)
}

#
# This function generates a null distribution for state/cell-type classification by sampling cells and shuffling the
# expression values, thus generating random cells. These random cells are then scored for the supplied set of signatures.
# In this approach the cells are scored within each samples.
#
score_within_samples <- function(umi_data_list, md, sigs) {
  
  samples <- unique(md$Sample)
  
  res_all <- lapply(samples, function(sname) {
    
    print(sname)
    
    d <- md %>% filter(Sample == sname)
    
    ud <- umi_data_list[[sname]][, d$CellID]
    
    npcs <- 100
    
    if(nrow(d) < npcs) {
      npcs <- nrow(d) - 1
      print(paste0("npcs set to ", npcs))    
    }
    
    tmp <- new_seurat_object(m = ud, md = d, project = sname, min.cells = 0, min.features = 0, npcs = npcs,
                             run_pca = F, run_clustering = F, run_umap = F, verbose = T)
    
    tmp <- AddModuleScore(object = tmp, features = sigs)  
    
    scores <- tmp@meta.data[, grep("Cluster", colnames(tmp@meta.data))]
    colnames(scores) <- names(sigs)
    
    scores <- as_tibble(scores, rownames = "CellID")
    
    res <- d %>%
      select(CellID, Patient, Timepoint, Sample) %>%
      left_join(scores, by = "CellID")
    
    return(res)
  })
  res_all <- do.call(rbind, res_all)
  dim(res_all)
  
  res_all <- as_tibble(res_all)
  
  return(res_all)
}
