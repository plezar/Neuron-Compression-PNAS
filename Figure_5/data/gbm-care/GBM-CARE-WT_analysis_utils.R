library(BiocParallel)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Global definitions
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

celltype_color_vec <- c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Macrophage" = "#80B1D3", # 4, 7, 5
                        "Astrocyte" = "purple", "Excitatory neuron" = "#BC80BD", "Inhibitory neuron" = "#FFED6F", "Other neuron" = "#CCEBC5", # 6, 10, 12
                        "Endothel" = "#FCCDE5", "Pericyte" = "#FFFFB3", "Tcell" = "#8DD3C7", "Bcell" = "#BEBADA", "OPC" = "#FDB462", # 8, 2, 1, 3, 11
                        "Other normal" = "#D9D9D9", "Unresolved" = "grey") # 10

celltype_color_vec_reduced <- c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "TAM" = "#80B1D3",
                                "Astrocyte" = "#BEBADA", "Excitatory neuron" = "#BC80BD", "Inhibitory neuron" = "#FFED6F",
                                "Endothel" = "#FCCDE5", "Pericyte" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "OPC" = "#CCEBC5",
                                "Other" = "grey") # 10

timepoint_color_vec <- c("T1" = "#66C2A5", "T2" = "#FC8D62", "T3" = "#8DA0CB")

mgmt_color_vec <- c("MET" = "red", "UM" = "black", "NOS" = "grey")

scp_color_vec <- c("SCP-ECM" = "#41AB5D", "SCP-Neuronal" = "#807DBA", "SCP-Glial" = "#EF3B2C", "Mixed" = "grey")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Functions
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

`%ni%` <- Negate(`%in%`)

mem_usage <- function() {
  sort(sapply(ls(),function(x){object.size(get(x))})) 
}

load_sample <- function(sname, path) {
  readRDS(paste0(path, "sds_", sname, ".RDS"))
}

plapply <- function(X, FUN, ..., pack = F, profile = F, n_workers = 8) {
  
  start_time <- Sys.time()
  res <- bplapply(X = X, FUN = FUN, ... = ..., BPPARAM = MulticoreParam(workers = n_workers))
  end_time <- Sys.time()
  end_time - start_time
  print(end_time - start_time)
  
  if(isTRUE(pack))
    res <- unlist(lapply(res, function(x) x), recursive = F)
  
  res
}

umi2upm <- function(m) {
  count_sum <- colSums(m)
  upm_data <- (t(t(m)/count_sum)) * 1e+06
  upm_data
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

select_genes <- function(umi_data_list, # named list of UMI count matrices, one per sample
                         md, # metadata data.frame with at least columns Sample and CellID
                         exp_th = 5, # log₂(mean UPM + 1) threshold for “highly expressed” in step 1
                         freq_th1 = 25, # minimum number of samples in which a gene must exceed exp_th
                         umi_th = 10, # UMI‐count threshold per cell in step 2
                         cells_th = 10, # minimum cells per sample with UMI > umi_th
                         freq_th2 = 10, # minimum number of samples in which a gene must exceed cells_th
                         verbose = F,
                         plot = F) {
  
  if(verbose == T)
    print("Detecting HE genes - step 1")

  ### Part 1
  # find genes with average expression above exp_th
  he_genes <- lapply(X = unique(md$Sample),
                     FUN = function(sname) {
                       m <- umi_data_list[[sname]]
                       m <- m[, md$CellID[md$Sample == sname]]
                       m <- umi2upm(m)
                       heg <- log2(rowMeans(m) + 1)
                       rownames(m)[heg > exp_th]
                     })
  he_genes <- unname(unlist(he_genes, recursive = F))
  
  # in how many samples each gene eas HE
  freq <- as.data.frame(table(he_genes))
  
  if(plot == T)
    print(gghistogram(data = freq, x = "Freq", bins = 25, fill = "slateblue", alpha = .75, color = "black") +
            geom_vline(xintercept = freq_th1, linetype = "dashed", color = "black", size = 1))
  
  # keep only genes that are HE in at least freq_th1 samples, 5 by default
  qc_genes <- as.character(freq$he_genes[freq$Freq >= freq_th1])
  length(qc_genes)
  
  if(verbose == T)
    print("Detecting HE genes - step 2")
  
  ### Part 2
  # find genes that are detected robustly (>= umi_th) in at least cells_th per sample
  ext_genes <- lapply(X = unique(md$Sample),
                      FUN = function(sname) {
                        m <- umi_data_list[[sname]]
                        m <- m[, md$CellID[md$Sample == sname]]
                        m <- m > umi_th
                        rowSums(m)
                      })
  he_genes2 <- unname(unlist(lapply(ext_genes, function(eg) names(eg)[eg >= cells_th]), recursive = F))
  table(he_genes2) %>% sort(decreasing = T)
  
  freq <- as.data.frame(table(he_genes2))
  
  if(plot == T)
    print(gghistogram(data = freq, x = "Freq", bins = 25, fill = "slateblue", alpha = .75, color = "black") +
            geom_vline(xintercept = freq_th2, linetype = "dashed", color = "black", size = 1))
  
  # Do so in >= freq_th2 samples
  qc_genes2 <- as.character(freq$he_genes2[freq$Freq >= freq_th2])
  
  # merge
  genes <- c(qc_genes, qc_genes2)
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(umi_data_list[[1]])]
  length(genes)
  
  if(verbose == T)
    print(paste0("Detected ", length(genes), " HE genes"))
  
  return(genes)
}

new_seurat_object <- function(m, md, project, min.cells = 3, min.features = 0, cell_subset = NULL, gene_subset = NULL,
                              scale.factor = 10000, nfeatures = 5000, npcs = 100, resolution = 0.5,
                              run_pca = T, run_clustering = T, run_umap = T, metric = "cosine", verbose = F) {
  
  stopifnot(!is.null(rownames(md)))
  
  if(verbose == T)
    print(paste0("Generating new seurat object for project ", project))
  
  if(!is.null(cell_subset)) {
    
    m <- m[, cell_subset]
    md <- md[, cell_subset]
  }
  
  if(!is.null(gene_subset))
    m <- m[gene_subset, ]
  
  tmp <- CreateSeuratObject(counts = m,
                            project = project,
                            min.cells = min.cells, min.features = min.features, 
                            meta.data = md)
  dim(tmp)
  gc()
  
  tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = scale.factor)
  gc()
  
  tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = nfeatures)
  gc()
  
  tmp <- ScaleData(tmp, features = rownames(tmp))
  gc()
  
  if(run_pca) {
    tmp <- RunPCA(tmp, features = VariableFeatures(object = tmp), npcs = npcs)
    gc()
  }
  
  # Elbow plot? Not sure why 100 PCs are needed
  # options(future.globals.maxSize = 4000 * 1024^2)
  
  if(run_clustering) {
    tmp <- FindNeighbors(tmp, dims = 1:npcs)
    tmp <- FindClusters(tmp, resolution = resolution)
    gc()
  }
  
  if(run_umap) {
    n_neighbors <- 30L
    if(n_neighbors < nrow(tmp))
      n_neighbors <- 5
    
    tmp <- RunUMAP(tmp, dims = 1:npcs, n.neighbors = n_neighbors, metric = metric)
    gc()
  }
  
  return(tmp)
}

doublet_detection_null_dist <- function(m, cluster_means, frac_sampled = .1, n_iter = 100, verbose = F) {
  
  n_sampled <- ceiling(ncol(m) * frac_sampled)
  
  if(verbose == T)
    print(paste0("Sampling ", n_sampled, " cells ", n_iter, " times"))
  
  null_dist <- lapply(1:n_iter, function(i) {
    
    print(paste0("Iteration ", i))
    
    sampled_cells <- sample(colnames(m), n_sampled)
    
    m_perm <- m[, sampled_cells]
    
    m_perm <- t(apply(m_perm, 1, gtools::permute))
    colnames(m_perm) <- sampled_cells; rownames(m_perm) <- rownames(m)
    
    rm <- rowMeans(m)
    
    cr <- sapply(1:ncol(m_perm), function(j) {
      m_perm[, j] #- rm
    })
    dim(cr)
    colnames(cr) <- sampled_cells
    
    cc <- sapply(1:ncol(m_perm), function(j) {
      cor(cr[, j], cluster_means)
    })
    colnames(cc) <- sampled_cells
    
    mc <- apply(cc, 2, function(x) max(x))
    names(mc) <- paste0(names(mc), "_", i)
    
    return(mc)
  })
  null_dist <- unname(unlist(null_dist, recursive = F))
  
  return(null_dist)  
}

patient_pairs <- function(data) {
  
  pt_pairs <- data %>%
    dplyr::filter(Timepoint != "T3") %>%
    group_by(Patient) %>%
    filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
    group_by(Patient, Timepoint, Sample, IDHstatus) %>%
    summarise(n = n()) %>%
    ungroup()
  
  genotype <- lapply(unique(pt_pairs$Patient), function(pt) {
    d <- pt_pairs %>% filter(Patient == pt)
    res <- "WT"
    if(d$IDHstatus[d$Timepoint == "T1"] == "NOS")
      res <- d$IDHstatus[d$Timepoint == "T2"]
    else
      res <- d$IDHstatus[d$Timepoint == "T1"]
    return(tibble(Patient = pt, Genotype = res))
  })
  genotype <- do.call(rbind, genotype)
  
  pt_pairs <- pt_pairs %>%
    left_join(genotype, by = "Patient")
  
  return(pt_pairs)
}

plot_samples_butterfly_density <- function(data, samples, nrows = 1) {
  
  d <- data %>% filter(Sample %in% samples)
  
  d <- d %>% arrange(Timepoint)
  
  d$ID <- factor(d$ID, unique(d$ID))
  
  ggplot(data = d, mapping = aes(x = Dx, y = Dy)) +
    facet_wrap(facets = ~ID, nrow = nrows) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(name = "Number of cells", option = "plasma",
                         breaks = seq(0, 15, 5), labels = c("0", "5", "10", "15"), limits = c(0, 15)) +
    # scale_fill_viridis_c(name = "Number of cells", option = "plasma",
    #                      breaks = seq(0, 30, 5), labels = c("0", "", "10", "", "20", "", "30"), limits = c(0, 30)) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
    theme_pubr() +
    xlab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
    ylab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
    theme(axis.title = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.ticks.x = element_line(), axis.ticks.y = element_line(),
          legend.title = element_text(size = 30), legend.text = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 20)) +
    theme(panel.spacing = unit(.0, "lines"),
          panel.border = element_rect(color = "black", fill = NA), 
          strip.background = element_rect(color = "black", fill = "grey90"))
}

plot_samples_butterfly <- function(data, samples, nrows = 1) {
  
  d <- data %>% filter(Sample %in% samples)
  
  d <- d %>% arrange(Timepoint)
  
  d$ID <- factor(d$ID, unique(d$ID))
  
  ggplot(data = d, mapping = aes(x = Dx, y = Dy)) +
    # facet_grid(cols = vars(Patient_factor), rows = vars(Timepoint)) +
    # facet_grid(cols = vars(ID)) +
    facet_wrap(facets = ~ID, nrow = nrows) +
    geom_point(mapping = aes(color = State_factor, fill = State_factor, alpha = State_factor), shape = 21, color = "black", size = 5,
               show.legend = c("color" = TRUE, "fill" = TRUE, "alpha" = FALSE)) +
    # scale_fill_manual(name = "State", values = c("AC-like" = "#F8766D","MES-like" = "#A3A500", "Synaptic" = "darkgreen",
    #                                              "NPC-like" = "#00BF7D", "OPC-like" = "#00B0F6", Hybrid = "#E76BF3", Undifferentiated = "black")) +
    scale_fill_discrete(name = "State") +
    scale_alpha_manual(values = c(AC = 1, MES = 1, NPC = 1, OPC = 1, Hybrid = .5, Undifferentiated = .25)) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
    theme_pubr() +
    xlab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
    ylab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
    theme(axis.title = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.ticks.x = element_line(), axis.ticks.y = element_line(),
          legend.title = element_text(size = 30), legend.text = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 20)) +
    theme(panel.spacing = unit(.0, "lines"),
          panel.border = element_rect(color = "black", fill = NA), 
          strip.background = element_rect(color = "black", fill = "grey90"))
}

plot_patients_butterfly <- function(data, patients) {
  
  ggplot(data = data %>% filter(Patient %in% patients), mapping = aes(x = Dx, y = Dy)) +
    facet_grid(cols = vars(Patient_factor), rows = vars(Timepoint)) +
    # facet_grid(cols = vars(ID)) +
    geom_point(mapping = aes(color = State, fill = State, alpha = State), shape = 21, color = "black", size = 5,
               show.legend = c("color" = TRUE, "fill" = TRUE, "alpha" = FALSE)) +
    scale_fill_manual(name = "State", values = c("AC-like" = "#F8766D","MES-like" = "#A3A500",
                                                 "NPC-like" = "#00BF7D", "OPC-like" = "#00B0F6", Hybrid = "#E76BF3", Undifferentiated = "black")) +
    scale_alpha_manual(values = c(AC = 1, MES = 1, NPC = 1, OPC = 1, Hybrid = .5, Undifferentiated = .25)) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 2) +
    theme_pubr() +
    xlab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
    ylab("Relative meta-module score\n[log2(|SC1 - SC2| + 1)]") +
    theme(axis.title = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.ticks.x = element_line(), axis.ticks.y = element_line(),
          legend.title = element_text(size = 30), legend.text = element_text(size = 20),
          plot.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 20)) +
    theme(panel.spacing = unit(.0, "lines"),
          panel.border = element_rect(color = "black", fill = NA), 
          strip.background = element_rect(color = "black", fill = "grey90"))
}

sample2patient <- function(data, samples) {
  data %>% filter(Sample %in% samples) %>% pull(Patient) %>% unique()
}

theme_gbm_pvsr <- function(theme = theme_pubr(),
                           axis.title = element_text(size = 30), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
                           axis.ticks.x = element_line(), axis.ticks.y = element_line(),
                           legend.title = element_text(size = 30), legend.text = element_text(size = 20), legend.position = "top",
                           plot.title = element_text(size = 20, face = "bold"), strip.text = element_text(size = 20),
                           panel.spacing = unit(.0, "lines"),  panel.border = element_rect(color = "black", fill = NA), 
                           strip.background = element_rect(color = "black", fill = "grey90")) {
  
  res <- theme +
    theme(axis.title = axis.title,
          axis.text.x = axis.text.x,
          axis.text.y = axis.text.y,
          axis.ticks.x = axis.ticks.x, axis.ticks.y = axis.ticks.y,
          legend.title = legend.title, legend.text = legend.text, legend.position = legend.position,
          plot.title = plot.title,
          strip.text = strip.text) +
    theme(panel.spacing = panel.spacing,
          panel.border = panel.border, 
          strip.background = strip.background)
  
  return(res)
}

add_glioma_sigs <- function(pathways) {
  pathways$IDH_O_AC <- scandal::SCANDAL_IDH_O_AC_MARKERS
  pathways$IDH_O_OC <- scandal::SCANDAL_IDH_O_OC_MARKERS
  pathways$IDH_O_NPC <- scandal::SCANDAL_IDH_O_STEMNESS_MARKERS
  pathways$IDH_A_AC <- scandal::SCANDAL_IDH_A_AC_MARKERS
  pathways$IDH_A_OC <- scandal::SCANDAL_IDH_A_OC_MARKERS
  pathways$IDH_A_NPC <- scandal::SCANDAL_IDH_A_STEMNESS_MARKERS
  pathways$GLIOMA_G1S <- unique(c(scandal::SCANDAL_IDH_O_G1S_MARKERS, scandal::SCANDAL_IDH_A_G1S_MARKERS, scalop::Signatures_GBM$G1S))
  pathways$GLIOMA_G2M <- unique(c(scandal::SCANDAL_IDH_O_G2M_MARKERS, scandal::SCANDAL_IDH_A_G2M_MARKERS, scalop::Signatures_GBM$G2M))
  pathways$GBM_AC <- scalop::Signatures_GBM$AC
  pathways$GBM_MES1 <- scalop::Signatures_GBM$MES1
  pathways$GBM_MES2 <- scalop::Signatures_GBM$MES2
  pathways$GBM_NPC1 <- scalop::Signatures_GBM$NPC1
  pathways$GBM_NPC2 <- scalop::Signatures_GBM$NPC2
  pathways$GBM_OPC <- scalop::Signatures_GBM$OPC
  
  return(pathways)
}

add_normal_brain_sigs <- function(pathways) {
  
  sigs <- load(file = "data/BrainSignatures.rda")
  
  pathways$`AC-camp17` <- BrainSignatures$`AC-camp17`
  pathways$`AC-habib16` <- BrainSignatures$`AC-habib16`
  pathways$`AC-nowak17` <- BrainSignatures$`AC-nowak17`
  pathways$`AC-schaum18` <- BrainSignatures$`AC-schaum18`
  pathways$`AC.FB-velm19` <- BrainSignatures$`AC.FB-velm19`
  pathways$`AC.PP-velm19` <- BrainSignatures$`AC.PP-velm19`
  
  pathways$`OC.1-camp17` <- BrainSignatures$`OC.1-camp17`
  pathways$`OC.2-camp17` <- BrainSignatures$`OC.2-camp17`
  pathways$`OC.3-camp17` <- BrainSignatures$`OC.3-camp17`
  pathways$`OC-habib16` <- BrainSignatures$`OC-habib16`
  pathways$`OC-schaum18` <- BrainSignatures$`OC-schaum18`
  pathways$`OC-velm19` <- BrainSignatures$`OC-velm19`
  
  pathways$`OPC-camp17` <- BrainSignatures$`OPC-camp17`
  pathways$`OPC-habib16` <- BrainSignatures$`OPC-habib16`
  pathways$`OPC-nowak17` <- BrainSignatures$`OPC-nowak17`
  pathways$`OPC-poli19` <- BrainSignatures$`OPC-poli19`
  pathways$`OPC-schaum18` <- BrainSignatures$`OPC-schaum18`
  pathways$`OPC-velm19` <- BrainSignatures$`OPC-velm19`
  
  return(pathways)
}

load_c8_pathways <- function(path = "GBM-CARE-WT/data/msigdb/") {
  
  return(fgsea::gmtPathways(paste0(path, "c8.all.v7.4.symbols.gmt")))
}

default_pathways <- function(path = "/data/msigdb/") {
  
  load(file = "~/Datasets/BrainSignatures.rda")
  
  pathways <- c()
  
  pathways <- c(pathways, BrainSignatures)
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/h.all.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c5.go.bp.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c5.go.cc.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c5.go.mf.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c2.cp.reactome.v7.4.symbols.gmt"))
  
  # pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c8.all.v7.4.symbols.gmt"))
  
  pathways <- add_glioma_sigs(pathways)
  
  pathways <- c(pathways, list(Synaptic = readRDS("integrated_analysis3/data/syn_sig.RDS")))
  
  # pathways <- c(pathways, read.csv("integrated_analysis4/MP_named.csv", header = T, stringsAsFactors = F) %>% as.list())
  
  return(pathways)  
}

functional_pathways <- function() {
  
  pathways <- c()
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/h.all.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c5.go.bp.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c5.go.cc.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c5.go.mf.v7.4.symbols.gmt"))
  
  pathways <- c(pathways, fgsea::gmtPathways("data/msigdb/c2.cp.reactome.v7.4.symbols.gmt"))
  
  return(pathways)  
}

normal_brain_signatures <- function() {
  
  load(file = "~/Datasets/BrainSignatures.rda")
  
  pathways <- c()
  
  pathways <- c(pathways, BrainSignatures)
  
  return(pathways)
}

as_bulk <- function(m, groups) {
  m_bulk <- lapply(groups, function(g) rowMeans(m[, g]))
  m_bulk <- do.call(cbind, m_bulk)
  rownames(m_bulk) <- rownames(m)
  colnames(m_bulk) <- names(groups)
  
  return(m_bulk)
}

umi_as_bulk <- function(umi_data, groups, genes_subset = NULL, to_upm = F) {
  
  m_bulk <- lapply(names(groups), function(g) {
    cellids <- groups[[g]]
    m <- umi_data[[g]]
    if(!is.null(genes_subset))
      m <- m[genes_subset, ]
    rowMeans(m[, cellids])
  })
  m_bulk <- do.call(cbind, m_bulk)
  
  if(to_upm == T)
    m_bulk <- umi2upm(m_bulk)
  
  # rownames(m_nulk) <- rownames(m)
  colnames(m_bulk) <- names(groups)
  
  return(m_bulk)
}

get_sname <- function(x, split = "_") {
  sapply(strsplit(x, split = split), function(y) y[[1]][1])
}

shannon_index <- function(x) {
  -sum(x[x > 0] * log(x[x > 0]))
}

derive_markers <- function(pathways, ct, n_rep = 5, max_genes = 75) {
  res <- unname(unlist(pathways[grep(ct, names(pathways))], recursive = F))
  res <- table(res) %>% sort(decreasing = T)
  res <- res[res >= n_rep]
  length(res)
  res <- res %>% head(n = max_genes) %>% names()
  return(res)
}

load_brain_signatures <- function() {
  
  load(file = "~/Datasets/BrainSignatures.rda")
  
  sigs_list <- list(Macrophage = derive_markers(load_c8_pathways(), "MACROPHAGE|MICROGLIA"),
                    Tcell = unique(c(derive_markers(load_c8_pathways(), "_T_CELL|_NK"),
                                     "CD2", "SKAP1", "ITK", "ICOS", "CD247", "CD96", "GZMA", "CD3E", "GZMK", "CD6", "CD3G", "IL7R", "CD3D", "TRAT1",
                                     "ZAP70", "KLRC1", "IL2RB", "TRBC1", "TRBC2", "CCL5", "CD8A")),
                    Bcell = unique(c(derive_markers(load_c8_pathways(), "_B_CELL|_PLASMA", n_rep = 3),
                                     "IGKC", "IGHM", "FCRL1", "FCRL5" , "MS4A1", "BLK", "IGHG1", "IGHG4" ,"IGHA1", "CLNK", "IGHGP", "IGHG3", "JCHAIN",
                                     "IGLC2", "IDO1", "BANK1", "IGHG2", "CXCL10", "MS4A2", "IGLC3", "BTLA", "IDO2")),
                    Endothel = derive_markers(load_c8_pathways(), "ENDOTHEL"),
                    Pericyte = derive_markers(load_c8_pathways(), "PERICYTE", n_rep = 2),
                    Oligodendrocyte = c("MOBP", "OPALIN", "CLDN11", "PLP1", "TF", "MBP", "MOG", "MAG", "EDIL3", "ST18", "CNTN2", "ANLN", "CNTN2",
                                        "DOCK5", "TMEM144", "ENPP2", "SPOCK3", "PLEKHH1", "PIP4K2A", "BCAS1", "CERCAM",
                                        "CNDP1", "CNTNAP4", "NKAIN2", "FRMD4B", "CLMN", "ZNF536", "AK5", "CTNNA3",
                                        "ELMO1", "RNF220", "PEX5L", "SLC24A2", "PCSK6", "ST6GALNAC3", "IL1RAPL1", "MAP7",
                                        "PLCL1", "UNC5C", "SLC5A11", "KLHL32", "UGT8", "CDK18", "ANO4", "FAM107B",
                                        "GRM3", "MYRF", "TMEFF2", "FOLH1", "HHIP", "KCNH8", "PLD1"))
  lengths(sigs_list)
  
  BrainSignatures <- c(BrainSignatures, sigs_list)
  
  return(BrainSignatures)
}

run_scDblFinder <- function(sample_path, md, samples, verbose = F) {
  
  require(scDblFinder)
  
  start_time <- Sys.time()
  res <- lapply(samples, function(sname) {
    
    if(verbose == T)
      print(paste0("Loading sample ", sname))
    
    tmp <- readRDS(paste0(sample_path, "sds_", sname, ".RDS"))
    dim(tmp)
    
    d <- md %>% filter(Sample == sname)
    dim(d)
    
    celltypes <- setNames(d$CellType, d$CellID)
    celltypes <- celltypes[colnames(tmp)]
    length(celltypes)
    
    Idents(tmp) <- celltypes[colnames(tmp)]
    
    if(verbose == T)
      print("Converting to SingleCellExperiment")
    
    tmp <- Seurat::as.SingleCellExperiment(tmp)
    
    if(verbose == T)
      print("Running scDblFinder")
    
    sc_dbl_f_res <- scDblFinder(sce = tmp, clusters = celltypes)
    
    sc_dbl_f_res <- sc_dbl_f_res@colData %>%
      as.data.frame() %>%
      select(CellID, Sample, starts_with("scDblFinder"))
    
    return(sc_dbl_f_res)
  })
  end_time <- Sys.time()
  end_time - start_time
  
  if(verbose == T)
    print(end_time - start_time)
  
  res <- do.call(rbind, res)
  
  return(res)  
}

dT_state <- function(md) {
  
  res <- md %>%
    group_by(ID, Patient, Timepoint, State) %>%
    summarise(n = n()) %>%
    group_by(ID, Patient, Timepoint) %>%
    mutate(N = sum(n), Freq = n / N) %>%
    ungroup()
  
  res_m <- acast(data = res, formula = ID ~ State, value.var = "Freq")
  res_m[is.na(res_m)] <- 0
  dim(res_m)
  
  res_m <- melt(data = res_m)
  colnames(res_m) <- c("ID", "State", "Freq")
  res_m <- as_tibble(res_m)
  dim(res_m)
  
  res <- res %>%
    ungroup() %>%
    dplyr::select(-State, -n, -N, -Freq) %>%
    filter(!duplicated(ID))
  
  res_m <- res_m %>%
    ungroup() %>%
    left_join(res %>%
                dplyr::select(Patient, Timepoint, ID),
              by = "ID")
  
  pts <- unique(res_m$Patient)
  states <- unique(res_m$State)
  
  t1_m <- lapply(states, function(s) res_m %>% filter(Timepoint == "T1", State == s) %>% pull(Freq))
  t1_m <- t(do.call(rbind, t1_m))
  colnames(t1_m) <- paste0("T1_ABS_State_", states)
  rownames(t1_m) <- pts
  
  t2_m <- lapply(states, function(s) t2 <- res_m %>% filter(Timepoint == "T2", State == s) %>% pull(Freq))
  t2_m <- t(do.call(rbind, t2_m))
  colnames(t2_m) <- paste0("T2_ABS_State_", states)
  rownames(t2_m) <- pts
  
  t2_vs_t1 <- lapply(1:nrow(t1_m), function(i) t2_m[i, ] - t1_m[i, ])
  t2_vs_t1 <- do.call(rbind, t2_vs_t1)
  colnames(t2_vs_t1) <- paste0("dT_State_", states)
  rownames(t2_vs_t1) <- pts
  
  res_m <- list(dT = t2_vs_t1,
                T1_ABS = t1_m,
                T2_ABS = t2_m)
  
  return(res_m)
}

dT_CC <- function(md) {
  
  t2_vs_t1_cc <- md %>%
    group_by(Patient, Timepoint, isCC) %>%
    summarise(n = n()) %>%
    group_by(Patient, Timepoint) %>%
    mutate(N = sum(n), Freq = n / N) %>%
    filter(isCC == T)
  
  t2_vs_t1_cc <- acast(data = t2_vs_t1_cc, formula = Patient ~ Timepoint, value.var = "Freq")
  t2_vs_t1_cc[is.na(t2_vs_t1_cc)] <- 0
  
  t1 <- t2_vs_t1_cc[, 1]
  t2 <- t2_vs_t1_cc[, 2]
  
  t2_vs_t1_cc_m <- t2 - t1
  t2_vs_t1_cc_m <- matrix(data = t2_vs_t1_cc_m, ncol = 1)
  colnames(t2_vs_t1_cc_m) <- "dT_CC"
  rownames(t2_vs_t1_cc_m) <- rownames(t2_vs_t1_cc)
  
  t1_abs_cc_m <- t1
  t1_abs_cc_m <- matrix(data = t1_abs_cc_m, ncol = 1)
  colnames(t1_abs_cc_m) <- "T1_ABS_CC"
  rownames(t1_abs_cc_m) <- rownames(t2_vs_t1_cc)
  
  t2_abs_cc_m <- t2
  t2_abs_cc_m <- matrix(data = t2_abs_cc_m, ncol = 1)
  colnames(t2_abs_cc_m) <- "T2_ABS_CC"
  rownames(t2_abs_cc_m) <- rownames(t2_vs_t1_cc)
  
  res_m <- list(dT = t2_vs_t1_cc_m,
                T1_ABS = t1_abs_cc_m,
                T2_ABS = t2_abs_cc_m)
  
  return(res_m)
}

dTME_celltype <- function(md, celltype_ext = F) {
  
  if(celltype_ext)
    md$CellType <- md$CellType_ext
  
  res <- md %>%
    group_by(ID, Patient, Timepoint, CellType) %>%
    summarise(n = n()) %>%
    group_by(ID, Patient, Timepoint) %>%
    mutate(N = sum(n), Freq = n / N) %>%
    ungroup()
  
  res_m <- acast(data = res, formula = ID ~ CellType, value.var = "Freq")
  res_m[is.na(res_m)] <- 0
  dim(res_m)
  
  res_m <- melt(data = res_m)
  colnames(res_m) <- c("ID", "CellType", "Freq")
  res_m <- as_tibble(res_m)
  dim(res_m)
  
  res <- res %>%
    ungroup() %>%
    dplyr::select(-CellType, -n, -N, -Freq) %>%
    filter(!duplicated(ID))
  
  res_m <- res_m %>%
    ungroup() %>%
    left_join(res %>%
                dplyr::select(Patient, Timepoint, ID),
              by = "ID")
  
  pts <- unique(res_m$Patient)
  states <- unique(res_m$CellType)
  
  t1_m <- lapply(states, function(s) res_m %>% filter(Timepoint == "T1", CellType == s) %>% pull(Freq))
  t1_m <- t(do.call(rbind, t1_m))
  colnames(t1_m) <- paste0("T1_ABS_TME_", states)
  rownames(t1_m) <- pts
  
  t2_m <- lapply(states, function(s) t2 <- res_m %>% filter(Timepoint == "T2", CellType == s) %>% pull(Freq))
  t2_m <- t(do.call(rbind, t2_m))
  colnames(t2_m) <- paste0("T2_ABS_TME_", states)
  rownames(t2_m) <- pts
  
  t2_vs_t1 <- lapply(1:nrow(t1_m), function(i) t2_m[i, ] - t1_m[i, ])
  t2_vs_t1 <- do.call(rbind, t2_vs_t1)
  colnames(t2_vs_t1) <- paste0("dTME_CellType_", states)
  rownames(t2_vs_t1) <- pts
  
  res_m <- list(dT = t2_vs_t1,
                T1_ABS = t1_m,
                T2_ABS = t2_m)
  
  return(res_m)
}

dTME_states <- function(md, states, pt_ord) {
  
  md$CellType <- md$CellType_ext
  
  res <- md %>%
    filter(CellType %in% states) %>%
    group_by(ID, Patient, Timepoint, CellType) %>%
    summarise(n = n()) %>%
    group_by(ID, Patient, Timepoint) %>%
    mutate(N = sum(n), Freq = n / N) %>%
    ungroup()
  
  res_m <- acast(data = res, formula = ID ~ CellType, value.var = "Freq")
  res_m[is.na(res_m)] <- 0
  dim(res_m)
  
  missing_ids <- setdiff(unique(md$ID), rownames(res_m))
  missing_ids_m <- matrix(data = 0, nrow = length(missing_ids), ncol = ncol(res_m),
                          dimnames = list(missing_ids, colnames(res_m)))
  
  res_m <- rbind(res_m, missing_ids_m)
  
  res_m <- melt(data = res_m)
  colnames(res_m) <- c("ID", "CellType", "Freq")
  res_m <- as_tibble(res_m)
  dim(res_m)
  
  res <- md %>%
    group_by(ID) %>%
    summarise(Sample = first(Sample), Patient = first(Patient), Timepoint = first(Timepoint), .groups = "drop")
  
  res_m <- res_m %>%
    ungroup() %>%
    left_join(res %>%
                dplyr::select(Patient, Timepoint, ID),
              by = "ID")
  
  pts <- unique(res_m$Patient)
  states <- unique(res_m$CellType)
  
  t1_m <- lapply(states, function(s) res_m %>% filter(Timepoint == "T1", CellType == s) %>% pull(Freq))
  t1_m <- t(do.call(rbind, t1_m))
  colnames(t1_m) <- paste0("T1_ABS_TME_", states)
  rownames(t1_m) <- pts
  t1_m <- t1_m[pt_ord, ]
  
  t2_m <- lapply(states, function(s) t2 <- res_m %>% filter(Timepoint == "T2", CellType == s) %>% pull(Freq))
  t2_m <- t(do.call(rbind, t2_m))
  colnames(t2_m) <- paste0("T2_ABS_TME_", states)
  rownames(t2_m) <- pts
  t2_m <- t2_m[pt_ord, ]
  
  t2_vs_t1 <- lapply(1:nrow(t1_m), function(i) t2_m[i, ] - t1_m[i, ])
  t2_vs_t1 <- do.call(rbind, t2_vs_t1)
  colnames(t2_vs_t1) <- paste0("dTME_CellType_", states)
  rownames(t2_vs_t1) <- pt_ord
  
  res_m <- list(dT = t2_vs_t1,
                T1_ABS = t1_m,
                T2_ABS = t2_m)
  
  return(res_m)
}
