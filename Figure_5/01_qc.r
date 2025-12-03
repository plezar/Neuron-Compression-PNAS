library(dplyr)
library(tidyr)
library(stringr)

source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_analysis_utils.R")
source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_CNA_utils.R")
#source("GBM-CARE-WT/R/GBM-CARE-WT_NMF.R")
source("/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/GBM-CARE-WT/R/GBM-CARE-WT_State_utils.R")

ANALYSIS_ROOT <- "/afs/crc.nd.edu/user/m/mzarodn2/Private/GSE274546/"
DATA_ROOT <- paste0(ANALYSIS_ROOT, "/data/")
TABLES_ROOT <- paste0(ANALYSIS_ROOT, "/tables/")
SDS_SAMPLE_PATH <- paste0(DATA_ROOT, "/sample_sds/")
NMSDS_SAMPLE_PATH <- paste0(DATA_ROOT, "/sample_nmsds/")

# Load the per-sample meta-data ----------------------------------------------------------------------------

sample_data <- as_tibble(read.csv(paste0(DATA_ROOT, "sample_data.csv"), stringsAsFactors = F))
sample_data$Gender[sample_data$Gender==""] <- NA
sample_data %>% fill(Gender, Age) -> sample_data

sample_data$Patient <- gsub(".*(P[0-9]+).*", "\\1", sample_data$Sample)
sample_data$Patient_factor <- factor(sample_data$Patient, unique(sample_data$Patient))
sample_data$Timepoint_factor <- factor(sample_data$Timepoint, unique(sample_data$Timepoint))
sample_data$Timepoint <- trimws(sample_data$Timepoint)
sample_data$Timepoint_num <- recode(sample_data$Timepoint, "Primary" = 0, "1st Recurrent" = 1, "2nd Recurrent" = 2)
sample_data$ID <- paste0(sample_data$Patient, "_", sample_data$Timepoint_num)
#sample_data$ID_extended <- paste0(sample_data$ID, "_", sample_data$Sample)
sample_data$ID_factor <- factor(sample_data$ID, unique(sample_data$ID))
sample_data$Sample <- factor(sample_data$Sample, unique(sample_data$Sample))

table(sample_data$ID_factor)
saveRDS(object = sample_data, file = paste0(DATA_ROOT, "sample_data.RDS"))


# Create per-sample list -------------------------------------------------------------------------------------

dirs_list <- list.files(DATA_ROOT, full.names = T)
dirs_list <- dirs_list[grepl("umi_counts.RDS", dirs_list)]
#dirs_list <- dirs_list[1:2]

umi_data_all <- lapply(dirs_list, function(s) {
  message(basename(s))
  umi_data <- readRDS(s)
  #colnames(umi_data) <- paste0(s, "-", colnames(umi_data))
  umi_data
})
gc()

names(umi_data_all) <- str_split(unlist(lapply(dirs_list, basename)), "_", Inf, T)[,2]
saveRDS(object = umi_data_all, file = paste0(DATA_ROOT, "umi_data_list.RDS"))


# Load cell-level meta-data ----------------------------------------------------------------------------

cellids <- unname(unlist(lapply(umi_data_all, colnames), recursive = F))
message(paste0("Number of cells: ", length(cellids)))

meta_data <- tibble(CellID = cellids,
                    Sample = rep(names(umi_data_all), times = sapply(umi_data_all, ncol)))

message("Per-sample cell numbers:")
print(table(meta_data$Sample))

meta_data <- meta_data %>% left_join(sample_data, by = "Sample")
meta_data <- as.data.frame(meta_data)
rownames(meta_data) <- meta_data$CellID
dim(meta_data)

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data.RDS"))

rm(umi_data_all)
gc()

# Basic QC ------------------------------------------------------------------------
#registered()

## Compute cell QC metrics (complexity, MT expression) ------------------------------------------------


qc_df <- data.frame()
for (s in dirs_list) {
  message(basename(s))
  umi_data <- readRDS(s)
  c <- umi_data != 0
  c <- colSums(c)

  mt_genes <- rownames(umi_data)[grep("^MT-", rownames(umi_data))]
  mt_sum <- colSums(umi_data[mt_genes, ])
  all_sum <- colSums(umi_data)
  mt_pct <- setNames((mt_sum / all_sum) * 100, colnames(umi_data))

  qc_df <- rbind(qc_df, cbind(c, mt_pct))
  gc()
}


meta_data$Complexity <- qc_df[meta_data$CellID, ]$c
meta_data$Percent.MT <- qc_df[meta_data$CellID, ]$mt_pct

qc_stats <- meta_data %>%
  group_by(Patient_factor, Sample) %>%
  summarise(Ncells_PreQC = n(),
            MeanComplexity = mean(Complexity),
            MedianComplexity = median(Complexity),
            Q005Complexity = quantile(Complexity, .005),
            Q025Complexity = quantile(Complexity, .025),
            Q25Complexity = quantile(Complexity, .25),
            Q75Complexity = quantile(Complexity, .75),
            Q975Complexity = quantile(Complexity, .975))

meta_data <- meta_data %>% filter(Complexity >= 200 & Complexity <= 10000, Percent.MT <= 3)

saveRDS(object = meta_data, file = paste0(DATA_ROOT, "meta_data_qc.RDS"))

qc_stats <- left_join(x = meta_data %>% group_by(Sample) %>% summarise(Ncells_PostQC = n()),
                      y = qc_stats,
                      by = "Sample")
qc_stats$PercentCellsDropped <- round((1 - qc_stats$Ncells_PostQC / qc_stats$Ncells_PreQC) * 100, 1)

qc_stats <- qc_stats %>% select(Sample,
                                Ncells_PreQC, Ncells_PostQC, PercentCellsDropped,
                                MeanComplexity, MedianComplexity,
                                Q025Complexity, Q25Complexity, Q75Complexity, Q975Complexity)

qc_stats %>%
  arrange(Sample) %>%
  write.csv(file = paste0(TABLES_ROOT, "qc_stats_tbl.tsv"))

message("01_QC.r is complete")
