#!/usr/bin/env Rscript

library(ggpubr)
library(cowplot)
library(tidyverse)
library(rstatix)

gene_data_file = "/Users/plezar/Library/CloudStorage/Box-Box/TIME Lab_Shared Folder/Personnel/Graduate Students/Maksym Zarodniuk/GBM ECM paper/qPCR/2025-01-11 IHA compression/targets.csv"
rq_results_file = "/Users/plezar/Library/CloudStorage/Box-Box/TIME Lab_Shared Folder/Personnel/Graduate Students/Maksym Zarodniuk/GBM ECM paper/qPCR/2025-01-11 IHA compression/rq_results.csv"
sample_data_file = "/Users/plezar/Library/CloudStorage/Box-Box/TIME Lab_Shared Folder/Personnel/Graduate Students/Maksym Zarodniuk/GBM ECM paper/qPCR/2025-01-11 IHA compression/samples.csv"

sample_data = read.csv(sample_data_file)
rq_results = read.csv(rq_results_file, comment.char = "#")
gene_data = read.csv(gene_data_file)

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
  delta_delta_eqcq[[i]] = delta_delta_eqcq[[i]] - mean(delta_delta_eqcq[[i]][delta_delta_eqcq$Refgroup==1])
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

stat.test_rq = rq %>%
  group_by(Target) %>%
  t_test(RQ ~ Biogroup) %>%
  adjust_pvalue(method="fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Biogroup", fun = "mean") %>%
  mutate(y.position = y.position)

stat.test_ddct$y.position = stat.test_rq$y.position

# TODO: add a parameter specifying whether want to display only near-significant p values
#stat.test_ddct = stat.test_ddct %>%
#  filter(p < 0.1)

rq_summarised_IHA <- rq_summarised %>% filter(Biogroup == "Comp") %>% mutate(RQ_geomean = log2(RQ_geomean)) %>% ungroup() %>% rename("IHA 24h" = RQ_geomean) %>% select(-c(Biogroup, RQ_SE))
rq_summarised_IHA <- as.data.frame(rq_summarised_IHA)
rownames(rq_summarised_IHA) <- rq_summarised_IHA$Target
rq_summarised_IHA$Target <- NULL
