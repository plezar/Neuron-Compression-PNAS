library(RRHO2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)

# Neurons--------

res_Wbo2 <- read.csv("Figure_3/results/iN1_neuron.csv", row.names = 1)
res_I27  <- read.csv("Figure_3/results/iN2_neuron.csv", row.names = 1)
res_1019 <- read.csv("Figure_3/results/iN3_neuron.csv", row.names = 1)

res_Wbo2$pval_lfc <- -log10(res_Wbo2$pvalue)*res_Wbo2$log2FoldChange
res_I27$pval_lfc <- -log10(res_I27$pvalue)*res_I27$log2FoldChange
res_1019$pval_lfc <- -log10(res_1019$pvalue)*res_1019$log2FoldChange

# RRHO analysis------

run_RRHO2 <- function(df1, df2, label1 = "Sample1", label2 = "Sample2", rank_col) {

  df1$symbol <- rownames(df1)
  df2$symbol <- rownames(df2)
  
  df1 <- df1[, c("symbol", rank_col)]
  df2 <- df2[, c("symbol", rank_col)]
  
  shared_genes <- intersect(df1$symbol, df2$symbol)
  
  df1_sub <- df1[df1$symbol %in% shared_genes, ]
  df2_sub <- df2[df2$symbol %in% shared_genes, ]
  
  df1_sub <- df1_sub[match(shared_genes, df1_sub$symbol), ]
  df2_sub <- df2_sub[match(shared_genes, df2_sub$symbol), ]
  
  rr_up <- RRHO2_initialize(df1_sub,
                            df2_sub,
                            labels = c(label1, label2),
                            stepsize = nrow(df1_sub) / 100,
                            method = "hyper")
  
  return(rr_up)
}

spectral_continuous <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(256)


rr <- run_RRHO2(res_Wbo2, res_I27, "iN #1", "iN #2", "log2FoldChange")
pdf("Figure_3/figures/S3B_iN.pdf", width = 6, height = 5.5)
RRHO2_heatmap(rr, colorGradient = spectral_continuous)
dev.off()


# Scatterplot analysis------

compare_logFC_plot <- function(df1, df2,
                               label1 = "Dataset X",
                               label2 = "Dataset Y",
                               fc_thresh = 1,
                               return_stats = FALSE) {
  # Join by gene (rownames expected to be genes)
  df_joined <- df1 %>%
    tibble::rownames_to_column("gene") %>%
    inner_join(df2 %>% tibble::rownames_to_column("gene"), by = "gene")
  
  # Spearman correlation
  ct <- suppressWarnings(
    stats::cor.test(
      df_joined$log2FoldChange.x,
      df_joined$log2FoldChange.y,
      method = "spearman",
      exact = FALSE
    )
  )
  rho <- unname(ct$estimate)
  pval <- ct$p.value
  n <- nrow(df_joined)
  corr_label <- sprintf("rho = %.2f\np = %.2e", rho, pval)
  
  # Categorize genes
  df_plot <- df_joined %>%
    mutate(
      cat = case_when(
        log2FoldChange.x >  fc_thresh & log2FoldChange.y >  fc_thresh ~ "Concordant up",
        log2FoldChange.x < -fc_thresh & log2FoldChange.y < -fc_thresh ~ "Concordant down",
        (log2FoldChange.x >  fc_thresh & log2FoldChange.y < -fc_thresh) |
          (log2FoldChange.x < -fc_thresh & log2FoldChange.y >  fc_thresh) ~ "Discordant",
        TRUE ~ "Other"
      ),
      label = ifelse(cat %in% c("Concordant up", "Concordant down", "Discordant"),
                     gene, NA_character_)
    )
  
  pal <- c(
    "Concordant up"   = "#1b9e77",
    "Concordant down" = "#d95f02",
    "Discordant"      = "#7570b3",
    "Other"           = "grey80"
  )
  
  # Plot
  p <- ggplot(df_plot, aes(x = log2FoldChange.x, y = log2FoldChange.y, color = cat)) +
    geom_point(alpha = 0.65, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4, color = "grey55") +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", linewidth = 0.4, color = "grey70") +
    geom_hline(yintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", linewidth = 0.4, color = "grey70") +
    geom_text_repel(
      data = dplyr::filter(df_plot, !is.na(label)),
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0,
      seed = 42
    ) +
    annotate(
      "text",
      x = min(df_plot$log2FoldChange.x, na.rm = TRUE),
      y = max(df_plot$log2FoldChange.y, na.rm = TRUE),
      label = corr_label,
      hjust = 0, vjust = 1,
      size = 4,
      fontface = "italic"
    ) +
    scale_color_manual(values = pal, name = "Category") +
    labs(
      x = paste0("log2 fold change (", label1, ")"),
      y = paste0("log2 fold change (", label2, ")")
    ) +
    theme_cowplot(12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  if (isTRUE(return_stats)) {
    return(list(plot = p, spearman_rho = rho, p_value = pval, n = n))
  } else {
    return(p)
  }
}

p1 <- compare_logFC_plot(res_Wbo2, res_I27, label1 = "iN #1", label2 = "iN #2")
#p2 <- compare_logFC_plot(res_Wbo2, res_1019, label1 = "iN #1", label2 = "iN #3") + theme(legend.position = "none")
#p3 <- compare_logFC_plot(res_1019, res_I27, label1 = "iN #2", label2 = "iN #3") + theme(legend.position = "none")
#p <- p1 + p2 + p3 + plot_layout(widths = c(1, 1, 1))
ggsave(p1, path="Figure_3/figures", filename= "S3C_iN.pdf", width=4.5, height=3, dpi = 700, units = "in")


# Glia--------

res_Wbo2 <- read.csv("Figure_3/results/iN1_glia.csv", row.names = 1)
res_I27  <- read.csv("Figure_3/results/iN2_glia.csv", row.names = 1)
res_1019 <- read.csv("Figure_3/results/iN3_glia.csv", row.names = 1)

res_Wbo2$pval_lfc <- -log10(res_Wbo2$pvalue)*res_Wbo2$log2FoldChange
res_I27$pval_lfc <- -log10(res_I27$pvalue)*res_I27$log2FoldChange
res_1019$pval_lfc <- -log10(res_1019$pvalue)*res_1019$log2FoldChange

# RRHO analysis------

rr <- run_RRHO2(res_Wbo2, res_I27, "mG #1", "mG #2", "log2FoldChange")
pdf("Figure_3/figures/S3B_mG.pdf", width = 6, height = 5.5)
RRHO2_heatmap(rr, colorGradient = spectral_continuous)
dev.off()


# Scatterplot analysis------

p1 <- compare_logFC_plot(res_Wbo2, res_I27, label1 = "mG #1", label2 = "mG #2")
#p2 <- compare_logFC_plot(res_Wbo2, res_1019, label1 = "mG #1", label2 = "mG #3") + theme(legend.position = "none")
#p3 <- compare_logFC_plot(res_1019, res_I27, label1 = "mG #2", label2 = "mG #3") + theme(legend.position = "none")
#p <- p1 + p2 + p3 + plot_layout(widths = c(1, 1, 1))
ggsave(p1, path="Figure_3/figures", filename= "S3C_mG.pdf", width=4.5, height=3, dpi = 700, units = "in")

