library(tidyverse)
library(rstatix)
library(ggpubr)
library(cowplot)
library(patchwork)

scale_factor <- 2

# ==== Prepare data =====
data_path <- "Figure_1/data"
csv_files <- list.files(data_path, full.names = F)

all_d <- data.frame()
for (i in 1:length(csv_files)) {
  cfile <- paste0(data_path, "/", csv_files[i])
  
  d <- read.csv(cfile)
  d$id <- 1:nrow(d)
  
  d <- pivot_longer(d, !id)
  d$stat <- str_split(csv_files[i], "\\.", Inf, T)[1]
  d$donor <- str_split(d$name, "_", Inf, T)[,1]
  d$donor <- factor(d$donor)
  levels(d$donor) <- c("iN #2", "iN #1", "iN #3")
  d$condition <- str_split(d$name, "_", Inf, T)[,2]
  d$name <- NULL
  
  all_d <- rbind(all_d, d)
}

all_d <- all_d %>% pivot_wider(names_from = "stat", values_from = "value")
colnames(all_d) <- c("id", "donor", "condition", "GFP", "DAPI", "PI.norm", "GFP_AND_PI.norm", "PI_MINUS_GFP_freq", "PI")
all_d$condition <- factor(all_d$condition, levels = c("Ctrl", "Comp"))
all_d$group <- paste(all_d$donor, all_d$condition, sep=" ")
all_d$group <- factor(all_d$group, levels = c("iN #1 Ctrl", "iN #1 Comp", "iN #2 Ctrl", "iN #2 Comp", "iN #3 Ctrl", "iN #3 Comp"))

# ===== Figure 1 =====
##
##
## ==== GFP count =====

stat.test <- all_d %>%
  group_by(donor) %>%
  wilcox_test(GFP ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", fun = "max", step.increase = 0.2)

gfp_count <- ggbarplot(all_d,
                       x = "group",
                       y = "GFP",
                       add = c("mean_se"),
                       fill = "group",
                       linewidth = 0.5,
                       error.plot = "upper_errorbar") +
  geom_point(data = all_d, aes(x=group, y = GFP, fill=group), shape=21, alpha = 0.4, size=.3*scale_factor, stroke = 0.5, position=position_jitter(height=0, width=0)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005, size=4*scale_factor, hide.ns=T) +
  theme_cowplot(10*scale_factor) +
  scale_y_continuous(limits = c(0, 270)) +
  labs(y="GFP+ count", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values = c("#79CFE2", "#5862AC", "#F2775F", "#DF2C26", "#AAD05A", "#098C45"))


## ==== GFP and PI frac. =====

stat.test <- all_d %>%
  group_by(donor) %>%
  wilcox_test(GFP_AND_PI.norm ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", fun = "max", step.increase = 0.2)

stat.test$y.position <- c(1.05, 1.05, 0.52)

pi_gfp_freq <- ggbarplot(all_d,
                         x = "group",
                         y = "GFP_AND_PI.norm",
                         add = c("mean_se"),
                         fill = "group",
                         linewidth = 0.5,
                         error.plot = "upper_errorbar") +
  geom_point(data = all_d, aes(x=group, y = GFP_AND_PI.norm, fill=group), shape=21, alpha = 0.4, stroke = 0.5, size=.3*scale_factor, position=position_jitter(height=0, width=0)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005,  size=4*scale_factor, hide.ns=T) +
  theme_cowplot(10*scale_factor) +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(y="PI+ frac. of GFP+", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values = c("#79CFE2", "#5862AC", "#F2775F", "#DF2C26", "#AAD05A", "#098C45"))


## ==== PI counts =====

stat.test <- all_d %>%
  group_by(donor) %>%
  wilcox_test(PI ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", fun = "max", step.increase = 0.2)

pi_count <- ggbarplot(all_d,
                      x = "group",
                      y = "PI",
                      add = c("mean_se"),
                      fill = "group",
                      linewidth = 0.5,
                      error.plot = "upper_errorbar") +
  geom_point(data = all_d, aes(x=group, y = PI, fill=group), shape=21, alpha = 0.4, stroke = 0.5, size=.3*scale_factor, position=position_jitter(height=0, width=0)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005,  size=4*scale_factor, hide.ns=T) +
  theme_cowplot(10*scale_factor) +
  scale_y_continuous(limits = c(0, 1000)) +
  labs(y="PI+ count", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values = c("#79CFE2", "#5862AC", "#F2775F", "#DF2C26", "#AAD05A", "#098C45"))


## ==== PI frac. =====

stat.test <- all_d %>%
  group_by(donor) %>%
  wilcox_test(PI.norm ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", fun = "max", step.increase = 0.2)

pi_freq <- ggbarplot(all_d,
                     x = "group",
                     y = "PI.norm",
                     add = c("mean_se"),
                     fill = "group",
                     linewidth = 0.5,
                     error.plot = "upper_errorbar") +
  geom_point(data = all_d, aes(x=group, y = PI.norm, fill=group), shape=21, alpha = 0.4, stroke = 0.5, size=.3*scale_factor, position=position_jitter(height=0, width=0)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005,  size=4*scale_factor, hide.ns=T) +
  theme_cowplot(10*scale_factor) +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(y="PI+ frac. of Hoechst+", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values = c("#79CFE2", "#5862AC", "#F2775F", "#DF2C26", "#AAD05A", "#098C45"))

p_D <- plot_spacer() + pi_count + plot_spacer() + pi_freq + plot_spacer() + gfp_count + plot_spacer() + pi_gfp_freq + plot_layout(nrow = 1, widths = c(0.05, 1, 0.1, 1, 0.1, 1, 0.1, 1))
ggsave(p_D, path="Figure_1/results/", filename= "Fig1_D.pdf", width=18*scale_factor, height=5*scale_factor, dpi = 700, units = "cm")

# ==== Figure S1 =====
##
##
## ===== PI+GFP- frac. =====

stat.test <- all_d %>%
  group_by(donor) %>%
  wilcox_test(PI_MINUS_GFP_freq ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", fun = "max", step.increase = 0.2)

pi_minus_gfp_freq <- ggbarplot(all_d,
                               x = "group",
                               y = "PI_MINUS_GFP_freq",
                               add = c("mean_se"),
                               fill = "group",
                               linewidth = 0.5,
                               error.plot = "upper_errorbar") +
  geom_point(data = all_d, aes(x=group, y = PI_MINUS_GFP_freq, fill=group), shape=21, alpha = 0.4, stroke = 0.5, size=.3*scale_factor, position=position_jitter(height=0, width=0)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005,  size=4*scale_factor, hide.ns=T) +
  theme_cowplot(10*scale_factor) +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(y="PI+GFP- frac. of Hoechst+", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y = element_text(size = 16), legend.position = "none") +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values = c("#79CFE2", "#5862AC", "#F2775F", "#DF2C26", "#AAD05A", "#098C45"))

## ===== Total counts =====

stat.test <- all_d %>%
  group_by(donor) %>%
  wilcox_test(DAPI ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "group", fun = "max", step.increase = 0.2)

hoechst_count <- ggbarplot(all_d,
                           x = "group",
                           y = "DAPI",
                           add = c("mean_se"),
                           fill = "group",
                           linewidth = 0.5,
                           error.plot = "upper_errorbar") +
  geom_point(data = all_d, aes(x=group, y = DAPI, fill=group), shape=21, alpha = 0.4, stroke = 0.5, size=.3*scale_factor, position=position_jitter(height=0, width=0)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005,  size=4*scale_factor, hide.ns=T) +
  theme_cowplot(10*scale_factor) +
  scale_y_continuous(limits = c(0, 1200)) +
  labs(y="Hoechst+ count", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  scale_alpha_discrete(range = c(0.4,1)) +
  scale_fill_manual(values = c("#79CFE2", "#5862AC", "#F2775F", "#DF2C26", "#AAD05A", "#098C45"))


p_A <- pi_minus_gfp_freq + plot_spacer() + hoechst_count + plot_layout(widths = c(1, 0.5, 1))
ggsave(p_A, path="Figure_1/results/", filename= "Fig1S.pdf", width=11*scale_factor, height=5*scale_factor, dpi = 700, units = "cm")