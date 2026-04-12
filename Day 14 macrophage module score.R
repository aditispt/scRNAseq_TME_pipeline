##Quantitiative analysis
# ================================================== #
# MACROPHAGE MODULE SCORING + COMPARISON 
# Uses existing objects: integrated, sna, snahigh
# ================================================== #

library(Seurat)
library(dplyr)
library(ggplot2)

# -------------------------------------------------- #
# 1. DEFINE GENE MODULES (based on YOUR DEG results)
# -------------------------------------------------- #

# ================================================== #
# EXPANDED MACROPHAGE GENE MODULES (ROBUST SET)
# ================================================== #

# IFN / ANTIVIRAL
IFN_genes <- c(
  "Ifitm1","Ifitm2","Ifitm3",
  "Rsad2","Ifit1","Ifit2","Ifit3",
  "Isg15","Irf7","Stat1",
  "Oas1a","Oas2","Mx1",
  "Usp18","Herc6"
)

# ECM / REMODELING
ECM_genes <- c(
  "Col1a1","Col1a2","Col3a1",
  "Mmp9","Mmp12","Mmp14",
  "Timp1","Pcolce",
  "Fn1","Spp1","Mgp",
  "Thbs1","Lox","Loxl2"
)

# M2 / IMMUNOSUPPRESSIVE
M2_genes <- c(
  "Retnla","Arg1","Mrc1",
  "Il10","C1qa","C1qb","C1qc",
  "Cd163","Trem2","Lgals1",
  "Msr1","Chil3","Ear2"
)

# METABOLIC (glycolysis / mitochondrial adaptation)
Metabolic_genes <- c(
  "Tars","Slc25a44","Gstt2",
  "Gapdh","Pgk1","Ldha",
  "Acly","Idh1","Pdk1",
  "Cs","Sdha","Ndufa1"
)

# HYPOXIA / STRESS (tumor-relevant)
Hypoxia_genes <- c(
  "Hif1a","Vegfa","Hilpda",
  "Slc2a1","Pdk1","Bnip3",
  "Egln3","Ndrg1"
)

# INFLAMMATORY (non-IFN activation)
Inflammatory_genes <- c(
  "Il1b","Tnf","Cxcl2","Cxcl3",
  "Nfkbia","Nfkbiz","Ptgs2",
  "Ccl3","Ccl4"
)

# -------------------------------------------------- #
# 2. ADD MODULE SCORES (DAY 14 object: integrated)
# -------------------------------------------------- #

integrated <- AddModuleScore(
  object = integrated,
  features = list(
    IFN_genes,
    ECM_genes,
    M2_genes,
    Metabolic_genes,
    Hypoxia_genes,
    Inflammatory_genes
  ),
  name = c("IFN", "ECM", "M2", "Metabolic", "Hypoxia", "Inflammatory")
)

# Clean naming (important)
integrated$IFN_score          <- integrated$IFN1
integrated$ECM_score          <- integrated$ECM2
integrated$M2_score           <- integrated$M23
integrated$Metabolic_score    <- integrated$Metabolic4
integrated$Hypoxia_score      <- integrated$Hypoxia5
integrated$Inflammatory_score <- integrated$Inflammatory6

# -------------------------------------------------- #
# 3. SUBSET MACROPHAGES (same naming as your code)
# -------------------------------------------------- #

macro_d14 <- subset(
  integrated,
  idents = c("Macrophage_1", "Macrophage_2", "Macrophage_3")
)

# -------------------------------------------------- #
# 4. SNAhigh vs SNAlow COMPARISON (Day 14)
# -------------------------------------------------- #

# Violin plots
VlnPlot(
  macro_d14,
  features = c(
    "IFN_score",
    "ECM_score",
    "M2_score",
    "Metabolic_score",
    "Hypoxia_score",
    "Inflammatory_score"
  ),
  group.by = "sample",
  pt.size = 0
)

# -------------------------------------------------- #
# 5. STATISTICS (Wilcoxon test)
# -------------------------------------------------- #

meta_d14 <- macro_d14@meta.data

wilcox_IFN  <- wilcox.test(IFN_score ~ sample, data = meta_d14)
wilcox_ECM  <- wilcox.test(ECM_score ~ sample, data = meta_d14)
wilcox_M2   <- wilcox.test(M2_score ~ sample, data = meta_d14)
wilcox_MET  <- wilcox.test(Metabolic_score ~ sample, data = meta_d14)
wilcox_HYP  <- wilcox.test(Hypoxia_score ~ sample, data = meta_d14)
wilcox_INF  <- wilcox.test(Inflammatory_score ~ sample, data = meta_d14)

# Print results
wilcox_IFN
wilcox_ECM
wilcox_M2
wilcox_MET
wilcox_HYP
wilcox_INF

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# -------------------------------------------------- #
# Prepare data
# -------------------------------------------------- #

plot_data <- macro_d14@meta.data %>%
  select(
    sample,
    IFN_score,
    ECM_score,
    M2_score,
    Metabolic_score,
    Hypoxia_score,
    Inflammatory_score
  ) %>%
  pivot_longer(
    cols = -sample,
    names_to = "module",
    values_to = "score"
  )

# -------------------------------------------------- #
# Violin plot with significance
# -------------------------------------------------- #
scale_fill_manual(values = c("855" = "blue", "857" = "red"))

ggplot(plot_data, aes(x = sample, y = score, fill = sample)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~module, scales = "free_y") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    label.y.npc = 0.9,
    label.x = 1.5,
    hide.ns = TRUE,
    size = 6
  ) +
  scale_fill_manual(values = c("855" = "blue", "857" = "red")) +
  theme_classic() +
  theme(strip.text = element_text(face = "bold")) +
  ggtitle("Macrophage Functional Programs (Day 14)")
