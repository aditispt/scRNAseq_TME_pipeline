# ============================================================ #
#   CELL-TYPE-SPECIFIC DEG + MODULE SCORING + TEMPORAL ANALYSIS
#   Cell types: T cells (+ NK), MDSC, Epithelial, Fibroblasts
#   Timepoints: Day 14 (integrated), Day 25 (sna)
#   SNAhigh vs SNAlow comparison at D14 & D25
#   NOTE: Day 40 excluded — no SNAlow control available
# ============================================================ #

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(broom)

# ------------------------------------------------------------ #
# STEP 0: ADD CONDITION COLUMN TO DAY 14 OBJECT
# Maps sample IDs to SNAhigh / SNAlow
# ------------------------------------------------------------ #

integrated$condition <- ifelse(integrated$sample == "855", "SNAhigh", "SNAlow")


# ============================================================ #
#   SECTION 1: T CELLS (including NK_Cells)
# ============================================================ #

# ------------------------------------------------------------ #
# 1A. DEFINE T CELL GENE MODULES
# ------------------------------------------------------------ #

# Cytotoxicity / Effector function
Tcell_Cytotoxic_genes <- c(
  "Gzmb", "Gzma", "Gzmk", "Prf1",
  "Ifng", "Tnf", "Fasl",
  "Nkg7", "Klrg1", "Cx3cr1"
)

# Exhaustion
Tcell_Exhaustion_genes <- c(
  "Pdcd1",    # PD-1
  "Havcr2",   # TIM-3
  "Lag3",
  "Tigit",
  "Ctla4",
  "Tox",
  "Nr4a1",
  "Entpd1",   # CD39
  "Cd244"
)

# Memory / Stemness
Tcell_Memory_genes <- c(
  "Tcf7",     # TCF1
  "Sell",     # CD62L
  "Ccr7",
  "Il7r",
  "Lef1",
  "Bcl6",
  "Slamf6",
  "Cxcr5"
)

# Regulatory T cell (Treg)
Tcell_Treg_genes <- c(
  "Foxp3",
  "Il2ra",    # CD25
  "Ikzf2",    # Helios
  "Ctla4",
  "Il10",
  "Tgfb1",
  "Tnfrsf18"  # GITR
)

# NK / Innate-like
NK_Innate_genes <- c(
  "Ncr1",     # NKp46
  "Klrb1c",   # NK1.1
  "Klrd1",    # CD94
  "Eomes",
  "Tbx21",    # T-bet
  "Xcl1",
  "Klrc2",
  "Ly6c2"
)

# Proliferation
Tcell_Proliferation_genes <- c(
  "Mki67",
  "Top2a",
  "Pcna",
  "Cdk1",
  "Birc5",
  "Stmn1"
)

# ------------------------------------------------------------ #
# 1B. SUBSET T CELLS FROM EACH OBJECT
# ------------------------------------------------------------ #

# Day 14 — T_cells_1, T_cells_2, NK_Cells
tcell_d14 <- subset(
  integrated,
  idents = c("T_cells_1", "T_cells_2", "NK_Cells")
)

# Day 25 — T_cells, NK_Cells (adjust if your object has different names)
tcell_d25 <- subset(
  sna,
  idents = c("T_cells")
)

# ------------------------------------------------------------ #
# 1C. DEG: SNAhigh vs SNAlow (Day 14)
# ------------------------------------------------------------ #
DefaultAssay(tcell_d14) <- "RNA"

tcell_d14 <- NormalizeData(tcell_d14)
tcell_d14 <- FindVariableFeatures(tcell_d14)
tcell_d14 <- ScaleData(tcell_d14)

Idents(tcell_d14) <- "condition"

deg_tcell_d14 <- FindMarkers(
  tcell_d14,
  ident.1 = "SNAhigh",
  ident.2 = "SNAlow",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

deg_tcell_d14$gene <- rownames(deg_tcell_d14)
deg_tcell_d14$timepoint <- "Day14"

write.csv(deg_tcell_d14, "DEG_Tcells_Day14_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

# ------------------------------------------------------------ #
# 1D. DEG: SNAhigh vs SNAlow (Day 25)
# ------------------------------------------------------------ #

Idents(tcell_d25) <- "SNA_level_D25"

deg_tcell_d25 <- FindMarkers(
  tcell_d25,
  ident.1 = "SNAhigh",
  ident.2 = "SNAlow",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

deg_tcell_d25$gene <- rownames(deg_tcell_d25)
deg_tcell_d25$timepoint <- "Day25"

write.csv(deg_tcell_d25, "DEG_Tcells_Day25_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

# ------------------------------------------------------------ #
# 1E. VOLCANO PLOTS — T Cells (Day 14 & Day 25)
# ------------------------------------------------------------ #

plot_volcano <- function(deg_df, title = "Volcano Plot", 
                         fc_cut = 0.5, p_cut = 0.05) {
  deg_df <- deg_df %>%
    mutate(
      sig = case_when(
        avg_log2FC >  fc_cut & p_val_adj < p_cut ~ "Up in SNAhigh",
        avg_log2FC < -fc_cut & p_val_adj < p_cut ~ "Up in SNAlow",
        TRUE ~ "NS"
      )
    )
  
  top_labels <- deg_df %>%
    filter(sig != "NS") %>%
    arrange(p_val_adj) %>%
    slice_head(n = 20)
  
  ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text(
      data = top_labels,
      aes(label = gene),
      size = 3, hjust = -0.1, vjust = 0.5, check_overlap = TRUE
    ) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c(
      "Up in SNAhigh" = "#2166AC",
      "Up in SNAlow"  = "#D6604D",
      "NS"            = "grey70"
    )) +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    ) +
    labs(title = title, x = "avg log2FC", y = "-log10(adj. p-value)")
}

volcano_tcell_d14 <- plot_volcano(
  deg_tcell_d14,
  title = "T Cells — Day 14: SNAhigh vs SNAlow"
)

volcano_tcell_d25 <- plot_volcano(
  deg_tcell_d25,
  title = "T Cells — Day 25: SNAhigh vs SNAlow"
)

ggsave("Volcano_Tcell_D14.png", volcano_tcell_d14, width = 8, height = 6, dpi = 300)
ggsave("Volcano_Tcell_D25.png", volcano_tcell_d25, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------ #
# 1F. MODULE SCORING — T Cells (Day 14 & Day 25)
# ------------------------------------------------------------ #

tcell_modules <- list(
  Cytotoxic      = Tcell_Cytotoxic_genes,
  Exhaustion     = Tcell_Exhaustion_genes,
  Memory         = Tcell_Memory_genes,
  Treg           = Tcell_Treg_genes,
  NK_Innate      = NK_Innate_genes,
  Proliferation  = Tcell_Proliferation_genes
)

# Day 14
tcell_d14 <- AddModuleScore(
  object   = tcell_d14,
  features = tcell_modules,
  name     = names(tcell_modules)
)

tcell_d14$Cytotoxic_score     <- tcell_d14$Cytotoxic1
tcell_d14$Exhaustion_score    <- tcell_d14$Exhaustion2
tcell_d14$Memory_score        <- tcell_d14$Memory3
tcell_d14$Treg_score          <- tcell_d14$Treg4
tcell_d14$NK_Innate_score     <- tcell_d14$NK_Innate5
tcell_d14$Proliferation_score <- tcell_d14$Proliferation6

# Day 25
tcell_d25 <- AddModuleScore(
  object   = tcell_d25,
  features = tcell_modules,
  name     = names(tcell_modules)
)

tcell_d25$Cytotoxic_score     <- tcell_d25$Cytotoxic1
tcell_d25$Exhaustion_score    <- tcell_d25$Exhaustion2
tcell_d25$Memory_score        <- tcell_d25$Memory3
tcell_d25$Treg_score          <- tcell_d25$Treg4
tcell_d25$NK_Innate_score     <- tcell_d25$NK_Innate5
tcell_d25$Proliferation_score <- tcell_d25$Proliferation6

# ------------------------------------------------------------ #
# 1G. VIOLIN PLOTS WITH STATS — T Cells
# ------------------------------------------------------------ #

score_cols_tcell <- c(
  "Cytotoxic_score", "Exhaustion_score", "Memory_score",
  "Treg_score", "NK_Innate_score", "Proliferation_score"
)

plot_module_violin <- function(meta_df, group_col, score_cols, 
                               group_levels, group_colors, title = "") {
  
  plot_data <- meta_df %>%
    select(all_of(c(group_col, score_cols))) %>%
    pivot_longer(cols = -all_of(group_col),
                 names_to = "module", values_to = "score")
  
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = group_levels)
  
  ggplot(plot_data, aes(x = .data[[group_col]], y = score,
                        fill = .data[[group_col]])) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    facet_wrap(~module, scales = "free_y") +
    stat_compare_means(
      method     = "wilcox.test",
      label      = "p.signif",
      label.y.npc = 0.9,
      label.x    = 1.5,
      hide.ns    = TRUE,
      size       = 6
    ) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    ) +
    labs(title = title, x = NULL, y = "Module Score")
}

condition_colors <- c("SNAhigh" = "#2166AC", "SNAlow" = "#D6604D")

vln_tcell_d14 <- plot_module_violin(
  meta_df       = tcell_d14@meta.data,
  group_col     = "condition",
  score_cols    = score_cols_tcell,
  group_levels  = c("SNAhigh", "SNAlow"),
  group_colors  = condition_colors,
  title         = "T Cell Programs — Day 14"
)

vln_tcell_d25 <- plot_module_violin(
  meta_df       = tcell_d25@meta.data,
  group_col     = "SNA_level_D25",
  score_cols    = score_cols_tcell,
  group_levels  = c("SNAhigh", "SNAlow"),
  group_colors  = condition_colors,
  title         = "T Cell Programs — Day 25"
)

ggsave("ModuleScore_Tcell_D14.png", vln_tcell_d14, width = 12, height = 8, dpi = 300)
ggsave("ModuleScore_Tcell_D25.png", vln_tcell_d25, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------ #
# 1H. TEMPORAL TRAJECTORY — T Cells (D14 → D25)
# ------------------------------------------------------------ #

meta_tcell_d14 <- tcell_d14@meta.data %>%
  select(condition, all_of(score_cols_tcell)) %>%
  dplyr::rename(SNA_level = condition) %>%
  mutate(timepoint = "Day14")

meta_tcell_d25 <- tcell_d25@meta.data %>%
  select(SNA_level_D25, all_of(score_cols_tcell)) %>%
  dplyr::rename(SNA_level = SNA_level_D25) %>%
  mutate(timepoint = "Day25")

combined_tcell <- bind_rows(meta_tcell_d14, meta_tcell_d25) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Day14", "Day25")),
    SNA_level = factor(SNA_level, levels = c("SNAhigh", "SNAlow"))
  ) %>%
  pivot_longer(
    cols      = all_of(score_cols_tcell),
    names_to  = "module",
    values_to = "score"
  )

# Summary for line plot
summary_tcell <- combined_tcell %>%
  dplyr::group_by(module, SNA_level, timepoint) %>%
  dplyr::summarise(
    mean_score = mean(score, na.rm = TRUE),
    se         = sd(score, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )

traj_tcell <- ggplot(summary_tcell,
                     aes(x = timepoint, y = mean_score,
                         group = SNA_level, color = SNA_level)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                width = 0.1) +
  facet_wrap(~module, scales = "free_y") +
  scale_color_manual(values = c("SNAhigh" = "#2166AC", "SNAlow" = "#D6604D")) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "T Cell Program Trajectories: Day 14 → Day 25",
       x = "Timepoint", y = "Mean Module Score", color = "SNA Level")

ggsave("Trajectory_Tcell_D14_D25.png", traj_tcell, width = 14, height = 8, dpi = 300)

# Interaction model — does SNA level modulate temporal change?
interaction_tcell <- combined_tcell %>%
  group_by(module) %>%
  do(tidy(lm(score ~ timepoint * SNA_level, data = .))) %>%
  filter(term == "timepointDay25:SNA_levelSNAlow")

print(interaction_tcell %>% select(module, estimate, p.value))


# ============================================================ #
#   SECTION 2: MDSCs
# ============================================================ #

# ------------------------------------------------------------ #
# 2A. DEFINE MDSC GENE MODULES
# ------------------------------------------------------------ #

# Immunosuppression
MDSC_Suppression_genes <- c(
  "Arg1",
  "Nos2",
  "Il10",
  "Tgfb1",
  "S100a8", "S100a9",
  "Ptgs2",
  "Cybb",
  "Ido1"
)

# Myeloid / Granulocytic identity
MDSC_Myeloid_genes <- c(
  "Cd33",
  "Itgam",   # CD11b
  "Ly6g",
  "Ly6c1",
  "Cebpb",
  "Irf8",
  "Csf1r",
  "Fcgr3"
)

# ROS / Oxidative stress
MDSC_Oxidative_genes <- c(
  "Cyba",
  "Ncf1",
  "Ncf2",
  "Sod2",
  "Gpx1",
  "Txnrd1",
  "Hmox1",
  "Nfe2l2"
)

# Recruitment / Migration
MDSC_Recruitment_genes <- c(
  "Cxcr2",
  "Cxcr4",
  "Ccr2",
  "Cxcl2",
  "S100a8",
  "Mmp8",
  "Cd44"
)

# Differentiation block (immature markers)
MDSC_Immature_genes <- c(
  "Cd34",
  "Kit",
  "Flt3",
  "Gfi1",
  "Cebpa",
  "Mpo",
  "Elane"
)

mdsc_modules <- list(
  Suppression  = MDSC_Suppression_genes,
  Myeloid      = MDSC_Myeloid_genes,
  Oxidative    = MDSC_Oxidative_genes,
  Recruitment  = MDSC_Recruitment_genes,
  Immature     = MDSC_Immature_genes
)

score_cols_mdsc <- c(
  "Suppression_score", "Myeloid_score", "Oxidative_score",
  "Recruitment_score", "Immature_score"
)

# ------------------------------------------------------------ #
# 2B. SUBSET MDSCs
# ------------------------------------------------------------ #

mdsc_d14 <- subset(integrated, idents = "MDSC")
mdsc_d25 <- subset(sna,        idents = "MDSC")

# ------------------------------------------------------------ #
# 2C. DEG: SNAhigh vs SNAlow
# ------------------------------------------------------------ #

# Switch to RNA assay
DefaultAssay(mdsc_d14) <- "RNA"

# Normalize (safe even if already normalized)
mdsc_d14 <- NormalizeData(mdsc_d14)
mdsc_d14 <- FindVariableFeatures(mdsc_d14)
mdsc_d14 <- ScaleData(mdsc_d14)

# Set identities
Idents(mdsc_d14) <- "condition"

# Run DEG
deg_mdsc_d14 <- FindMarkers(
  mdsc_d14,
  ident.1 = "SNAhigh",
  ident.2 = "SNAlow",
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_mdsc_d14$gene <- rownames(deg_mdsc_d14)
write.csv(deg_mdsc_d14, "DEG_MDSC_Day14_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

Idents(mdsc_d25) <- "SNA_level_D25"
deg_mdsc_d25 <- FindMarkers(mdsc_d25, ident.1 = "SNAhigh", ident.2 = "SNAlow",
                            min.pct = 0.25, logfc.threshold = 0.25)
deg_mdsc_d25$gene <- rownames(deg_mdsc_d25)
write.csv(deg_mdsc_d25, "DEG_MDSC_Day25_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

# Volcano plots
volcano_mdsc_d14 <- plot_volcano(deg_mdsc_d14, "MDSC — Day 14: SNAhigh vs SNAlow")
volcano_mdsc_d25 <- plot_volcano(deg_mdsc_d25, "MDSC — Day 25: SNAhigh vs SNAlow")
ggsave("Volcano_MDSC_D14.png", volcano_mdsc_d14, width = 8, height = 6, dpi = 300)
ggsave("Volcano_MDSC_D25.png", volcano_mdsc_d25, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------ #
# 2D. MODULE SCORING — MDSCs
# ------------------------------------------------------------ #

mdsc_d14 <- AddModuleScore(mdsc_d14, features = mdsc_modules, name = names(mdsc_modules))
mdsc_d14$Suppression_score  <- mdsc_d14$Suppression1
mdsc_d14$Myeloid_score      <- mdsc_d14$Myeloid2
mdsc_d14$Oxidative_score    <- mdsc_d14$Oxidative3
mdsc_d14$Recruitment_score  <- mdsc_d14$Recruitment4
mdsc_d14$Immature_score     <- mdsc_d14$Immature5

mdsc_d25 <- AddModuleScore(mdsc_d25, features = mdsc_modules, name = names(mdsc_modules))
mdsc_d25$Suppression_score  <- mdsc_d25$Suppression1
mdsc_d25$Myeloid_score      <- mdsc_d25$Myeloid2
mdsc_d25$Oxidative_score    <- mdsc_d25$Oxidative3
mdsc_d25$Recruitment_score  <- mdsc_d25$Recruitment4
mdsc_d25$Immature_score     <- mdsc_d25$Immature5

# Violin plots
vln_mdsc_d14 <- plot_module_violin(mdsc_d14@meta.data, "condition",
                                   score_cols_mdsc, c("SNAhigh","SNAlow"),
                                   condition_colors, "MDSC Programs — Day 14")
vln_mdsc_d25 <- plot_module_violin(mdsc_d25@meta.data, "SNA_level_D25",
                                   score_cols_mdsc, c("SNAhigh","SNAlow"),
                                   condition_colors, "MDSC Programs — Day 25")
ggsave("ModuleScore_MDSC_D14.png", vln_mdsc_d14, width = 12, height = 8, dpi = 300)
ggsave("ModuleScore_MDSC_D25.png", vln_mdsc_d25, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------ #
# 2E. TEMPORAL TRAJECTORY — MDSCs
# ------------------------------------------------------------ #

meta_mdsc_d14 <- mdsc_d14@meta.data %>%
  select(condition, all_of(score_cols_mdsc)) %>%
  dplyr::rename(SNA_level = condition) %>%
  mutate(timepoint = "Day14")

meta_mdsc_d25 <- mdsc_d25@meta.data %>%
  select(SNA_level_D25, all_of(score_cols_mdsc)) %>%
  dplyr::rename(SNA_level = SNA_level_D25) %>%
  mutate(timepoint = "Day25")

combined_mdsc <- bind_rows(meta_mdsc_d14, meta_mdsc_d25) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Day14", "Day25")),
    SNA_level = factor(SNA_level, levels = c("SNAhigh", "SNAlow"))
  ) %>%
  pivot_longer(cols = all_of(score_cols_mdsc),
               names_to = "module", values_to = "score")

summary_mdsc <- combined_mdsc %>%
  dplyr::group_by(module, SNA_level, timepoint) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE),
            se = sd(score, na.rm = TRUE) / sqrt(n()), .groups = "drop")

traj_mdsc <- ggplot(summary_mdsc,
                    aes(x = timepoint, y = mean_score,
                        group = SNA_level, color = SNA_level)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.1) +
  facet_wrap(~module, scales = "free_y") +
  scale_color_manual(values = c("SNAhigh" = "#2166AC", "SNAlow" = "#D6604D")) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "MDSC Program Trajectories: Day 14 → Day 25",
       x = "Timepoint", y = "Mean Module Score", color = "SNA Level")

ggsave("Trajectory_MDSC_D14_D25.png", traj_mdsc, width = 12, height = 8, dpi = 300)

interaction_mdsc <- combined_mdsc %>%
  group_by(module) %>%
  do(tidy(lm(score ~ timepoint * SNA_level, data = .))) %>%
  filter(term == "timepointDay25:SNA_levelSNAlow")
print(interaction_mdsc %>% select(module, estimate, p.value))


# ============================================================ #
#   SECTION 3: EPITHELIAL CELLS
# ============================================================ #

# ------------------------------------------------------------ #
# 3A. DEFINE EPITHELIAL GENE MODULES
# ------------------------------------------------------------ #

# Epithelial identity / differentiation
Epi_Identity_genes <- c(
  "Epcam",
  "Krt8", "Krt18", "Krt19",
  "Cldn3", "Cldn4",
  "Cdh1",    # E-cadherin
  "Ocln",
  "Tjp1"
)

# EMT (Epithelial-to-Mesenchymal Transition)
Epi_EMT_genes <- c(
  "Vim",
  "Cdh2",    # N-cadherin
  "Fn1",
  "Snai1", "Snai2",
  "Twist1", "Twist2",
  "Zeb1", "Zeb2",
  "Mmp2", "Mmp9"
)

# Proliferation / Cell cycle
Epi_Proliferation_genes <- c(
  "Mki67",
  "Top2a",
  "Pcna",
  "Cdk1",
  "Ccnb1",
  "Ccnd1",
  "E2f1"
)

# Stress / Senescence
Epi_Stress_genes <- c(
  "Cdkn1a",  # p21
  "Cdkn2a",  # p16
  "Trp53",
  "Hmga1",
  "Sqstm1",
  "Hspa5",
  "Ddit3"
)

# Tumor / Oncogenic signaling
Epi_Tumor_genes <- c(
  "Myc",
  "Kras",
  "Egfr",
  "Erbb2",
  "Igf1r",
  "Birc5",   # Survivin
  "Vegfa",
  "Hif1a"
)

# Secretory / Inflammatory signaling
Epi_Secretory_genes <- c(
  "Il6",
  "Il8",
  "Cxcl1",
  "Cxcl5",
  "Csf2",
  "Tgfa",
  "Areg",
  "Ereg"
)

epi_modules <- list(
  Identity      = Epi_Identity_genes,
  EMT           = Epi_EMT_genes,
  Proliferation = Epi_Proliferation_genes,
  Stress        = Epi_Stress_genes,
  Tumor         = Epi_Tumor_genes,
  Secretory     = Epi_Secretory_genes
)

score_cols_epi <- c(
  "Identity_score", "EMT_score", "Proliferation_score",
  "Stress_score", "Tumor_score", "Secretory_score"
)

# ------------------------------------------------------------ #
# 3B. SUBSET EPITHELIAL CELLS
# ------------------------------------------------------------ #

epi_d14 <- subset(integrated, idents = "Epithelial")
epi_d25 <- subset(sna,        idents = "Epithelial")

# ------------------------------------------------------------ #
# 3C. DEG: SNAhigh vs SNAlow
# ------------------------------------------------------------ #
DefaultAssay(epi_d14) <- "RNA"

# Normalize (safe even if already done)
epi_d14 <- NormalizeData(epi_d14)
epi_d14 <- FindVariableFeatures(epi_d14)
epi_d14 <- ScaleData(epi_d14)

Idents(epi_d14) <- "condition"
deg_epi_d14 <- FindMarkers(
  epi_d14,
  ident.1 = "SNAhigh",
  ident.2 = "SNAlow",
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_epi_d14$gene <- rownames(deg_epi_d14)
write.csv(deg_epi_d14, "DEG_Epithelial_Day14_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

Idents(epi_d25) <- "SNA_level_D25"
deg_epi_d25 <- FindMarkers(epi_d25, ident.1 = "SNAhigh", ident.2 = "SNAlow",
                           min.pct = 0.25, logfc.threshold = 0.25)
deg_epi_d25$gene <- rownames(deg_epi_d25)
write.csv(deg_epi_d25, "DEG_Epithelial_Day25_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

volcano_epi_d14 <- plot_volcano(deg_epi_d14, "Epithelial — Day 14: SNAhigh vs SNAlow")
volcano_epi_d25 <- plot_volcano(deg_epi_d25, "Epithelial — Day 25: SNAhigh vs SNAlow")
ggsave("Volcano_Epithelial_D14.png", volcano_epi_d14, width = 8, height = 6, dpi = 300)
ggsave("Volcano_Epithelial_D25.png", volcano_epi_d25, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------ #
# 3D. MODULE SCORING — Epithelial
# ------------------------------------------------------------ #

epi_d14 <- AddModuleScore(epi_d14, features = epi_modules, name = names(epi_modules))
epi_d14$Identity_score      <- epi_d14$Identity1
epi_d14$EMT_score           <- epi_d14$EMT2
epi_d14$Proliferation_score <- epi_d14$Proliferation3
epi_d14$Stress_score        <- epi_d14$Stress4
epi_d14$Tumor_score         <- epi_d14$Tumor5
epi_d14$Secretory_score     <- epi_d14$Secretory6

epi_d25 <- AddModuleScore(epi_d25, features = epi_modules, name = names(epi_modules))
epi_d25$Identity_score      <- epi_d25$Identity1
epi_d25$EMT_score           <- epi_d25$EMT2
epi_d25$Proliferation_score <- epi_d25$Proliferation3
epi_d25$Stress_score        <- epi_d25$Stress4
epi_d25$Tumor_score         <- epi_d25$Tumor5
epi_d25$Secretory_score     <- epi_d25$Secretory6

vln_epi_d14 <- plot_module_violin(epi_d14@meta.data, "condition",
                                  score_cols_epi, c("SNAhigh","SNAlow"),
                                  condition_colors, "Epithelial Programs — Day 14")
vln_epi_d25 <- plot_module_violin(epi_d25@meta.data, "SNA_level_D25",
                                  score_cols_epi, c("SNAhigh","SNAlow"),
                                  condition_colors, "Epithelial Programs — Day 25")
ggsave("ModuleScore_Epithelial_D14.png", vln_epi_d14, width = 12, height = 8, dpi = 300)
ggsave("ModuleScore_Epithelial_D25.png", vln_epi_d25, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------ #
# 3E. TEMPORAL TRAJECTORY — Epithelial
# ------------------------------------------------------------ #

meta_epi_d14 <- epi_d14@meta.data %>%
  select(condition, all_of(score_cols_epi)) %>%
  dplyr::rename(SNA_level = condition) %>%
  mutate(timepoint = "Day14")

meta_epi_d25 <- epi_d25@meta.data %>%
  select(SNA_level_D25, all_of(score_cols_epi)) %>%
  dplyr::rename(SNA_level = SNA_level_D25) %>%
  mutate(timepoint = "Day25")

combined_epi <- bind_rows(meta_epi_d14, meta_epi_d25) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Day14", "Day25")),
    SNA_level = factor(SNA_level, levels = c("SNAhigh", "SNAlow"))
  ) %>%
  pivot_longer(cols = all_of(score_cols_epi),
               names_to = "module", values_to = "score")

summary_epi <- combined_epi %>%
  dplyr::group_by(module, SNA_level, timepoint) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE),
            se = sd(score, na.rm = TRUE) / sqrt(n()), .groups = "drop")

traj_epi <- ggplot(summary_epi,
                   aes(x = timepoint, y = mean_score,
                       group = SNA_level, color = SNA_level)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.1) +
  facet_wrap(~module, scales = "free_y") +
  scale_color_manual(values = c("SNAhigh" = "#2166AC", "SNAlow" = "#D6604D")) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Epithelial Program Trajectories: Day 14 → Day 25",
       x = "Timepoint", y = "Mean Module Score", color = "SNA Level")

ggsave("Trajectory_Epithelial_D14_D25.png", traj_epi, width = 14, height = 8, dpi = 300)

interaction_epi <- combined_epi %>%
  group_by(module) %>%
  do(tidy(lm(score ~ timepoint * SNA_level, data = .))) %>%
  filter(term == "timepointDay25:SNA_levelSNAlow")
print(interaction_epi %>% select(module, estimate, p.value))


# ============================================================ #
#   SECTION 4: FIBROBLASTS
#   Includes: Stromal_Fibroblasts + Proliferating_Fibroblasts (D14)
#             Fibroblasts (D25, D40)
# ============================================================ #

# ------------------------------------------------------------ #
# 4A. DEFINE FIBROBLAST GENE MODULES
# ------------------------------------------------------------ #

# Stromal / CAF identity
Fib_Stromal_genes <- c(
  "Col1a1", "Col1a2", "Col3a1",
  "Dcn", "Lum",
  "Pdgfra",
  "S100a4",
  "Fap",
  "Acta2",   # alpha-SMA
  "Thy1"     # CD90
)

# ECM remodeling
Fib_ECM_genes <- c(
  "Mmp2", "Mmp11", "Mmp14",
  "Timp1", "Timp2",
  "Lox", "Loxl2",
  "Fn1",
  "Thbs1", "Thbs2",
  "Postn"
)

# Inflammatory CAF (iCAF)
Fib_iCAF_genes <- c(
  "Il6",
  "Il1b",
  "Cxcl1", "Cxcl2", "Cxcl12",
  "Ccl2", "Ccl5",
  "Lif",
  "Tnfsf11",
  "Saa3"
)

# Myofibroblast / Contractile (myCAF)
Fib_myCAF_genes <- c(
  "Acta2",
  "Tagln",
  "Myl9",
  "Mylk",
  "Tpm1", "Tpm2",
  "Cnn1",
  "Des"
)

# Proliferation
Fib_Proliferation_genes <- c(
  "Mki67",
  "Top2a",
  "Pcna",
  "Cdk1",
  "Ccnb1",
  "Stmn1",
  "Birc5"
)

# Antigen presentation / immune crosstalk
Fib_ImmuneX_genes <- c(
  "H2-Aa", "H2-Ab1",
  "Cd74",
  "B2m",
  "Pdcd1lg2",  # PD-L2
  "Cd274",     # PD-L1
  "Icam1",
  "Vcam1"
)

fib_modules <- list(
  Stromal       = Fib_Stromal_genes,
  ECM           = Fib_ECM_genes,
  iCAF          = Fib_iCAF_genes,
  myCAF         = Fib_myCAF_genes,
  Proliferation = Fib_Proliferation_genes,
  ImmuneX       = Fib_ImmuneX_genes
)

score_cols_fib <- c(
  "Stromal_score", "ECM_score", "iCAF_score",
  "myCAF_score", "Proliferation_score", "ImmuneX_score"
)

# ------------------------------------------------------------ #
# 4B. SUBSET FIBROBLASTS
# Day 14 has Stromal_Fibroblasts + Proliferating_Fibroblasts
# ------------------------------------------------------------ #

fib_d14 <- subset(
  integrated,
  idents = c("Stromal_Fibroblasts", "Proliferating_Fibroblasts")
)

fib_d25 <- subset(sna, idents = "Fibroblasts")

# ------------------------------------------------------------ #
# 4C. DEG: SNAhigh vs SNAlow
# ------------------------------------------------------------ #
DefaultAssay(fib_d14) <- "RNA"

fib_d14 <- NormalizeData(fib_d14)
fib_d14 <- FindVariableFeatures(fib_d14)
fib_d14 <- ScaleData(fib_d14)
Idents(fib_d14) <- "condition"
deg_fib_d14 <- FindMarkers(fib_d14, ident.1 = "SNAhigh", ident.2 = "SNAlow",
                           min.pct = 0.25, logfc.threshold = 0.25)
deg_fib_d14$gene <- rownames(deg_fib_d14)
write.csv(deg_fib_d14, "DEG_Fibroblasts_Day14_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

Idents(fib_d25) <- "SNA_level_D25"
deg_fib_d25 <- FindMarkers(fib_d25, ident.1 = "SNAhigh", ident.2 = "SNAlow",
                           min.pct = 0.25, logfc.threshold = 0.25)
deg_fib_d25$gene <- rownames(deg_fib_d25)
write.csv(deg_fib_d25, "DEG_Fibroblasts_Day25_SNAhigh_vs_SNAlow.csv", row.names = TRUE)

volcano_fib_d14 <- plot_volcano(deg_fib_d14, "Fibroblasts — Day 14: SNAhigh vs SNAlow")
volcano_fib_d25 <- plot_volcano(deg_fib_d25, "Fibroblasts — Day 25: SNAhigh vs SNAlow")
ggsave("Volcano_Fibroblasts_D14.png", volcano_fib_d14, width = 8, height = 6, dpi = 300)
ggsave("Volcano_Fibroblasts_D25.png", volcano_fib_d25, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------ #
# 4D. MODULE SCORING — Fibroblasts
# ------------------------------------------------------------ #

fib_d14 <- AddModuleScore(fib_d14, features = fib_modules, name = names(fib_modules))
fib_d14$Stromal_score       <- fib_d14$Stromal1
fib_d14$ECM_score           <- fib_d14$ECM2
fib_d14$iCAF_score          <- fib_d14$iCAF3
fib_d14$myCAF_score         <- fib_d14$myCAF4
fib_d14$Proliferation_score <- fib_d14$Proliferation5
fib_d14$ImmuneX_score       <- fib_d14$ImmuneX6

fib_d25 <- AddModuleScore(fib_d25, features = fib_modules, name = names(fib_modules))
fib_d25$Stromal_score       <- fib_d25$Stromal1
fib_d25$ECM_score           <- fib_d25$ECM2
fib_d25$iCAF_score          <- fib_d25$iCAF3
fib_d25$myCAF_score         <- fib_d25$myCAF4
fib_d25$Proliferation_score <- fib_d25$Proliferation5
fib_d25$ImmuneX_score       <- fib_d25$ImmuneX6

vln_fib_d14 <- plot_module_violin(fib_d14@meta.data, "condition",
                                  score_cols_fib, c("SNAhigh","SNAlow"),
                                  condition_colors, "Fibroblast Programs — Day 14")
vln_fib_d25 <- plot_module_violin(fib_d25@meta.data, "SNA_level_D25",
                                  score_cols_fib, c("SNAhigh","SNAlow"),
                                  condition_colors, "Fibroblast Programs — Day 25")
ggsave("ModuleScore_Fibroblasts_D14.png", vln_fib_d14, width = 12, height = 8, dpi = 300)
ggsave("ModuleScore_Fibroblasts_D25.png", vln_fib_d25, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------ #
# 4E. TEMPORAL TRAJECTORY — Fibroblasts
# ------------------------------------------------------------ #

meta_fib_d14 <- fib_d14@meta.data %>%
  select(condition, all_of(score_cols_fib)) %>%
  dplyr::rename(SNA_level = condition) %>%
  mutate(timepoint = "Day14")

meta_fib_d25 <- fib_d25@meta.data %>%
  select(SNA_level_D25, all_of(score_cols_fib)) %>%
  dplyr::rename(SNA_level = SNA_level_D25) %>%
  mutate(timepoint = "Day25")

combined_fib <- bind_rows(meta_fib_d14, meta_fib_d25) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("Day14", "Day25")),
    SNA_level = factor(SNA_level, levels = c("SNAhigh", "SNAlow"))
  ) %>%
  pivot_longer(cols = all_of(score_cols_fib),
               names_to = "module", values_to = "score")

summary_fib <- combined_fib %>%
  dplyr::group_by(module, SNA_level, timepoint) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE),
            se = sd(score, na.rm = TRUE) / sqrt(n()), .groups = "drop")

traj_fib <- ggplot(summary_fib,
                   aes(x = timepoint, y = mean_score,
                       group = SNA_level, color = SNA_level)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.1) +
  facet_wrap(~module, scales = "free_y") +
  scale_color_manual(values = c("SNAhigh" = "#2166AC", "SNAlow" = "#D6604D")) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Fibroblast Program Trajectories: Day 14 → Day 25",
       x = "Timepoint", y = "Mean Module Score", color = "SNA Level")

ggsave("Trajectory_Fibroblasts_D14_D25.png", traj_fib, width = 14, height = 8, dpi = 300)

interaction_fib <- combined_fib %>%
  group_by(module) %>%
  do(tidy(lm(score ~ timepoint * SNA_level, data = .))) %>%
  filter(term == "timepointDay25:SNA_levelSNAlow")
print(interaction_fib %>% select(module, estimate, p.value))


# ============================================================ #
# END OF SCRIPT
# All output files saved to working directory:
#
# DEGs:
#   DEG_Tcells_Day14/25_SNAhigh_vs_SNAlow.csv
#   DEG_MDSC_Day14/25_SNAhigh_vs_SNAlow.csv
#   DEG_Epithelial_Day14/25_SNAhigh_vs_SNAlow.csv
#   DEG_Fibroblasts_Day14/25_SNAhigh_vs_SNAlow.csv
#
# Volcano plots:
#   Volcano_[CellType]_D14/25.png
#
# Module score violins (D14, D25):
#   ModuleScore_[CellType]_D14/25.png
#
# Temporal trajectories (D14 → D25):
#   Trajectory_[CellType]_D14_D25.png
#
# NOTE: Day 40 excluded from all analyses —
#       no SNAlow control available for valid comparison
# ============================================================ #