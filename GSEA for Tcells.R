# ============================================================ #
#   GSEA ANALYSIS — T CELLS
#   Per-subcluster ranked gene lists (SNAhigh vs SNAlow)
#   Day 14: T_cells_1, T_cells_2, T_cells_3 (3 subclusters)
#   Day 25: T_cells (1 cluster)
#   Dot plots, heatmaps, bar plots
# ============================================================ #

library(Seurat)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(tidyr)
library(tibble)
library(forcats)

# ============================================================ #
# STEP 1: GENE SETS
# ============================================================ #

hallmark_sets <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::mutate(gs_name = gsub("^HALLMARK_", "", gs_name)) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

gobp_sets <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::mutate(gs_name = gsub("^GOBP_", "", gs_name)) %>%
  dplyr::group_by(gs_name) %>%
  dplyr::filter(dplyr::n() >= 15 & dplyr::n() <= 500) %>%
  dplyr::ungroup() %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

kegg_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::mutate(gs_name = gsub("^KEGG_", "", gs_name)) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# ============================================================ #
# STEP 2: GSEA FUNCTION (Seurat v5 compatible)
# ============================================================ #

run_gsea_for_subcluster <- function(seurat_obj,
                                    cluster_id,
                                    condition_col,
                                    gene_sets,
                                    ident_1 = "SNAhigh",
                                    ident_2 = "SNAlow") {
  
  message("  Subsetting: ", cluster_id)
  sub_obj <- subset(seurat_obj, idents = cluster_id)
  
  if (ncol(sub_obj) < 20) {
    message("  Skipping ", cluster_id, " — too few cells")
    return(NULL)
  }
  
  meta <- sub_obj@meta.data
  
  if (!all(c(ident_1, ident_2) %in% meta[[condition_col]])) {
    message("  Skipping ", cluster_id, " — missing condition levels")
    return(NULL)
  }
  
  Idents(sub_obj) <- meta[[condition_col]]
  
  n_per_group <- table(Idents(sub_obj))
  if (any(n_per_group[c(ident_1, ident_2)] < 10)) {
    message("  Skipping ", cluster_id, " — <10 cells in one condition")
    return(NULL)
  }
  
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- JoinLayers(sub_obj, assay = "RNA")
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  
  deg <- tryCatch(
    FindMarkers(
      sub_obj,
      ident.1         = ident_1,
      ident.2         = ident_2,
      min.pct         = 0.1,
      logfc.threshold = 0,
      assay           = "RNA",
      verbose         = FALSE
    ),
    error = function(e) {
      message("  Error in DEG for ", cluster_id, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(deg) || nrow(deg) < 50) {
    message("  Skipping ", cluster_id, " — too few DEGs")
    return(NULL)
  }
  
  deg$gene <- rownames(deg)
  deg <- deg %>%
    dplyr::filter(!is.na(p_val) & p_val > 0) %>%
    dplyr::mutate(rank_metric = avg_log2FC * -log10(p_val)) %>%
    dplyr::arrange(desc(rank_metric))
  
  ranked_genes <- setNames(deg$rank_metric, deg$gene)
  
  set.seed(42)
  gsea_res <- fgsea(
    pathways    = gene_sets,
    stats       = ranked_genes,
    minSize     = 10,
    maxSize     = 500,
    nPermSimple = 10000
  )
  
  gsea_res$subcluster <- cluster_id
  return(as.data.frame(gsea_res))
}

# ============================================================ #
# STEP 3: RUN GSEA
# ============================================================ #

gene_set_collection <- hallmark_sets
collection_label    <- "Hallmark"

# --- Prepare integrated object ---
# Remove SCT if still present (comment out if already done)
# integrated[["SCT"]] <- NULL

DefaultAssay(integrated) <- "RNA"
integrated <- JoinLayers(integrated, assay = "RNA")
integrated <- NormalizeData(integrated, verbose = FALSE)

integrated$condition <- ifelse(integrated$sample == "855", "SNAhigh", "SNAlow")
Idents(integrated)   <- integrated$celltype

# --- Prepare sna object ---
Idents(sna) <- sna$cluster_id

# ------------------------------------------------------------ #
# Day 14 — T_cells_1, T_cells_2, T_cells_3
# ------------------------------------------------------------ #

tcell_clusters_d14 <- c("T_cells_1", "T_cells_2", "T_cells_3")

message("=== Running GSEA — Day 14 T cells ===")

gsea_tcell_d14_list <- lapply(tcell_clusters_d14, function(cl) {
  message("Processing: ", cl)
  run_gsea_for_subcluster(
    seurat_obj    = integrated,
    cluster_id    = cl,
    condition_col = "condition",
    gene_sets     = gene_set_collection
  )
})

gsea_tcell_d14 <- dplyr::bind_rows(gsea_tcell_d14_list)
gsea_tcell_d14$timepoint <- "Day14"

# ------------------------------------------------------------ #
# Day 25 — T_cells (single cluster)
# ------------------------------------------------------------ #

message("\n=== Running GSEA — Day 25 T cells ===")

gsea_tcell_d25 <- run_gsea_for_subcluster(
  seurat_obj    = sna,
  cluster_id    = "T_cells",
  condition_col = "SNA_level_D25",
  gene_sets     = gene_set_collection
)

if (!is.null(gsea_tcell_d25)) {
  gsea_tcell_d25$timepoint <- "Day25"
}

# ------------------------------------------------------------ #
# Combine and save
# ------------------------------------------------------------ #

gsea_tcell_all <- dplyr::bind_rows(gsea_tcell_d14, gsea_tcell_d25)

write.csv(gsea_tcell_all,
          paste0("GSEA_Tcells_", collection_label, "_all_results.csv"),
          row.names = FALSE)

# ============================================================ #
# STEP 4: DOT PLOT — Day 14 (3 subclusters on x-axis)
# ============================================================ #

plot_gsea_dotplot <- function(gsea_df,
                              timepoint_label,
                              cluster_order,
                              top_n = 15,
                              padj_cutoff = 0.25,
                              title = NULL) {
  
  df <- gsea_df %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    dplyr::mutate(
      direction      = ifelse(NES > 0, "Enriched in SNAhigh", "Enriched in SNAlow"),
      neg_log10_padj = -log10(padj),
      pathway_clean  = gsub("_", " ", pathway),
      subcluster     = factor(subcluster, levels = cluster_order)
    )
  
  top_pathways <- df %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::group_by(subcluster, direction) %>%
    dplyr::slice_min(order_by = padj, n = top_n, with_ties = FALSE) %>%
    dplyr::pull(pathway_clean) %>%
    unique()
  
  if (length(top_pathways) == 0) {
    message("No significant pathways at padj < ", padj_cutoff,
            " for ", timepoint_label)
    return(NULL)
  }
  
  plot_df <- df %>%
    dplyr::filter(pathway_clean %in% top_pathways)
  
  pathway_order <- plot_df %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_NES) %>%
    dplyr::pull(pathway_clean)
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = pathway_order)
  
  if (is.null(title)) {
    title <- paste0("GSEA: T Cells — ", timepoint_label,
                    " (", collection_label, ")")
  }
  
  p <- ggplot(plot_df, aes(x = subcluster, y = pathway_clean)) +
    geom_point(aes(size = neg_log10_padj, color = NES)) +
    scale_color_gradient2(
      low = "#D6604D", mid = "white", high = "#2166AC",
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(
      range = c(1, 8),
      name  = expression(-log[10](p[adj]))
    ) +
    facet_wrap(~direction, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y      = element_text(size = 9),
      strip.text        = element_text(face = "bold", size = 12),
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position   = "right",
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.3)
    ) +
    labs(title = title, x = "Subcluster", y = "Pathway")
  
  return(p)
}

gsea_dot_tcell_d14 <- plot_gsea_dotplot(
  gsea_tcell_all, "Day14",
  cluster_order = c("T_cells_1", "T_cells_2", "T_cells_3"),
  top_n = 15
)
if (!is.null(gsea_dot_tcell_d14)) {
  ggsave(paste0("GSEA_Dotplot_Tcell_D14_", collection_label, ".png"),
         gsea_dot_tcell_d14, width = 14, height = 10, dpi = 300)
}

# ============================================================ #
# STEP 5: HEATMAP — Day 14 (3 subclusters)
# ============================================================ #

plot_gsea_heatmap <- function(gsea_df,
                              timepoint_label,
                              cluster_order,
                              top_n = 25,
                              padj_cutoff = 0.25,
                              title = NULL) {
  
  df <- gsea_df %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    dplyr::mutate(
      pathway_clean = gsub("_", " ", pathway),
      NES_display   = NES,
      subcluster    = factor(subcluster, levels = cluster_order)
    )
  
  top_pathways <- df %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(
      max_abs_NES = max(abs(NES)),
      n_sig       = dplyr::n(),
      .groups     = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(max_abs_NES)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(pathway_clean)
  
  if (length(top_pathways) == 0) {
    message("No significant pathways for heatmap at ", timepoint_label)
    return(NULL)
  }
  
  plot_df <- df %>%
    dplyr::filter(pathway_clean %in% top_pathways)
  
  nes_mat <- plot_df %>%
    dplyr::select(pathway_clean, subcluster, NES) %>%
    tidyr::pivot_wider(names_from = subcluster, values_from = NES,
                       values_fill = 0) %>%
    tibble::column_to_rownames("pathway_clean")
  
  if (nrow(nes_mat) > 2) {
    pathway_order <- rownames(nes_mat)[hclust(dist(nes_mat))$order]
  } else {
    pathway_order <- rownames(nes_mat)
  }
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = pathway_order)
  
  if (is.null(title)) {
    title <- paste0("GSEA NES: T Cells — ", timepoint_label,
                    " (", collection_label, ")")
  }
  
  p <- ggplot(plot_df, aes(x = subcluster, y = pathway_clean, fill = NES_display)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, na.value = "grey90",
      name = "NES\n(SNAhigh \u2191 / SNAlow \u2193)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y  = element_text(size = 9),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid   = element_blank()
    ) +
    labs(title = title, x = "Subcluster", y = "Pathway")
  
  return(p)
}

heatmap_tcell_d14 <- plot_gsea_heatmap(
  gsea_tcell_all, "Day14",
  cluster_order = c("T_cells_1", "T_cells_2", "T_cells_3"),
  top_n = 25
)
if (!is.null(heatmap_tcell_d14)) {
  ggsave(paste0("GSEA_Heatmap_Tcell_D14_", collection_label, ".png"),
         heatmap_tcell_d14, width = 10, height = 12, dpi = 300)
}

# ============================================================ #
# STEP 6: Day 25 — BAR PLOT (single cluster)
# ============================================================ #

plot_gsea_barplot <- function(gsea_df,
                              timepoint_label,
                              top_n = 15,
                              padj_cutoff = 0.05) {
  
  df <- gsea_df %>%
    dplyr::filter(timepoint == timepoint_label, padj < padj_cutoff) %>%
    dplyr::mutate(pathway_clean = gsub("_", " ", pathway))
  
  if (nrow(df) == 0) {
    df <- gsea_df %>%
      dplyr::filter(timepoint == timepoint_label, padj < 0.25) %>%
      dplyr::mutate(pathway_clean = gsub("_", " ", pathway))
    if (nrow(df) == 0) return(NULL)
    padj_cutoff <- 0.25
  }
  
  top_up <- df %>%
    dplyr::filter(NES > 0) %>%
    dplyr::slice_max(order_by = NES, n = top_n, with_ties = FALSE)
  
  top_down <- df %>%
    dplyr::filter(NES < 0) %>%
    dplyr::slice_min(order_by = NES, n = top_n, with_ties = FALSE)
  
  bar_df <- dplyr::bind_rows(top_up, top_down) %>%
    dplyr::mutate(
      direction     = ifelse(NES > 0, "SNAhigh", "SNAlow"),
      pathway_clean = forcats::fct_reorder(pathway_clean, NES)
    )
  
  p <- ggplot(bar_df, aes(x = NES, y = pathway_clean, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("SNAhigh" = "#2166AC", "SNAlow" = "#D6604D")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y  = element_text(size = 9),
      legend.title = element_blank()
    ) +
    labs(
      title = paste0("GSEA: T Cells — ", timepoint_label,
                     " (", collection_label, ", padj < ", padj_cutoff, ")"),
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    )
  
  return(p)
}

bar_tcell_d25 <- plot_gsea_barplot(gsea_tcell_all, "Day25", top_n = 15)
if (!is.null(bar_tcell_d25)) {
  ggsave(paste0("GSEA_Barplot_Tcell_D25_", collection_label, ".png"),
         bar_tcell_d25, width = 10, height = 8, dpi = 300)
}

# ============================================================ #
# STEP 7: COMBINED VIEW — ALL 4 CLUSTERS ACROSS TIMEPOINTS
# x-axis: T_cells_1, T_cells_2, T_cells_3, T_cells_D25
# ============================================================ #

gsea_tcell_all_relabeled <- gsea_tcell_all %>%
  dplyr::mutate(
    subcluster_tp = case_when(
      timepoint == "Day25" ~ paste0(subcluster, "_D25"),
      TRUE                 ~ subcluster
    )
  )

combined_cluster_order <- c("T_cells_1", "T_cells_2", "T_cells_3", "T_cells_D25")

# Combined dot plot
plot_gsea_dotplot_combined <- function(gsea_df,
                                       cluster_order,
                                       top_n = 15,
                                       padj_cutoff = 0.25) {
  
  df <- gsea_df %>%
    dplyr::mutate(
      direction      = ifelse(NES > 0, "Enriched in SNAhigh", "Enriched in SNAlow"),
      neg_log10_padj = -log10(padj),
      pathway_clean  = gsub("_", " ", pathway),
      subcluster_tp  = factor(subcluster_tp, levels = cluster_order)
    )
  
  top_pathways <- df %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::group_by(subcluster_tp, direction) %>%
    dplyr::slice_min(order_by = padj, n = top_n, with_ties = FALSE) %>%
    dplyr::pull(pathway_clean) %>%
    unique()
  
  if (length(top_pathways) == 0) return(NULL)
  
  plot_df <- df %>%
    dplyr::filter(pathway_clean %in% top_pathways)
  
  pathway_order <- plot_df %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_NES) %>%
    dplyr::pull(pathway_clean)
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = pathway_order)
  
  p <- ggplot(plot_df, aes(x = subcluster_tp, y = pathway_clean)) +
    geom_point(aes(size = neg_log10_padj, color = NES)) +
    scale_color_gradient2(
      low = "#D6604D", mid = "white", high = "#2166AC",
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(range = c(1, 8),
                          name = expression(-log[10](p[adj]))) +
    facet_wrap(~direction, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y      = element_text(size = 9),
      strip.text        = element_text(face = "bold", size = 12),
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position   = "right",
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.3)
    ) +
    labs(
      title = paste0("GSEA: T Cells — D14 + D25 Combined (",
                     collection_label, ")"),
      x = "Subcluster (Timepoint)",
      y = "Pathway"
    )
  
  return(p)
}

gsea_dot_combined <- plot_gsea_dotplot_combined(
  gsea_tcell_all_relabeled, combined_cluster_order, top_n = 15
)
if (!is.null(gsea_dot_combined)) {
  ggsave(paste0("GSEA_Dotplot_Tcell_Combined_", collection_label, ".png"),
         gsea_dot_combined, width = 16, height = 12, dpi = 300)
}

# Combined heatmap
plot_gsea_heatmap_combined <- function(gsea_df,
                                       cluster_order,
                                       top_n = 25,
                                       padj_cutoff = 0.25) {
  
  df <- gsea_df %>%
    dplyr::mutate(
      pathway_clean = gsub("_", " ", pathway),
      NES_display   = NES,
      subcluster_tp = factor(subcluster_tp, levels = cluster_order)
    )
  
  top_pathways <- df %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(
      max_abs_NES = max(abs(NES)),
      n_sig       = dplyr::n(),
      .groups     = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(max_abs_NES)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(pathway_clean)
  
  if (length(top_pathways) == 0) return(NULL)
  
  plot_df <- df %>%
    dplyr::filter(pathway_clean %in% top_pathways)
  
  nes_mat <- plot_df %>%
    dplyr::select(pathway_clean, subcluster_tp, NES) %>%
    tidyr::pivot_wider(names_from = subcluster_tp, values_from = NES,
                       values_fill = 0) %>%
    tibble::column_to_rownames("pathway_clean")
  
  if (nrow(nes_mat) > 2) {
    pathway_order <- rownames(nes_mat)[hclust(dist(nes_mat))$order]
  } else {
    pathway_order <- rownames(nes_mat)
  }
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = pathway_order)
  
  p <- ggplot(plot_df, aes(x = subcluster_tp, y = pathway_clean,
                           fill = NES_display)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, na.value = "grey90",
      name = "NES\n(SNAhigh \u2191 / SNAlow \u2193)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y  = element_text(size = 9),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid   = element_blank()
    ) +
    labs(
      title = paste0("GSEA NES: T Cells — D14 + D25 Combined (",
                     collection_label, ")"),
      x = "Subcluster (Timepoint)",
      y = "Pathway"
    )
  
  return(p)
}

heatmap_combined <- plot_gsea_heatmap_combined(
  gsea_tcell_all_relabeled, combined_cluster_order, top_n = 25
)
if (!is.null(heatmap_combined)) {
  ggsave(paste0("GSEA_Heatmap_Tcell_Combined_", collection_label, ".png"),
         heatmap_combined, width = 12, height = 12, dpi = 300)
}

# ============================================================ #
# OUTPUT FILES:
#   GSEA_Tcells_Hallmark_all_results.csv
#   GSEA_Dotplot_Tcell_D14_Hallmark.png       (3 subclusters)
#   GSEA_Heatmap_Tcell_D14_Hallmark.png       (3 subclusters)
#   GSEA_Barplot_Tcell_D25_Hallmark.png       (single cluster)
#   GSEA_Dotplot_Tcell_Combined_Hallmark.png  (all 4 together)
#   GSEA_Heatmap_Tcell_Combined_Hallmark.png  (all 4 together)
#
# To switch gene sets, change these and re-run from Step 3:
#   gene_set_collection <- gobp_sets
#   collection_label    <- "GO_BP"
# ============================================================ #