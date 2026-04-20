# ============================================================ #
#   GSEA ANALYSIS — MACROPHAGE SUBCLUSTERS
#   Per-subcluster ranked gene lists (SNAhigh vs SNAlow)
#   Dot plots: pathways (y) × subclusters (x), split by SNA level
#   Day 14 (3 subclusters) and Day 25 (4 subclusters)
# ============================================================ #

library(Seurat)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(tidyr)
library(tibble)

# ============================================================ #
# STEP 1: RETRIEVE GENE SETS (MSigDB — MOUSE NATIVE)
# ============================================================ #

# --- Hallmark ---
hallmark_sets <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  # Clean up names: remove "HALLMARK_" prefix, replace underscores
  dplyr::mutate(gs_name = gsub("^HALLMARK_", "", gs_name)) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# --- GO BP (filtered to reasonably sized sets) ---
gobp_sets <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::mutate(gs_name = gsub("^GOBP_", "", gs_name)) %>%
  dplyr::group_by(gs_name) %>%
  dplyr::filter(dplyr::n() >= 15 & dplyr::n() <= 500) %>%
  dplyr::ungroup() %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# --- KEGG ---
kegg_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::mutate(gs_name = gsub("^KEGG_", "", gs_name)) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# ============================================================ #
# STEP 2: FIXED FUNCTION (ONLY CHANGE IS HERE)
# ============================================================ #

run_gsea_for_subcluster <- function(seurat_obj,
                                    cluster_id,
                                    condition_col,
                                    gene_sets,
                                    ident_1 = "SNAhigh",
                                    ident_2 = "SNAlow") {
  
  # Subset using Idents (fix)
  sub_obj <- subset(seurat_obj, idents = cluster_id)
  
  if (ncol(sub_obj) < 20) {
    message("  Skipping ", cluster_id, " — too few cells")
    return(NULL)
  }
  
  meta <- sub_obj@meta.data
  
  # Check conditions
  if (!all(c(ident_1, ident_2) %in% meta[[condition_col]])) {
    message("  Skipping ", cluster_id, " — missing condition levels")
    return(NULL)
  }
  
  # Assign identities by condition
  Idents(sub_obj) <- meta[[condition_col]]
  
  # Check cell counts
  n_per_group <- table(Idents(sub_obj))
  if (any(n_per_group[c(ident_1, ident_2)] < 10)) {
    message("  Skipping ", cluster_id, " — <10 cells in one condition")
    return(NULL)
  }
  
  # DEG
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  
  deg <- tryCatch(
    FindMarkers(
      sub_obj,
      ident.1        = ident_1,
      ident.2        = ident_2,
      min.pct        = 0.1,
      logfc.threshold = 0,
      assay          = "RNA",
      verbose        = FALSE
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
  
  # Ranking
  deg$gene <- rownames(deg)
  
  deg <- deg %>%
    dplyr::filter(!is.na(p_val) & p_val > 0) %>%
    dplyr::mutate(rank_metric = avg_log2FC * -log10(p_val)) %>%
    dplyr::arrange(desc(rank_metric))
  
  ranked_genes <- setNames(deg$rank_metric, deg$gene)
  
  # FGSEA
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

mac_clusters_d14 <- c("Macrophage_1", "Macrophage_2", "Macrophage_3")

# Condition
integrated$condition <- ifelse(integrated$sample == "855", "SNAhigh", "SNAlow")

# CRITICAL
Idents(integrated) <- integrated$cluster_id

message("=== Running GSEA — Day 14 Macrophages ===")

gsea_d14_list <- lapply(mac_clusters_d14, function(cl) {
  message("Processing: ", cl)
  
  run_gsea_for_subcluster(
    seurat_obj    = integrated,
    cluster_id    = cl,
    condition_col = "condition",
    gene_sets     = gene_set_collection
  )
})

gsea_d14 <- dplyr::bind_rows(gsea_d14_list)
gsea_d14$timepoint <- "Day14"

# ------------------------------------------------------------ #
# Day 25 — 4 macrophage subclusters
# Adjust these names to match your sna object's cluster labels
# ------------------------------------------------------------ #

# Auto-detect macrophage clusters from D25 object
sna$cluster_id <- Idents(sna)
Idents(sna) <- sna$cluster_id   # CRITICAL

all_clusters_d25 <- levels(Idents(sna))

mac_clusters_d25 <- all_clusters_d25[
  grepl("Macrophage|Mac", all_clusters_d25, ignore.case = TRUE)
]

message("\nDetected D25 macrophage clusters: ", 
        paste(mac_clusters_d25, collapse = ", "))

message("\n=== Running GSEA — Day 25 Macrophages ===")

gsea_d25_list <- lapply(mac_clusters_d25, function(cl) {
  message("Processing: ", cl)
  
  run_gsea_for_subcluster(
    seurat_obj    = sna,
    cluster_id    = cl,
    condition_col = "SNA_level_D25",
    gene_sets     = gene_set_collection
  )
})

gsea_d25 <- dplyr::bind_rows(gsea_d25_list)
gsea_d25$timepoint <- "Day25"

# ============================================================ #
# STEP 4: PREPARE COMBINED DATA FOR PLOTTING
# ============================================================ #

gsea_all <- dplyr::bind_rows(gsea_d14, gsea_d25)

# Save full results
write.csv(gsea_all, paste0("GSEA_Macrophages_", collection_label, "_all_results.csv"),
          row.names = FALSE)

# ============================================================ #
# STEP 5: DOT PLOT — Pathways (y) × Subclusters (x)
#   Separate panels for SNAhigh-enriched vs SNAlow-enriched
#   Color = NES, Size = -log10(padj)
# ============================================================ #

plot_gsea_dotplot <- function(gsea_df,
                              timepoint_label,
                              top_n = 15,
                              padj_cutoff = 0.25,
                              title_suffix = "") {
  
  df <- gsea_df %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    dplyr::mutate(
      direction = ifelse(NES > 0, "Enriched in SNAhigh", "Enriched in SNAlow"),
      neg_log10_padj = -log10(padj),
      pathway_clean  = gsub("_", " ", pathway)
    )
  
  # Select top pathways per subcluster per direction (by significance)
  top_pathways <- df %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::group_by(subcluster, direction) %>%
    dplyr::slice_min(order_by = padj, n = top_n, with_ties = FALSE) %>%
    dplyr::pull(pathway_clean) %>%
    unique()
  
  if (length(top_pathways) == 0) {
    message("No significant pathways at padj < ", padj_cutoff, " for ", timepoint_label)
    return(NULL)
  }
  
  plot_df <- df %>%
    dplyr::filter(pathway_clean %in% top_pathways)
  
  # Order pathways by mean NES across subclusters
  pathway_order <- plot_df %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_NES) %>%
    dplyr::pull(pathway_clean)
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = pathway_order)
  
  # Dot plot
  p <- ggplot(plot_df, aes(x = subcluster, y = pathway_clean)) +
    geom_point(aes(size = neg_log10_padj, color = NES)) +
    scale_color_gradient2(
      low      = "#D6604D",   # SNAlow direction (negative NES)
      mid      = "white",
      high     = "#2166AC",   # SNAhigh direction (positive NES)
      midpoint = 0,
      name     = "NES"
    ) +
    scale_size_continuous(
      range = c(1, 8),
      name  = expression(-log[10](p[adj]))
    ) +
    facet_wrap(~direction, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y     = element_text(size = 9),
      strip.text       = element_text(face = "bold", size = 12),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position  = "right",
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
    ) +
    labs(
      title = paste0("GSEA: Macrophage Subclusters — ", timepoint_label,
                     " (", collection_label, ")", title_suffix),
      x = "Macrophage Subcluster",
      y = "Pathway"
    )
  
  return(p)
}

# ------------------------------------------------------------ #
# Generate plots
# ------------------------------------------------------------ #

# Day 14
gsea_dot_d14 <- plot_gsea_dotplot(gsea_all, "Day14", top_n = 15)
if (!is.null(gsea_dot_d14)) {
  ggsave(paste0("GSEA_Dotplot_Macrophage_D14_", collection_label, ".png"),
         gsea_dot_d14, width = 14, height = 10, dpi = 300)
}

# Day 25
gsea_dot_d25 <- plot_gsea_dotplot(gsea_all, "Day25", top_n = 15)
if (!is.null(gsea_dot_d25)) {
  ggsave(paste0("GSEA_Dotplot_Macrophage_D25_", collection_label, ".png"),
         gsea_dot_d25, width = 16, height = 10, dpi = 300)
}

##Combined dotplot
# ============================================================ #
# STEP 6: COMBINED HEATMAP-STYLE NES PLOT (ALTERNATIVE VIEW)
#   Pathways (y) × Subclusters (x), tile colored by NES
#   Non-significant tiles greyed out
# ============================================================ #

plot_gsea_heatmap <- function(gsea_df,
                              timepoint_label,
                              top_n = 20,
                              padj_cutoff = 0.25) {
  
  df <- gsea_df %>%
    dplyr::filter(timepoint == timepoint_label) %>%
    dplyr::mutate(
      pathway_clean = gsub("_", " ", pathway),
      NES_display   =  NES 
    )
  
  # Top pathways (same logic)
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
  
  # Order pathways via clustering
  nes_mat <- plot_df %>%
    dplyr::select(pathway_clean, subcluster, NES) %>%
    tidyr::pivot_wider(names_from = subcluster, values_from = NES, values_fill = 0) %>%
    tibble::column_to_rownames("pathway_clean")
  
  if (nrow(nes_mat) > 2) {
    pathway_order <- rownames(nes_mat)[hclust(dist(nes_mat))$order]
  } else {
    pathway_order <- rownames(nes_mat)
  }
  
  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = pathway_order)
  
  # ======================== #
  # FINAL PLOT (CLEAN)
  # ======================== #
  
  p <- ggplot(plot_df, aes(x = subcluster, y = pathway_clean, fill = NES_display)) +
    
    geom_tile(color = "white", linewidth = 0.5) +
    
    scale_fill_gradient2(
      low      = "#B2182B",   # SNAlow
      mid      = "white",
      high     = "#2166AC",   # SNAhigh
      midpoint = 0,
      na.value = "grey90",
      name     = "NES\n(SNAhigh ↑ / SNAlow ↓)"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y  = element_text(size = 9),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid   = element_blank()
    ) +
    
    labs(
      title = paste0("GSEA NES: Macrophage Subclusters — ", timepoint_label,
                     " (", collection_label, ")"),
      x     = "Macrophage Subcluster",
      y     = "Pathway"
    )
  
  return(p)
}

heatmap_d14 <- plot_gsea_heatmap(gsea_all, "Day14", top_n = 25)
if (!is.null(heatmap_d14)) {
  ggsave(paste0("GSEA_Heatmap_Macrophage_D14_", collection_label, ".png"),
         heatmap_d14, width = 10, height = 12, dpi = 300)
}

heatmap_d25 <- plot_gsea_heatmap(gsea_all, "Day25", top_n = 25)
if (!is.null(heatmap_d25)) {
  ggsave(paste0("GSEA_Heatmap_Macrophage_D25_", collection_label, ".png"),
         heatmap_d25, width = 12, height = 12, dpi = 300)
}


# ============================================================ #
# STEP 7 (OPTIONAL): RE-RUN WITH DIFFERENT GENE SET COLLECTIONS
# ============================================================ #
# To switch collections, change these two lines and re-run from Step 3:
#
#   gene_set_collection <- gobp_sets
#   collection_label    <- "GO_BP"
#
# or:
#   gene_set_collection <- kegg_sets
#   collection_label    <- "KEGG"
#
# All filenames auto-update with the collection label.


# ============================================================ #
# OUTPUT FILES:
#   GSEA_Macrophages_[Collection]_all_results.csv
#   GSEA_Dotplot_Macrophage_D14/D25_[Collection].png
#   GSEA_Heatmap_Macrophage_D14/D25_[Collection].png
# ============================================================ #