##Final figures for SNAhigh and SNAlow paper
##SNAhigh Day 40 ##3 sites
DimPlot(snahigh, label = TRUE, split.by = "organ")

##SNAhigh vs low #Day 25
DimPlot(sna, reduction = "umap", split.by = "SNA_level_D25", label = TRUE, label.size = 4)

#SNAhigh vs SNAlow # Day 14
DimPlot(integrated, reduction = "umap", split.by = "sample", label = TRUE, repel = TRUE, ncol = 2)


##QC metrics check
# SNAhigh Day 40
snahigh[["percent.mt"]] <- PercentageFeatureSet(snahigh, pattern = "^mt-")

# Day 25
sna[["percent.mt"]] <- PercentageFeatureSet(sna, pattern = "^mt-")

# Day 14
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^mt-")

library(Seurat)
library(ggplot2)
library(patchwork)

qc_plot <- function(seurat_obj, group_var = NULL, title = "QC Metrics") {
  
  # Violin plots
  vln <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = group_var,
    pt.size = 0.1,
    ncol = 3
  ) + NoLegend()
  
  # Scatter: nCount vs nFeature
  scatter1 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
  
  # Scatter: nCount vs percent.mt
  scatter2 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )
  
  # Combine
  combined <- (vln / (scatter1 | scatter2)) +
    plot_annotation(title = title)
  
  return(combined)
}

qc_snahigh <- qc_plot(
  snahigh,
  group_var = "organ",
  title = "SNAhigh Day 40 QC (by organ)"
)

qc_snahigh

qc_sna25 <- qc_plot(
  sna,
  group_var = "SNA_level_D25",
  title = "Day 25 QC (SNAhigh vs SNAlow)"
)

qc_sna25

qc_sna14 <- qc_plot(
  integrated,
  group_var = "sample",
  title = "Day 14 QC (Sample-level)"
)

qc_sna14

##Post applyign QC metrics
# Day 40
snahigh_filtered <- subset(
  snahigh,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    nCount_RNA < 25000 &
    percent.mt < 5
)

# Day 25
sna_filtered <- subset(
  sna,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    nCount_RNA < 25000 &
    percent.mt < 5
)

# Day 14
integrated_filtered <- subset(
  integrated,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    nCount_RNA < 25000 &
    percent.mt < 5
)

library(Seurat)
library(ggplot2)
library(patchwork)

qc_plot <- function(seurat_obj, group_var = NULL, title = "") {
  
  vln <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = group_var,
    pt.size = 0.1,
    ncol = 3
  ) + NoLegend()
  
  scatter1 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
  
  scatter2 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )
  
  combined <- (vln / (scatter1 | scatter2)) +
    plot_annotation(title = title)
  
  return(combined)
}

qc_day14_post <- qc_plot(
  integrated_filtered,
  group_var = "sample",
  title = "Day 14 QC (Post-filtering)"
)

qc_day25_post <- qc_plot(
  sna_filtered,
  group_var = "SNA_level_D25",
  title = "Day 25 QC (Post-filtering)"
)

qc_day40_post <- qc_plot(
  snahigh_filtered,
  group_var = "organ",
  title = "Day 40 QC (Post-filtering)"
)

ggsave("QC_D14_post.png", qc_day14_post, width = 12, height = 8, dpi = 300)
ggsave("QC_D25_post.png", qc_day25_post, width = 12, height = 8, dpi = 300)
ggsave("QC_D40_post.png", qc_day40_post, width = 12, height = 8, dpi = 300)

##UMAPs
library(ggplot2)

# consistent colors (important: must match your identities)
celltype_colors <- c(
  "B_Cells_1" = "#F8766D",
  "B_Cells_2" = "#D89000",
  "B_cells_1" = "#D89000",
  "B_cells_2" = "#F8766D",
  "B_cells_3" = "#A3A500",
  "T_cell" = "#B79F00",
  "T_cells" = "#00BFC4",
  "T_cells_1" = "#00BFC4",
  "T_cells_2" = "#C77CFF",
  "NK_Cells" = "#FF61C3",
  "Macrophage_1" = "#7CAE00",
  "Macrophage_2" = "#00BA38",
  "Macrophage_3" = "#00A9FF",
  "Macrophage_4" = "#C77CFF",
  "MDSC" = "#00C08B",
  "Fibroblasts" = "#00B0F6",
  "Stromal_Fibroblasts" = "#9590FF",
  "Proliferating_Fibroblasts" = "#39B600",
  "Epithelial" = "#F564E3",
  "epithelial cells" = "#F564E3",
  "Endothelial" = "#FF66B2",
  "endothelial cells" = "#FF66B2",
  "Plasma_cells" = "#00BF7D",
  "DC_1" = "#00B4D8",
  "DC_2" = "#F15BB5",
  "pancreatic cells" = "#C77CFF"
)

# universal theme
umap_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    strip.text = element_text(size = 16, face = "bold") 
  )

#Day 40 UMAP
umap_day40 <- DimPlot(
  snahigh,
  reduction = "umap",
  split.by = "organ",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4,
  label.size = 5,
  cols = celltype_colors
) +
  ggtitle("Day 40: SNAhigh across organs") +
  umap_theme


##Day 25
umap_day25 <- DimPlot(
  sna,
  reduction = "umap",
  split.by = "SNA_level_D25",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5,
  label.size = 4,
  cols = celltype_colors
) +
  ggtitle("Day 25: SNAhigh vs SNAlow") +
  umap_theme

umap_day14 <- DimPlot(
  integrated,
  reduction = "umap",
  split.by = "sample",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5,
  label.size = 4,
  cols = celltype_colors,
  ncol = 2
) +
  ggtitle("Day 14: SNAhigh vs SNAlow") +
  umap_theme


###Frequency distribution Tables
library(dplyr)
library(ggplot2)

plot_composition <- function(seurat_obj, group_var, cluster_var = "ident", title = "") {
  
  df <- seurat_obj@meta.data
  
  # Add cluster info
  df$cluster <- Idents(seurat_obj)
  
  # Compute proportions
  comp <- df %>%
    group_by(.data[[group_var]], cluster) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[group_var]]) %>%
    dplyr::mutate(percent = (n / sum(n)) * 100)
  
  # Plot
  p <- ggplot(comp, aes(x = cluster, y = percent, fill = .data[[group_var]])) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = "Cell Type",
      y = "Percentage (%)",
      title = title,
      fill = NULL
    )
  
  return(p)
}

condition_colors <- c(
  "SNAhigh" = "blue",  # blue
  "SNAlow" = "red"    # red
)

#Day 25
comp_day25 <- plot_composition(
  sna,
  group_var = "SNA_level_D25",
  title = "Day 25: Cell Composition"
) +
  scale_fill_manual(values = condition_colors)

##Day 14
# adjust based on your mapping
integrated$condition <- ifelse(integrated$sample == "855", "SNAhigh", "SNAlow")

comp_day14 <- plot_composition(
  integrated,
  group_var = "condition",
  title = "Day 14: Cell Composition"
) +
  scale_fill_manual(values = condition_colors)

##Day 40
organ_colors <- c(
  "mesentery" = "#08306B",  # dark navy
  "omentum"   = "#2171B5",  # medium blue
  "ovary"     = "#6BAED6"   # light blue
)

comp_day40 <- plot_composition(
  snahigh,
  group_var = "organ",
  title = "Day 40: Cell Composition across organs"
) +
  scale_fill_manual(values = organ_colors)

##DotPlot with top 5 marker genes
library(Seurat)
library(dplyr)

get_top_markers <- function(seurat_obj) {
  
  markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.25
  )
  
  top5 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5)
  
  return(top5)
}

markers_day14 <- get_top_markers(integrated)
markers_day25 <- get_top_markers(sna)
markers_day40 <- get_top_markers(snahigh)

genes_day14 <- unique(markers_day14$gene)
genes_day25 <- unique(markers_day25$gene)
genes_day40 <- unique(markers_day40$gene)

dot_day14 <- DotPlot(
  integrated,
  features = genes_day14
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#0072B2") +
  ggtitle("Day 14: Top 5 markers per cluster") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dot_day25 <- DotPlot(
  sna,
  features = genes_day25
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#0072B2") +
  ggtitle("Day 25: Top 5 markers per cluster") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dot_day40 <- DotPlot(
  snahigh,
  features = genes_day40
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#0072B2") +
  ggtitle("Day 40: Top 5 markers per cluster") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("DotPlot_Day14.png", dot_day14, width = 14, height = 6, dpi = 300)
ggsave("DotPlot_Day25.png", dot_day25, width = 14, height = 6, dpi = 300)
ggsave("DotPlot_Day40.png", dot_day40, width = 14, height = 6, dpi = 300)
