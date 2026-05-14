# ============================================================
# T-cell subset analysis - Day 14
# Keep original UMAP, split all views by sample (855 vs 857)
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ---- 1. Subset T cells ------------------------------------------------------
Idents(integrated) <- "celltype"   # or whichever column holds T_cells_1/2/3

t_cells <- subset(integrated, idents = c("T_cells_1", "T_cells_2", "T_cells_3"))

# Always plot expression from RNA, not the integrated assay
DefaultAssay(t_cells) <- "RNA"
t_cells <- NormalizeData(t_cells)   # safe to re-run; ensures RNA data slot is filled

# Sanity check
table(Idents(t_cells), t_cells$sample)

# ---- 2. T-cell UMAP, split by sample ---------------------------------------
p_umap <- DimPlot(
  t_cells, reduction = "umap",
  split.by = "sample", label = TRUE, repel = TRUE
) + ggtitle("T cells - 855 vs 857")

# ---- 3. Confirm T-cell identity: Cd3d, Cd3e, Cd3g (split by sample) --------
cd3_markers <- c("Cd3d", "Cd3e", "Cd3g")

p_cd3_feat <- FeaturePlot(
  t_cells, features = cd3_markers,
  split.by = "sample", order = TRUE,
  min.cutoff = "q10", max.cutoff = "q90"
) & theme(legend.position = "right")

# Violin: group by cluster, split by sample so both samples sit side-by-side
p_cd3_vln <- VlnPlot(
  t_cells, features = cd3_markers,
  split.by = "sample", pt.size = 0, ncol = 3
)

# ---- 4. CD4 (helper) vs CD8 (cytotoxic), split by sample -------------------
p_cd4_feat <- FeaturePlot(
  t_cells, features = "Cd4",
  split.by = "sample", order = TRUE,
  min.cutoff = "q10", max.cutoff = "q90"
)

p_cd8_feat <- FeaturePlot(
  t_cells, features = c("Cd8a", "Cd8b1"),
  split.by = "sample", order = TRUE,
  min.cutoff = "q10", max.cutoff = "q90"
)

p_cd4_cd8_vln <- VlnPlot(
  t_cells, features = c("Cd4", "Cd8a", "Cd8b1"),
  split.by = "sample", pt.size = 0, ncol = 3
)

# ---- 5. Foxp3 / Treg detection, split by sample ----------------------------
p_foxp3_feat <- FeaturePlot(
  t_cells, features = "Foxp3",
  split.by = "sample", order = TRUE,
  min.cutoff = "q10"
)

p_foxp3_vln <- VlnPlot(
  t_cells, features = "Foxp3",
  split.by = "sample", pt.size = 0
)

# Note on blend: Seurat's blend = TRUE doesn't combine with split.by.
# Workaround: split the object by sample, run blended FeaturePlot on each,
# then stitch them together with patchwork so 855 and 857 sit side-by-side.
t_855 <- subset(t_cells, subset = sample == "855")
t_857 <- subset(t_cells, subset = sample == "857")

p_cd4_foxp3 <-
  (FeaturePlot(t_855, features = c("Cd4", "Foxp3"),
               blend = TRUE, order = TRUE,
               cols = c("lightgrey", "red", "blue")) +
     plot_annotation(title = "855: Cd4 x Foxp3")) /
  (FeaturePlot(t_857, features = c("Cd4", "Foxp3"),
               blend = TRUE, order = TRUE,
               cols = c("lightgrey", "red", "blue")) +
     plot_annotation(title = "857: Cd4 x Foxp3"))

p_cd8_foxp3 <-
  (FeaturePlot(t_855, features = c("Cd8a", "Foxp3"),
               blend = TRUE, order = TRUE,
               cols = c("lightgrey", "red", "blue")) +
     plot_annotation(title = "855: Cd8a x Foxp3")) /
  (FeaturePlot(t_857, features = c("Cd8a", "Foxp3"),
               blend = TRUE, order = TRUE,
               cols = c("lightgrey", "red", "blue")) +
     plot_annotation(title = "857: Cd8a x Foxp3"))

# ---- 6. Quantify lineages per cluster, per sample --------------------------
expr <- FetchData(
  t_cells,
  vars = c("Cd4", "Cd8a", "Cd8b1", "Foxp3", "seurat_clusters", "sample")
)

thr <- 0   # raise (e.g. 0.5 or 1) if violins show high background

expr <- expr %>%
  mutate(
    is_CD4  = Cd4  > thr,
    is_CD8  = (Cd8a > thr) | (Cd8b1 > thr),
    is_Treg = Cd4 > thr & Foxp3 > thr,
    lineage = case_when(
      is_Treg                  ~ "CD4_Treg",
      is_CD4 & !is_CD8         ~ "CD4_conv",
      is_CD8 & !is_CD4         ~ "CD8",
      is_CD4 &  is_CD8         ~ "DP",
      !is_CD4 & !is_CD8         ~ "DN",
      TRUE                     ~ "other"
    )
  )

lineage_summary <- expr %>%
  count(sample, seurat_clusters, lineage) %>%
  group_by(sample, seurat_clusters) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

print(lineage_summary)

# Side-by-side stacked bars: clusters on x, faceted by sample
p_lineage_bar <- ggplot(lineage_summary,
                        aes(x = seurat_clusters, y = pct, fill = lineage)) +
  geom_col() +
  facet_wrap(~ sample) +
  labs(x = "T-cell cluster", y = "% of cells", fill = "Lineage call",
       title = "T-cell lineage composition: 855 vs 857") +
  theme_classic(base_size = 12)

# Same data, grouped by lineage with samples side-by-side for direct comparison
p_lineage_grouped <- ggplot(lineage_summary,
                            aes(x = lineage, y = pct, fill = sample)) +
  geom_col(position = "dodge") +
  facet_wrap(~ seurat_clusters) +
  labs(x = "Lineage call", y = "% of cells in cluster",
       title = "Lineage frequency by sample, per cluster") +
  theme_classic(base_size = 12)

##Removing DN adn DP

lineage_summary_clean <- lineage_summary %>%
  filter(!lineage %in% c("DN", "DP"))

ggplot(lineage_summary_clean,
       aes(x = seurat_clusters, y = pct, fill = lineage)) +
  geom_col() +
  facet_wrap(~ sample) +
  labs(x = "T-cell cluster", y = "% of cells", fill = "Lineage call",
       title = "T-cell lineage composition: 855 vs 857") +
  theme_classic(base_size = 12)


lineage_summary_clean <- expr %>%
  filter(!lineage %in% c("DN", "DP")) %>%
  count(sample, seurat_clusters, lineage) %>%
  group_by(sample, seurat_clusters) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# ---- 7. View / save --------------------------------------------------------
p_umap
p_cd3_feat
p_cd3_vln
p_cd4_feat
p_cd8_feat
p_cd4_cd8_vln
p_foxp3_feat
p_foxp3_vln
p_cd4_foxp3
p_cd8_foxp3
p_lineage_bar
p_lineage_grouped

# ggsave("Tcells_UMAP_split.png",   p_umap,            width = 10, height = 5, dpi = 300)
# ggsave("Tcells_CD3_split.png",    p_cd3_feat,        width = 10, height = 9, dpi = 300)
# ggsave("Tcells_CD4_split.png",    p_cd4_feat,        width = 10, height = 4, dpi = 300)
# ggsave("Tcells_CD8_split.png",    p_cd8_feat,        width = 10, height = 7, dpi = 300)
# ggsave("Tcells_Foxp3_split.png",  p_foxp3_feat,      width = 10, height = 4, dpi = 300)
# ggsave("Tcells_CD4_Foxp3.png",    p_cd4_foxp3,       width = 12, height = 8, dpi = 300)
# ggsave("Tcells_CD8_Foxp3.png",    p_cd8_foxp3,       width = 12, height = 8, dpi = 300)
# ggsave("Tcells_lineage_bar.png",  p_lineage_bar,     width = 10, height = 5, dpi = 300)
# ggsave("Tcells_lineage_grouped.png", p_lineage_grouped, width = 11, height = 6, dpi = 300)