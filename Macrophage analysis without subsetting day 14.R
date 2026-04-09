##MAcrophage analysis- top 20 genes for each type
library(dplyr)

# ----------------------------- #
# MACROPHAGE 1 (855 vs 857)
# ----------------------------- #
macro1 <- subset(integrated, idents = "Macrophage_1")

# IMPORTANT: re-run SCT prep AFTER subsetting
macro1 <- PrepSCTFindMarkers(macro1)

Idents(macro1) <- macro1$sample

deg_macro1 <- FindMarkers(
  macro1,
  ident.1 = "855",
  ident.2 = "857",
  logfc.threshold = 0,
  min.pct = 0.1
)

top_macro1 <- list(
  top_855 = deg_macro1 %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 20),
  
  top_857 = deg_macro1 %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = 20)
)


# ----------------------------- #
# MACROPHAGE 2 (855 vs 857)
# ----------------------------- #
macro2 <- subset(integrated, idents = "Macrophage_2")
macro2 <- PrepSCTFindMarkers(macro2)

Idents(macro2) <- macro2$sample

deg_macro2 <- FindMarkers(macro2, ident.1 = "855", ident.2 = "857",
                          logfc.threshold = 0, min.pct = 0.1)

top_macro2 <- list(
  top_855 = deg_macro2 %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20),
  top_857 = deg_macro2 %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
)


# ----------------------------- #
# MACROPHAGE 3 (855 vs 857)
# ----------------------------- #
macro3 <- subset(integrated, idents = "Macrophage_3")
macro3 <- PrepSCTFindMarkers(macro3)

Idents(macro3) <- macro3$sample

deg_macro3 <- FindMarkers(macro3, ident.1 = "855", ident.2 = "857",
                          logfc.threshold = 0, min.pct = 0.1)

top_macro3 <- list(
  top_855 = deg_macro3 %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20),
  top_857 = deg_macro3 %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
)

# Add gene names as a column for all DEG tables
deg_macro1 <- deg_macro1 %>% tibble::rownames_to_column("gene")
deg_macro2 <- deg_macro2 %>% tibble::rownames_to_column("gene")
deg_macro3 <- deg_macro3 %>% tibble::rownames_to_column("gene")

# Re-extract top lists (now gene column will be included)
top_macro1 <- list(
  top_855 = deg_macro1 %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 20),
  
  top_857 = deg_macro1 %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = 20)
)

top_macro2 <- list(
  top_855 = deg_macro2 %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20),
  top_857 = deg_macro2 %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
)

top_macro3 <- list(
  top_855 = deg_macro3 %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 20),
  top_857 = deg_macro3 %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 20)
)

# Combine all into one named list (each = one Excel sheet)
all_sheets <- c(top_macro1, top_macro2, top_macro3)

# Write to Excel (each element becomes a sheet)
write_xlsx(all_sheets, path = "Macrophage_DEGs_855_vs_857.xlsx")

##Dotplots for macrophages
library(Seurat)
library(dplyr)
library(patchwork)

# ----------------------------- #
# Function to generate dotplot with color by sample
# ----------------------------- #
plot_dot <- function(obj, cluster, genes, sample_id, title, color) {
  sub <- subset(obj, idents = cluster)
  sub <- subset(sub, subset = sample == sample_id)
  
  DotPlot(sub, features = genes) +
    RotatedAxis() +
    ggtitle(title) +
    scale_color_gradient(low = "grey90", high = color) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# ----------------------------- #
# Generate plots
# ----------------------------- #

p1 <- plot_dot(integrated, "Macrophage_1", genes_m1$up_855, "855", "Macro1 - 855", "blue")
p2 <- plot_dot(integrated, "Macrophage_1", genes_m1$up_857, "857", "Macro1 - 857", "red")

p3 <- plot_dot(integrated, "Macrophage_2", genes_m2$up_855, "855", "Macro2 - 855", "blue")
p4 <- plot_dot(integrated, "Macrophage_2", genes_m2$up_857, "857", "Macro2 - 857", "red")

p5 <- plot_dot(integrated, "Macrophage_3", genes_m3$up_855, "855", "Macro3 - 855", "blue")
p6 <- plot_dot(integrated, "Macrophage_3", genes_m3$up_857, "857", "Macro3 - 857", "red")

# ----------------------------- #
# Combine with shared legend
# ----------------------------- #

(p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

##create Dotplots for avg log2FC and pval
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)

# ----------------------------- #
# Function to prep data for each panel
# ----------------------------- #
prep_panel <- function(deg, direction, label) {
  
  if (!"gene" %in% colnames(deg)) {
    deg <- deg %>% rownames_to_column("gene")
  }
  
  if (direction == "855") {
    df <- deg %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 5)
  } else {
    df <- deg %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = 5)
  }
  
  df$panel <- label
  df$neglog10_p <- -log10(df$p_val_adj + 1e-300)
  
  return(df)
}

# ----------------------------- #
# Prepare all 6 panels
# ----------------------------- #

m1_855 <- prep_panel(deg_macro1, "855", "Macro1 - 855")
m1_857 <- prep_panel(deg_macro1, "857", "Macro1 - 857")

m2_855 <- prep_panel(deg_macro2, "855", "Macro2 - 855")
m2_857 <- prep_panel(deg_macro2, "857", "Macro2 - 857")

m3_855 <- prep_panel(deg_macro3, "855", "Macro3 - 855")
m3_857 <- prep_panel(deg_macro3, "857", "Macro3 - 857")

# ----------------------------- #
# Plot function
# ----------------------------- #
plot_panel <- function(df) {
  
  ggplot(df, aes(x = gene, y = panel)) +
    geom_point(aes(size = neglog10_p, color = avg_log2FC)) +
    
    scale_color_gradient2(
      low = "red", mid = "white", high = "blue", midpoint = 0
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    ) +
    
    labs(
      x = "Features",
      y = "Identity",
      color = "log2FC",
      size = "-log10(padj)"
    ) +
    
    ggtitle(unique(df$panel))
}

# ----------------------------- #
# Generate plots
# ----------------------------- #

p1 <- plot_panel(m1_855)
p2 <- plot_panel(m1_857)

p3 <- plot_panel(m2_855)
p4 <- plot_panel(m2_857)

p5 <- plot_panel(m3_855)
p6 <- plot_panel(m3_857)

# ----------------------------- #
# Combine into 6-panel layout
# ----------------------------- #

(p1 | p2) / (p3 | p4) / (p5 | p6)

(p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ----------------------------- #
# ----------------------------- #
##MAcrophage analysis for day 25
# ----------------------------- #
# ----------------------------- #
library(Seurat)
library(dplyr)
library(tibble)

# -------------------------------------------------- #
# Helper: get DEGs + top genes (SNAhigh vs SNAlow)
# -------------------------------------------------- #
run_deg <- function(obj, cluster_name) {
  
  # subset cluster
  sub <- subset(obj, idents = cluster_name)
  
  # IMPORTANT: set identity to condition
  Idents(sub) <- sub$SNA_level_D25
  
  # If using SCT, uncomment this:
  # sub <- PrepSCTFindMarkers(sub)
  
  # Run DEG
  deg <- FindMarkers(
    sub,
    ident.1 = "SNAhigh",
    ident.2 = "SNAlow",
    logfc.threshold = 0,
    min.pct = 0.1
  )
  
  # Add gene column safely
  if (!"gene" %in% colnames(deg)) {
    deg <- deg %>% rownames_to_column("gene")
  }
  
  # Extract top genes
  top <- list(
    top_SNAhigh = deg %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20),
    
    top_SNAlow = deg %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = 20)
  )
  
  return(list(deg = deg, top = top))
}

# -------------------------------------------------- #
# Run for all macrophage clusters
# -------------------------------------------------- #

res_macro1 <- run_deg(sna, "Macrophage_1")
res_macro2 <- run_deg(sna, "Macrophage_2")
res_macro3 <- run_deg(sna, "Macrophage_3")
res_macro4 <- run_deg(sna, "Macrophage_4")

# Access results like:
# res_macro1$deg
# res_macro1$top$top_SNAhigh
# res_macro1$top$top_SNAlow
library(writexl)

# ----------------------------- #
# Combine all top DEG tables into one list (each = one sheet)
# ----------------------------- #

all_sheets <- list(
  Macro1_SNAhigh = res_macro1$top$top_SNAhigh,
  Macro1_SNAlow  = res_macro1$top$top_SNAlow,
  
  Macro2_SNAhigh = res_macro2$top$top_SNAhigh,
  Macro2_SNAlow  = res_macro2$top$top_SNAlow,
  
  Macro3_SNAhigh = res_macro3$top$top_SNAhigh,
  Macro3_SNAlow  = res_macro3$top$top_SNAlow,
  
  Macro4_SNAhigh = res_macro4$top$top_SNAhigh,
  Macro4_SNAlow  = res_macro4$top$top_SNAlow
)

# ----------------------------- #
# Export to Excel
# ----------------------------- #

write_xlsx(all_sheets, path = "Macrophage_D25_SNAhigh_vs_SNAlow_DEGs.xlsx")

##Doptlots for day 25
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)

# ----------------------------- #
# Prep function (same logic, new labels)
# ----------------------------- #
prep_panel <- function(deg, direction, label) {
  
  if (!"gene" %in% colnames(deg)) {
    deg <- deg %>% rownames_to_column("gene")
  }
  
  if (direction == "SNAhigh") {
    df <- deg %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 5)
  } else {
    df <- deg %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = 5)
  }
  
  df$panel <- label
  df$neglog10_p <- -log10(df$p_val_adj + 1e-300)
  
  return(df)
}

# ----------------------------- #
# Prepare all 8 panels
# ----------------------------- #

m1_high <- prep_panel(res_macro1$deg, "SNAhigh", "Macro1 - SNAhigh")
m1_low  <- prep_panel(res_macro1$deg, "SNAlow",  "Macro1 - SNAlow")

m2_high <- prep_panel(res_macro2$deg, "SNAhigh", "Macro2 - SNAhigh")
m2_low  <- prep_panel(res_macro2$deg, "SNAlow",  "Macro2 - SNAlow")

m3_high <- prep_panel(res_macro3$deg, "SNAhigh", "Macro3 - SNAhigh")
m3_low  <- prep_panel(res_macro3$deg, "SNAlow",  "Macro3 - SNAlow")

m4_high <- prep_panel(res_macro4$deg, "SNAhigh", "Macro4 - SNAhigh")
m4_low  <- prep_panel(res_macro4$deg, "SNAlow",  "Macro4 - SNAlow")

# ----------------------------- #
# Plot function (unchanged)
# ----------------------------- #
plot_panel <- function(df) {
  
  ggplot(df, aes(x = gene, y = panel)) +
    geom_point(aes(size = neglog10_p, color = avg_log2FC)) +
    
    scale_color_gradient2(
      low = "red", mid = "white", high = "blue",
      midpoint = 0
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    ) +
    
    labs(
      x = "Features",
      y = "Identity",
      color = "log2FC",
      size = "-log10(padj)"
    ) +
    
    ggtitle(unique(df$panel))
}

# ----------------------------- #
# Generate plots
# ----------------------------- #

p1 <- plot_panel(m1_high)
p2 <- plot_panel(m1_low)

p3 <- plot_panel(m2_high)
p4 <- plot_panel(m2_low)

p5 <- plot_panel(m3_high)
p6 <- plot_panel(m3_low)

p7 <- plot_panel(m4_high)
p8 <- plot_panel(m4_low)

# ----------------------------- #
# Combine into 8-panel layout
# ----------------------------- #

(p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | p8) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# -------------------------------------------------- #
# Macrophage ## Day 40 analysis
# -------------------------------------------------- #
library(Seurat)
library(dplyr)
library(tibble)

# -------------------------------------------------- #
# Helper: DEG + top genes (Omentum vs Mesentery)
# -------------------------------------------------- #
run_deg_organ <- function(obj, cluster_name) {
  
  # Subset macrophage cluster
  sub <- subset(obj, idents = cluster_name)
  
  # Set identity to organ
  Idents(sub) <- sub$organ
  
  # If using SCT assay, uncomment:
  # sub <- PrepSCTFindMarkers(sub)
  
  # Run DEG (Omentum vs Mesentery)
  deg <- FindMarkers(
    sub,
    ident.1 = "omentum",
    ident.2 = "mesentery",
    logfc.threshold = 0,
    min.pct = 0.1
  )
  
  # Add gene column safely
  if (!"gene" %in% colnames(deg)) {
    deg <- deg %>% rownames_to_column("gene")
  }
  
  # Extract top genes
  top <- list(
    top_omentum = deg %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 20),
    
    top_mesentery = deg %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = 20)
  )
  
  return(list(deg = deg, top = top))
}

# -------------------------------------------------- #
# Run for macrophage clusters
# -------------------------------------------------- #

res_macro1_org <- run_deg_organ(snahigh, "Macrophage_1")
res_macro2_org <- run_deg_organ(snahigh, "Macrophage_2")

# Access:
# res_macro1_org$deg
# res_macro1_org$top$top_omentum
# res_macro1_org$top$top_mesentery

library(writexl)

# ----------------------------- #
# Combine all sheets
# ----------------------------- #

all_sheets <- list(
  Macro1_Omentum   = res_macro1_org$top$top_omentum,
  Macro1_Mesentery = res_macro1_org$top$top_mesentery,
  
  Macro2_Omentum   = res_macro2_org$top$top_omentum,
  Macro2_Mesentery = res_macro2_org$top$top_mesentery
)

# ----------------------------- #
# Export to Excel
# ----------------------------- #

write_xlsx(all_sheets, path = "Macrophage_D40_Omentum_vs_Mesentery_DEGs.xlsx")

library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)

# ----------------------------- #
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)

# ----------------------------- #
# Prepare DEG panels
# ----------------------------- #
prep_shared_genes <- function(deg, label1, label2) {
  
  if (!"gene" %in% colnames(deg)) {
    deg <- deg %>% rownames_to_column("gene")
  }
  
  # Get top 5 from BOTH directions
  genes <- c(
    deg %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 5) %>% pull(gene),
    deg %>% arrange(avg_log2FC) %>% slice_head(n = 5) %>% pull(gene)
  )
  
  df <- deg %>% filter(gene %in% genes)
  df$neglog10_p <- -log10(df$p_val_adj + 1e-300)
  
  df1 <- df; df1$panel <- label1
  df2 <- df; df2$panel <- label2
  
  return(list(df1 = df1, df2 = df2))
}
# ----------------------------- #
# Create 4 datasets
# ----------------------------- #

m1 <- prep_shared_genes(res_macro1_org$deg, "Macro1 - Omentum", "Macro1 - Mesentery")
m2 <- prep_shared_genes(res_macro2_org$deg, "Macro2 - Omentum", "Macro2 - Mesentery")

m1_omen <- m1$df1
m1_mesen <- m1$df2

m2_omen <- m2$df1
m2_mesen <- m2$df2
# ----------------------------- #
# Plot function (blue-only)
# ----------------------------- #
plot_panel <- function(df) {
  
  ggplot(df, aes(x = gene, y = panel)) +
    geom_point(aes(size = neglog10_p, color = abs(avg_log2FC))) +
    
    scale_color_gradient(
      low = "lightblue",
      high = "blue"
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    ) +
    
    labs(
      x = "Features",
      y = "Identity",
      color = "|log2FC|",
      size = "-log10(padj)"
    ) +
    
    ggtitle(unique(df$panel))
}

# ----------------------------- #
# Generate plots
# ----------------------------- #

p1 <- plot_panel(m1_omen)
p2 <- plot_panel(m1_mesen)

p3 <- plot_panel(m2_omen)
p4 <- plot_panel(m2_mesen)

# ----------------------------- #
# Combine (4-panel layout)
# ----------------------------- #

(p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

