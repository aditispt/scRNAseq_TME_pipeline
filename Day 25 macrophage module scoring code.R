# ================================================== #
# MACROPHAGE MODULE SCORING + COMPARISON + TRAJECTORY
# Uses existing objects: integrated, sna, snahigh
# ================================================== #

library(Seurat)
library(dplyr)
library(ggplot2)

# -------------------------------------------------- #
# 1. DEFINE GENE MODULES (based on YOUR DEG results)
# -------------------------------------------------- #

IFN_genes <- c("Ifitm1", "Rsad2", "Ifitm6", "Slfn4")
ECM_genes <- c("Ctrb1", "Prss2", "Reg3b", "Reg2", "Mgp", "Pcolce")
M2_genes <- c("Retnla", "C3", "Col3a1", "Prss2")
Metabolic_genes <- c("Tars", "Slc25a44", "Gstt2", "Acod9")

### IFN / ANTIVIRAL
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
# 2. ADD MODULE SCORES (DAY 25 object: sna)
# -------------------------------------------------- #

sna <- AddModuleScore(
  object = sna,
  features = list(
    IFN_genes,
    ECM_genes,
    M2_genes,
    Metabolic_genes
  ),
  name = c("IFN", "ECM", "M2", "Metabolic")
)

# Clean naming (important)
sna$IFN_score <- sna$IFN1
sna$ECM_score <- sna$ECM2
sna$M2_score <- sna$M23
sna$Metabolic_score <- sna$Metabolic4

# -------------------------------------------------- #
# 3. SUBSET MACROPHAGES (same naming as your code)
# -------------------------------------------------- #

macro_d25 <- subset(
  sna,
  idents = c("Macrophage_1", "Macrophage_2", "Macrophage_3", "Macrophage_4")
)

# -------------------------------------------------- #
# 4. SNAhigh vs SNAlow COMPARISON (Day 25)
# -------------------------------------------------- #

# Violin plots
VlnPlot(
  macro_d25,
  features = c("IFN_score", "ECM_score", "M2_score", "Metabolic_score"),
  group.by = "SNA_level_D25",
  pt.size = 0
)

# -------------------------------------------------- #
# 5. STATISTICS (Wilcoxon test)
# -------------------------------------------------- #

meta_d25 <- macro_d25@meta.data

wilcox_IFN <- wilcox.test(IFN_score ~ SNA_level_D25, data = meta_d25)
wilcox_ECM <- wilcox.test(ECM_score ~ SNA_level_D25, data = meta_d25)
wilcox_M2  <- wilcox.test(M2_score ~ SNA_level_D25, data = meta_d25)
wilcox_MET <- wilcox.test(Metabolic_score ~ SNA_level_D25, data = meta_d25)

# Print results
wilcox_IFN
wilcox_ECM
wilcox_M2
wilcox_MET

