##SNAhigh in Omentum at 3 different time points

##Start after line 95 to avoid re- extractign, and re pre processing! Just import and do downstream analysis
library(Seurat)
##Extarcting SNAhigh from different time points
### -------------------------
### Day 14 – SNAhigh (Sample 855)
### -------------------------

day14_snahigh <- subset(
  integrated,
  subset = sample == "855"
)

### -------------------------
### Day 25 – SNAhigh
### -------------------------

day25_snahigh <- subset(
  sna,
  subset = SNA_level_D25 == "SNAhigh"
)

### -------------------------
### Day 40 – SNAhigh Omentum
### -------------------------

day40_snahigh_omentum <- subset(
  snahigh,
  subset = organ == "omentum"
)

### -------------------------
### Quick sanity checks
### -------------------------

table(day14_snahigh$sample)
table(day25_snahigh$SNA_level_D25)
table(day40_snahigh_omentum$organ)

### Cell counts
ncol(day14_snahigh)
ncol(day25_snahigh)
ncol(day40_snahigh_omentum)

### Add timepoint metadata
day14_snahigh$timepoint <- "D14"
day25_snahigh$timepoint <- "D25"
day40_snahigh_omentum$timepoint <- "D40"

### Make cell names unique
day14_snahigh <- RenameCells(day14_snahigh, add.cell.id = "D14")
day25_snahigh <- RenameCells(day25_snahigh, add.cell.id = "D25")
day40_snahigh_omentum <- RenameCells(day40_snahigh_omentum, add.cell.id = "D40")

##Merge the SNAhigh together
sna_timecourse <- merge(
  x = day14_snahigh,
  y = list(day25_snahigh, day40_snahigh_omentum)
)

table(sna_timecourse$timepoint)

DefaultAssay(sna_timecourse) <- "RNA"

sna_timecourse <- NormalizeData(sna_timecourse)
sna_timecourse <- FindVariableFeatures(sna_timecourse)
sna_timecourse <- ScaleData(sna_timecourse)
sna_timecourse <- RunPCA(sna_timecourse)
sna_timecourse <- RunUMAP(sna_timecourse, dims = 1:20)

DimPlot(sna_timecourse, group.by = "timepoint")

table(sna_timecourse$seurat_clusters, sna_timecourse$timepoint)
library(ggplot2)

prop <- prop.table(
  table(sna_timecourse$seurat_clusters, sna_timecourse$timepoint),
  margin = 2
)

prop_df <- as.data.frame(prop)

colnames(prop_df) <- c("cluster", "timepoint", "proportion")

ggplot(prop_df, aes(x = timepoint, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Cell proportion") +
  xlab("Timepoint") +
  labs(fill = "Cluster")

saveRDS(sna_timecourse, file = "sna_timecourse_seurat.rds")
