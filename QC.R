## SCRIPT: QC of the seurat spatial objects

## 30.10.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

#Read objects
sham <- readRDS("./objects/initial/sham.sp.rds")
day1.1 <- readRDS("./objects/initial/day1.1.sp.rds")
day1.2 <- readRDS("./objects/initial/day1.2.sp.rds")
day7.1 <- readRDS("./objects/initial/day7.1.sp.rds")
day7.2 <- readRDS("./objects/initial/day7.2.sp.rds")

#Add number of genes per UMI for each cell to metadata
sham$log10GenesPerUMI <- log10(sham$nFeature_Spatial) / log10(sham$nCount_Spatial)
day1.1$log10GenesPerUMI <- log10(day1.1$nFeature_Spatial) / log10(day1.1$nCount_Spatial)
day1.2$log10GenesPerUMI <- log10(day1.2$nFeature_Spatial) / log10(day1.2$nCount_Spatial)
day7.1$log10GenesPerUMI <- log10(day7.1$nFeature_Spatial) / log10(day7.1$nCount_Spatial)
day7.2$log10GenesPerUMI <- log10(day7.2$nFeature_Spatial) / log10(day7.2$nCount_Spatial)

# Compute percent mito ratio
sham <- PercentageFeatureSet(sham, "^mt-", col.name = "percent_mito")
day1.1 <- PercentageFeatureSet(day1.1, "^mt-", col.name = "percent_mito")
day1.2 <- PercentageFeatureSet(day1.2, "^mt-", col.name = "percent_mito")
day7.1 <- PercentageFeatureSet(day7.1, "^mt-", col.name = "percent_mito")
day7.2 <- PercentageFeatureSet(day7.2, "^mt-", col.name = "percent_mito")

##add sample column
sham@meta.data$sample <- "sham"
day1.1@meta.data$sample <- "day1.1"
day1.2@meta.data$sample <- "day1.2"
day7.1@meta.data$sample <- "day7.1"
day7.2@meta.data$sample <- "day7.2"

##merge them
sp.combined <- merge(sham, y = c(day1.1, day1.2, day7.1,day7.2 ), 
                     add.cell.ids = c("sham","day1.1","day1.2", "day7.1","day7.2"), project = "Tokio")

##visualization
# Visualize the number of spots counts per sample
png(file.path("./results/QC",filename = "Number of spot count per sample.png"))
sp.combined@meta.data%>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

# Visualize the distribution of genes detected per spot via histogram
png(file.path("./results/QC",filename = "genes detected per spot histogram.png"))
sp.combined@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_Spatial, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
dev.off()
# Visualize the distribution of genes detected per spot via boxplot
png(file.path("./results/QC",filename = "genes detected per spot boxplot.png"))
sp.combined@meta.data %>% 
  ggplot(aes(x=sample, y=log10(nFeature_Spatial), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Ngenes vs Npots")
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per spot
png(file.path("./results/QC",filename = "mito percentage per spot.png"))
sp.combined@meta.data %>% 
  ggplot(aes(color=sample, x=percent_mito, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
dev.off()
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
sp.combined@meta.data %>% 
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
#QUE ES ESTO

png(file.path("./results/QC",filename = "nCount violin.png"))
VlnPlot(sp.combined, features = c("nCount_Spatial"),group.by = "sample", pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()
png(file.path("./results/QC",filename = "nFeature violin.png"))
VlnPlot(sp.combined, features = c("nFeature_Spatial"),group.by = "sample",pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()
png(file.path("./results/QC",filename = "mito violin.png"))
VlnPlot(sp.combined, features = c("percent_mito"),group.by = "sample", pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

##change sample column
sham@meta.data$sample <- "sham"
day1.1@meta.data$sample <- "day1"
day7.1@meta.data$sample <- "day7"

##merge them
sp.combined <- merge(sham, y = c(day1.1, day7.1), 
                     add.cell.ids = c("sham","day1", "day7"), project = "Tokio")
##brief filter and save the selected
saveRDS(sham, file = "./objects/processed/sham.rds")
saveRDS(day1.1, file = "./objects/processed/day1.rds")
saveRDS(day7.1, file = "./objects/processed/day7.rds")
saveRDS(sp.combined, file = "./objects/processed/sp.combined.rds")

combined <- readRDS("./objects/processed/sp.combined.rds")

##############SPATIAL QC

##mito percentage
mito <- SpatialFeaturePlot(combined, 
                         features = c("percent_mito"),
                         alpha = 0.6,
                         pt.size.factor = 80, 
                         combine = FALSE)
fix.sc.mito <- scale_fill_continuous(limits = c(0.0,75.0), breaks = c(0.0, 35.0, 75.0), type ="viridis")
mito. <- lapply(mito, function (x) x + fix.sc.mito)

pdf(file.path("./results/QC",filename = "mito spatial.pdf"))
CombinePlots(mito.)
dev.off()

##spatial features
feature <- SpatialFeaturePlot(combined, 
                           features = c("nFeature_Spatial"),
                           alpha = 0.6,
                           pt.size.factor = 80, 
                           combine = FALSE)
fix.sc.feature <- scale_fill_continuous(limits = c(0,9500), breaks = c(0,5000, 9500), type ="viridis")
feature. <- lapply(feature, function (x) x + fix.sc.feature)

pdf(file.path("./results/QC",filename = "nfeature_combined.pdf"))
CombinePlots(feature.)
dev.off()

##spatial count
count <- SpatialFeaturePlot(combined, 
                              features = c("nCount_Spatial"),
                              alpha = 0.6,
                              pt.size.factor = 80, 
                              combine = FALSE)
fix.sc.count <- scale_fill_continuous(limits = c(0,90000), breaks = c(0,50000, 90000), type ="viridis")
count. <- lapply(count, function (x) x + fix.sc.count)

pdf(file.path("./results/QC",filename = "nCount_combined.pdf"))
CombinePlots(count.)
dev.off()







