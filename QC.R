## SCRIPT: QC of the seurat spatial objects

## 30.10.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

#Read objects
sham <- readRDS("./objects/initial/sham.sp.rds")
day1.1 <- readRDS("./objects/initial/day1.1.sp.rds")
day1.2 <- readRDS("./objects/initial/day1.2.sp.rds")
day7.1 <- readRDS("./objects/initial/day7.2.sp.rds")
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
day7.1@meta.data$sample <- "day7.2"

