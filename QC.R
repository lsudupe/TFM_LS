## SCRIPT: QC of the seurat spatial objects

## 30.10.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

#read objects
sham <- readRDS("./objects/initial/sham.sp.rds")
day1.1 <- readRDS("./objects/initial/day1.1.sp.rds")
day1.2 <- readRDS("./objects/initial/day1.2.sp.rds")
day7.1 <- readRDS("./objects/initial/day7.2.sp.rds")
day7.2 <- readRDS("./objects/initial/day7.2.sp.rds")

# Create sample column
sham@meta.data$sample <- NA
metadata$sample[which(str_detect(metadata$protocol, "sham"))] <- "sham"
metadata$sample[which(str_detect(metadata$protocol, "day1"))] <- "day 1"
metadata$sample[which(str_detect(metadata$protocol, "day7"))] <- "day 7"

sp.list <- c(sham,day1.1,day1.2,day7.1,day7.2)

# Add number of genes per UMI for each spot to metadata
for (i in 1:length(sp.list)) {
  sp.list[[i]]$log10GenesPerUMI <- log10(sp.list[[i]]$nFeature_Spatial) / log10(sp.list[[i]]$nCount_Spatial)
}
# Create sample column
for (i in 1:length(sp.list)) {
  sp.list[[i]]@meta.data$sample <- NA
}
sp.list[[1]]@meta.data$sample <- "sham"
sp.list[[2]]@meta.data$sample <- "day1.1"
sp.list[[3]]@meta.data$sample <- "day1.2"
sp.list[[4]]@meta.data$sample <- "day7.1"
sp.list[[5]]@meta.data$sample <- "day7.1"
