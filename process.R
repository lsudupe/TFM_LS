## SCRIPT: Process our spatial object

## 28.11.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

###Read data

combined <- readRDS("./objects/processed/sp.combined.rds")


###Transformation

sc.combined.sct <- SCTransform(combined, assay = "Spatial",verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)