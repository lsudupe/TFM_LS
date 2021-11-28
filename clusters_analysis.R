## SCRIPT: Cluster analysis

## 28.11.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

###Read data

integrated <- readRDS("./objects/processed/integrated.sct.rds")
saveRDS(integrated, "./objects/processed/integrated.sct.rds")
combined <- readRDS("./objects/processed/combined.sct.rds")

###Plots

pdf(file.path("./results/clusters",filename = "umap_r0.4_combined.pdf"))
DimPlot(combined, group.by = c("seurat_clusters"), label = T) + ggtitle("UMAP_r0.4")
dev.off()

pdf(file.path("./results/clusters",filename = "spatial_integrated.pdf"))
SpatialDimPlot(object = integrated ,group.by = c("ident"), pt.size.factor = 80)
dev.off()



