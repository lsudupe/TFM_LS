## SCRIPT: Cluster analysis

## 28.11.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

###Read data

integrated <- readRDS("./objects/processed/integrated.sct.rds")
saveRDS(integrated, "./objects/processed/integrated.sct.rds")
combined <- readRDS("./objects/processed/combined.sct.rds")

###Plots

pdf(file.path("./results/clusters",filename = "umap_r0.4_combined_ident.pdf"))
DimPlot(combined, group.by = c("sample"), label = T) + ggtitle("UMAP_r0.4")
dev.off()

pdf(file.path("./results/clusters",filename = "spatial_integrated.pdf"))
SpatialDimPlot(object = integrated ,group.by = c("ident"), pt.size.factor = 80)
dev.off()


####Some FB genes

Aspn <- SpatialFeaturePlot(integrated, 
                            features = c("Aspn"),
                            alpha = 0.6,
                            pt.size.factor = 80, 
                            combine = FALSE)
fix.sc.Aspn<- scale_fill_continuous(limits = c(0,10), breaks = c(0,10), type ="viridis")
Aspn. <- lapply(Aspn, function (x) x + fix.sc.Aspn)

pdf(file.path("./results/genes",filename = "Aspn.pdf"))
CombinePlots(Aspn.)
dev.off()

pdf(file.path("./results/genes",filename = "prueba.pdf"))
SpatialFeaturePlot(integrated, 
                   features = c("Aspn"),
                   alpha = 0.6,
                   pt.size.factor = 80, 
                   combine = FALSE)
dev.off()
