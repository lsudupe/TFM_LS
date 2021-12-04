## SCRIPT: Label the clusters

## 04.12.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

###Read data

integrated <- readRDS("./objects/processed/integrated.sct.rds")


###Label the clusters

integrated$ident <- integrated@meta.data[["seurat_clusters"]]
Seurat::Idents(object = integrated) <- integrated@meta.data[["ident"]]

integrated.label <- RenameIdents(object = integrated, `0` = "CD",
                                    `1` = "CD",`2` = "Immune cells",`3` = "EC",
                                    `4` = "FB",`5` = "CD",`6` = "FB")
saveRDS(integrated.label, "./objects/processed/integrated.label.rds")


###Plot 

integrated.label$ident <- (integrated.label@active.ident)

pdf(file.path("./results/clusters",filename = "umap_r0.4_integrated_label.pdf"))
DimPlot(integrated.label, group.by = c("ident"), label = T) + ggtitle("UMAP_r0.4")
dev.off()

pdf(file.path("./results/clusters",filename = "spatial_integrated_label.pdf"))
SpatialDimPlot(object = integrated.label ,group.by = c("ident"), pt.size.factor = 80)
dev.off()








