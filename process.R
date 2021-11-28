## SCRIPT: Process our spatial object

## 28.11.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

###Read data

combined <- readRDS("./objects/processed/sp.combined.rds")

###Integration
###https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
list <- SplitObject(combined, split.by = "sample")
list <- lapply(X = list, FUN = SCTransform, assay="SCT")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)


anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                              anchor.features = features)
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

saveRDS(combined, "./objects/processed/integrated.rds")
integrated <- readRDS("./objects/processed/integrated.rds")


###Transformation

integrated <- SCTransform(integrated, assay = "Spatial",verbose = FALSE) %>%
  RunPCA(assay = "SCT",npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)%>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

saveRDS(integrated, "./objects/processed/integrated.sct.rds")


###Markers

markers <- Seurat::FindAllMarkers(object = integrated, 
                                          assay = "SCT",
                                          slot = "data",
                                          verbose = TRUE, 
                                          only.pos = TRUE)

saveRDS(markers, "./results/integrated_markers.rds")





