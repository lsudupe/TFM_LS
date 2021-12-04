## SCRIPT: Differential expressed genes

## 04.12.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

###Read data

de <- readRDS("./results/integrated_markers.rds")
integrate <- readRDS("./objects/processed/integrated.label.rds")

### Extract top 5 markers per cluster
top5 <- de %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_log2FC)

###Plots

pdf(file.path("./results/genes",filename = "violinplot_de_genes.pdf"))
VlnPlot(object = integrate, 
        features = c("Igfbp7", "Tpm2", "Sncg", "Mustn1", "Ppp1r14a"))
dev.off()


markers <- Seurat::FindAllMarkers(object = integrate, 
                                  assay = "integrated",
                                  slot = "data",
                                  logfc.threshold = 0.25,
                                  verbose = TRUE, 
                                  only.pos = TRUE)

saveRDS(markers, "./results/integrated_label_markers.rds")

### Extract top 5 markers per cluster

top5.label <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5,
        wt = avg_log2FC)

pdf(file.path("./results/genes",filename = "violinplot_de_genes.pdf"))
VlnPlot(object = integrate, 
        features = c("Igfbp7", "Tpm2", "Sncg", "Mustn1", "Ppp1r14a"))
dev.off()

# Plot interesting marker gene expression for cluster FB
pdf(file.path("./results/genes",filename = "umapplot_de_genes.pdf"))
FeaturePlot(object = integrate, 
            features = c("Igfbp7", "Tpm2", "Sncg", "Mustn1", "Ppp1r14a"),
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()

##plot spatially
pdf(file.path("./results/genes",filename = "spatial_Igfbp7.pdf"))
SpatialFeaturePlot(object = integrate, 
            features = c("Igfbp7"),
            alpha = 0.6,
            pt.size.factor = 80, 
            combine = FALSE)
dev.off()

pdf(file.path("./results/genes",filename = "spatial_Tpm2.pdf"))
SpatialFeaturePlot(object = integrate, 
                   features = c("Tpm2"),
                   alpha = 0.6,
                   pt.size.factor = 80, 
                   combine = FALSE)
dev.off()

pdf(file.path("./results/genes",filename = "spatial_Sncg.pdf"))
SpatialFeaturePlot(object = integrate, 
                   features = c("Sncg"),
                   alpha = 0.6,
                   pt.size.factor = 80, 
                   combine = FALSE)
dev.off()

pdf(file.path("./results/genes",filename = "spatial_Mustn1.pdf"))
SpatialFeaturePlot(object = integrate, 
                   features = c("Mustn1"),
                   alpha = 0.6,
                   pt.size.factor = 80, 
                   combine = FALSE)
dev.off()
