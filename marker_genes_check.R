## SCRIPT: marker genes check

## 28.11.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")
source("marker_genes.R")

###Read data

a<- readRDS("./objects/processed/integrated.sct.rds")

###################################################3

###Plot and save heatmap

###################################################3
png(file.path("./results/genes/marker_genes_check",filename = "Macrophagues feature.png"))
DoHeatmap(subset(a, downsample = 100), features = Macrophagues, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename = "Schwann.cells feature.png"))
DoHeatmap(subset(a, downsample = 100), features = Schwann.cells, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="DC feature.png"))
DoHeatmap(subset(a, downsample = 100), features = DC.like.cells, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Fibro feature.png"))
DoHeatmap(subset(a, downsample = 100), features = fibroblast, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Circulation fibro feature.png"))
DoHeatmap(subset(a, downsample = 100), features = circulation.fibroblast, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Cardiomyocites feature.png"))
DoHeatmap(subset(a, downsample = 100), features = cardiomyocite, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Endohtelial feature.png"))
DoHeatmap(subset(a, downsample = 100), features = endothelial, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Smooth muscle feature.png"))
DoHeatmap(subset(a, downsample = 100), features = smooth.muscle, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="B cell feature.png"))
DoHeatmap(subset(a, downsample = 100), features = B.cell, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Machofage/monocytes feature.png"))
DoHeatmap(subset(a, downsample = 100), features = machofage.monocytes, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="T cells feature.png"))
DoHeatmap(subset(a, downsample = 100), features = T.cell, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Neutrophils feature.png"))
DoHeatmap(subset(a, downsample = 100), features = neutrophils, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Dentritic feature.png"))
DoHeatmap(subset(a, downsample = 100), features = dentritic, size = 3)
dev.off()
png(file.path("./results/genes/marker_genes_check",filename ="Granulocytes feature.png"))
DoHeatmap(subset(a, downsample = 100), features = granulocytes, size = 3)
dev.off()

###################################################3

###Plot and save doheatmap

###################################################3

pdf(file.path("./results/genes/marker_genes_check",filename ="Dentritic feature umap.pdf"))
FeaturePlot(a,features = dentritic)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename = "Macrophagues feature umap.pdf"))
FeaturePlot(a, features = Macrophagues)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename = "Schwann.cells feature umap.pdf"))
FeaturePlot(a, features = Schwann.cells)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="DC feature umap.pdf"))
FeaturePlot(a,features = DC.like.cells)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Fibro feature umap.pdf"))
FeaturePlot(a, features = fibroblast)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Circulation fibro feature umap.pdf"))
FeaturePlot(a, features = circulation.fibroblast)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Cardiomyocites feature umap.pdf"))
FeaturePlot(a, features = cardiomyocite)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Endohtelial feature umap.pdf"))
FeaturePlot(a,features = endothelial)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Smooth muscle feature umap.pdf"))
FeaturePlot(a,features = smooth.muscle)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="B cell feature umap.pdf"))
FeaturePlot(a,features = B.cell)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Machofage/monocytes feature umap.pdf"))
FeaturePlot(a,features = machofage.monocytes)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="T cells feature umap.pdf"))
FeaturePlot(a,features = T.cell)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Neutrophils feature umap.pdf"))
FeaturePlot(a, features = neutrophils)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Dentritic feature umap.pdf"))
FeaturePlot(a,features = dentritic)
dev.off()
pdf(file.path("./results/genes/marker_genes_check",filename ="Granulocytes feature umap.pdf"))
FeaturePlot(a, features = granulocytes)
dev.off()

###################################################3

###Plot and save spatialplot

###################################################3

FB <- SpatialFeaturePlot(a, 
                           features = fibroblast,
                           alpha = 0.6,
                           pt.size.factor = 80, 
                           combine = FALSE)
fix.sc.FB<- scale_fill_continuous(limits = c(0,10), breaks = c(0,10), type ="viridis")
FB. <- lapply(FB, function (x) x + fix.sc.FB)

pdf(file.path("./results/genes/marker_genes_check/spatial",filename = "FB.pdf"))
CombinePlots(FB.)
dev.off()


