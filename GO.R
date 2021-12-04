## SCRIPT: Go analysis

## 04.12.21 Laura Sudupe , git @lsudupe


#Libraries-------------------------------------------------
source("packages.R")
library("clusterProfiler")
library("org.Hs.eg.db")
library("AnnotationHub")

###Read data
integrated <- readRDS("./objects/processed/integrated.label.rds")
de <- readRDS("./results/integrated_label_markers.rds")

df <- de[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)


dfsample$`CD` = bitr(dfsample$`CD`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`Immune cells` = bitr(dfsample$`Immune cells`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`EC` = bitr(dfsample$`EC`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`FB` = bitr(dfsample$`FB`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#do the same here, a line like below for each cluster
genelist <- list("CD" = dfsample$`CD`$ENTREZID, 
                 "Immune cells" = dfsample$`Immune cells`$ENTREZID,
                 "EC" = dfsample$`EC`$ENTREZID,
                 "FB" = dfsample$`FB`$ENTREZID)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")

y <- c("CD","Immune cells","EC","FB")
##Plot
pdf(file.path("./results/GO",filename = "dotplot.pdf"))
dotplot(GOclusterplot)
dev.off()




