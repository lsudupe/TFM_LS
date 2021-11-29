## SCRIPT: Differential expression analysis

## 29.11.21 Laura Sudupe , git @lsudupe
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html

#Libraries-------------------------------------------------
source("packages.R")

###Read data
integrated <- readRDS("./objects/processed/integrated.sct.rds")
de <- readRDS("./results/integrated_markers.rds")


top5 <- de %>% group_by(cluster) %>% top_n(5, avg_log2FC) 
DoHeatmap(
  object = a,
  features = top5$gene
)

#####Pathway
library(ReactomeGSA)
gsva_result <- analyse_sc_clusters(integrated, verbose = TRUE, assay= "integrated",slot = "data")

pathway_expression <- pathways(gsva_result)

# simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

# find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))


max_difference$diff <- max_difference$max - max_difference$min


# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]


head(max_difference)


##PLOTING THE RESULTS
#1
pdf(file.path("./results/pathway",filename = "1.pdf"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
dev.off()


#2
# Additional parameters are directly passed to gplots heatmap.2 function
pdf(file.path("./results/pathway",filename = "2.pdf"))
plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(6,20))
dev.off()

## limit to selected B cell related pathways
relevant_pathways <- c("R-HSA-2393930", "R-HSA-141333", "R-HSA-8964041")
pdf(file.path("./results/pathway",filename = "3.pdf"))
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,30), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left
dev.off()

pdf(file.path("./results/pathway",filename = "4.pdf"))
plot_gsva_pca(gsva_result)
dev.off()
