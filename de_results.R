## SCRIPT: Differential analysis results

## 08.12.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

#Read objects
de <- readRDS("./results/integrated_label_markers.rds")

#Filter
de_subset<- subset(de, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
top3.de <- de_subset %>%
  group_by(cluster) %>%
  top_n(n = 3,
        wt = avg_log2FC)

#Save
write.csv(top3.de,"./results/genes/top3de.csv", row.names = FALSE)


