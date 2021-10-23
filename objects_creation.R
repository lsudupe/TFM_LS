## SCRIPT: Spatial seurat object creation with visium data

## 23.10.21 Laura Sudupe , git @lsudupe

#Libraries-------------------------------------------------
source("packages.R")

#objects directory
dir.create("./objects/initial")

#Data-----------------------------------------------------
sham <- Read10X(
  data.dir = "./data/visium/sham/",
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
) 
day1.1 <- Read10X(
  data.dir = "./data/visium/day1/sample1",
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
) 
day1.2 <- Read10X(
  data.dir = "./data/visium/day1/sample2",
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
) 
day7.1 <- Read10X(
  data.dir = "./data/visium/day7/sample1",
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE
) 
day7.2 <- Read10X(
  data.dir = "./data/visium/day7/sample2",
  gene.column = 2,
  unique.features = TRUE,
  strip.suffix = FALSE)

###########################################################
#create seurat object 
sp.sham <- CreateSeuratObject(
  sham,
  project = "Tokio",
  assay = "Spatial",
  min.cells = 3,
  min.features = 200
)
sp.day1.1 <- CreateSeuratObject(
  day1.1,
  project = "Tokio",
  assay = "Spatial",
  min.cells = 3,
  min.features = 200
)
sp.day1.2 <- CreateSeuratObject(
  day1.2,
  project = "Tokio",
  assay = "Spatial",
  min.cells = 3,
  min.features = 200
)
sp.day7.1 <- CreateSeuratObject(
  day7.1,
  project = "Tokio",
  assay = "Spatial",
  min.cells = 3,
  min.features = 200
)
sp.day7.2 <- CreateSeuratObject(
  day7.2,
  project = "Tokio",
  assay = "Spatial",
  min.cells = 3,
  min.features = 200
)
###########################################################
#add image
img.sham <- Read10X_Image(image.dir ="./data/visium/sham/spatial", filter.matrix = FALSE)
img.sham <- img.sham[Cells(x = sp.sham)]
DefaultAssay(spatial.sham.object = img.sham) <- "Spatial"
sp.sham[['sham.slice']] <- img.sham

img.day1.1 <- Read10X_Image(image.dir ="./data/visium/day1/sample1/spatial", filter.matrix = FALSE)
img.day1.1 <- img.day1.1[Cells(x = sp.day1.1)]
DefaultAssay(sp.day1.1 = img.day1.1) <- "Spatial"
sp.day1.1[['day1.1.slice']] <- img.day1.1

img.day1.2 <- Read10X_Image(image.dir ="./data/visium/day1/sample2/spatial", filter.matrix = FALSE)
img.day1.2 <- img.day1.1[Cells(x = sp.day1.2)]
DefaultAssay(sp.day1.2 = img.day1.2) <- "Spatial"
sp.day1.2[['day1.2.slice']] <- img.day1.2

img.day7.1 <- Read10X_Image(image.dir ="./data/visium/day7/sample1/spatial", filter.matrix = FALSE)
img.day7.1 <- img.day7.1[Cells(x = sp.day7.1)]
DefaultAssay(sp.day7.1 = img.day7.1) <- "Spatial"
sp.day7.1[['day7.1.slice']] <- img.day7.1

img.day7.2 <- Read10X_Image(image.dir ="./data/visium/day7/sample2/spatial", filter.matrix = FALSE)
img.day7.2 <- img.day7.2[Cells(x = sp.day7.2)]
DefaultAssay(sp.day7.2 = img.day7.2) <- "Spatial"
sp.day7.2[['day7.2.slice']] <- img.day7.2
###########################################################
#save objects
saveRDS(sp.sham,"./objects/initial/sham.sp.rds")
saveRDS(sp.day1.1,"./objects/initial/day1.1.sp.rds")
saveRDS(sp.day1.2,"./objects/initial/day1.2.sp.rds")
saveRDS(sp.day7.1,"./objects/initial/day7.2.sp.rds")
saveRDS(sp.day7.2,"./objects/initial/day7.2.sp.rds")






