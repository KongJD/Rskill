## Reference Mapping

```properties
### 背景
单细胞数据分析的核心步骤到底是在分析基因表达矩阵协方差的来源
对细胞的cell embedding，图谱提供好的model，可看自己数据在其中的关系
```

#### 1.数据清洗

```R
library(Seurat)
library(GEOquery)
setwd(here::here())

# cached.object <- "input/SadeFeldman.seurat.rds"
cached.object <- "input/SadeFeldman.seurat.qs"

options(timeout = max(3600, getOption("timeout")))

geo_acc <- "GSE120575"
datadir <- "input/SadeFeldman"

series <- paste0(geo_acc, "_series_matrix.txt.gz")

system(paste0("mkdir -p ", datadir))
getGEOSuppFiles(geo_acc, baseDir = datadir)

## Load expression matrix and metadata
exp.mat <- read.delim(sprintf("%s/%s/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
                              datadir, geo_acc), header = F, sep = "\t")

genes <- exp.mat[c(-1, -2), 1]
cells <- as.vector(t(exp.mat[1, 2:16292]))
samples <- as.factor(t(exp.mat[2, 2:16292]))

exp.mat <- exp.mat[c(-1, -2), 2:16292]
colnames(exp.mat) <- cells
rownames(exp.mat) <- genes
exp.mat <- as.matrix(exp.mat)
class(exp.mat) <- "numeric"

meta <- read.delim(sprintf("%s/%s/GSE120575_patient_ID_single_cells.txt.gz",
                           datadir, geo_acc), header = T, sep = "\t", skip = 19, nrows = 16291)
meta <- meta[, 1:7]

treat <- factor(ifelse(grepl("Post", samples), "Post", "Pre"))
response <- factor(meta$characteristics..response)
therapy <- factor(meta$characteristics..therapy)

## Create Seurat object and add meta data
query.object <- CreateSeuratObject(counts = exp.mat, project = "SadeFeldman", min.cells = 10)
query.object@meta.data$Sample <- samples
query.object@meta.data$Time <- treat
query.object@meta.data$Response <- response
query.object@meta.data$Therapy <- therapy

# saveRDS(query.object, file = cached.object)
qs::qsave(query.object, file = cached.object)
```

```R
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
setwd(here::here())

#' Data download from: https://zenodo.org/record/5461803

#### 泛癌肿瘤浸润CD4-T细胞图谱 ####

color.set <- readRDS("input/ZhengLiangtao/panC.colSet.list.rds")
cd4.ref.sce <- readRDS("input/ZhengLiangtao/int.CD4.S35.sce.merged.rds")

cd4.ref.df <- reducedDim(cd4.ref.sce, "harmony.umap") %>% as.data.frame()
colnames(cd4.ref.df) <- paste0("UMAP_", 1:2)
cd4.ref.df$meta.cluster <- cd4.ref.sce$meta.cluster
cd4.ref.df$cluster.name <- plyr::mapvalues(cd4.ref.df$meta.cluster,
                                           from = names(color.set$meta.cluster)[1:41],
                                           to = names(color.set$cluster.name))
## a list of dgCMatrix object
expr.list <- readRDS("input/ZhengLiangtao/CD4.expr.list.rds")
sapply(expr.list, class) # expr.list[[26]] is a dgTMatrix
expr.list[[26]] <- as(expr.list[[26]], "dgCMatrix")

sapply(expr.list, dim) # ~5000 genes kept in each dataset. To be noted, the genes are not same between datasets.

## 创建Seurat对象
seu.list <- lapply(expr.list, CreateSeuratObject)
seu.ref.cd4 <- base::Reduce(merge, seu.list)

## 写入cellmeta信息
cellmeta <- readRDS("input/ZhengLiangtao/CD4.cellmeta.rds")
colnames(cellmeta)[22:23] <- paste0("UMAP_", 1:2)
cellmeta <- cellmeta[rownames(seu.ref.cd4@meta.data), ]

umap.emb <- cellmeta[, paste0("UMAP_", 1:2)] %>% as.matrix()
cellmeta <- cellmeta[, -c(22,23)]

seu.ref.cd4@meta.data <- cellmeta
seu.ref.cd4[["umap"]] <- CreateDimReducObject(umap.emb, assay = "RNA")

seu.ref.cd4$meta.cluster.name <- plyr::mapvalues(
  x = seu.ref.cd4$meta.cluster,
  from = names(color.set$meta.cluster[1:41]),
  to = names(color.set$cluster.name))

## 写入metacell的信息
text.pos.df <- data.frame(
  x = c(-4.5, 4, 3.5,  2, -4.5, -7, -6),
  y = c(  5,  4, -4,  -6, -6,   -2,  2),
  label = c("Tn", "Treg", "Tfh", "Tfh/Th1", "Temra", "Tem", "Tm")
)

seu.ref.cd4@misc$metacell <- list(
  meta.data = cd4.ref.df,
  colors = color.set$cluster.name[levels(cd4.ref.df$cluster.name)],
  text.pos = text.pos.df
)

## 写入marker基因, 基因按照显著性排序
top.genes <- readRDS("input/ZhengLiangtao/CD4.markers.rds")
seu.ref.cd4@misc$markers <- top.genes

qs::qsave(seu.ref.cd4, "input/ZhengLiangtao.CD4.seurat.qs")

#### 泛癌肿瘤浸润CD8-T细胞图谱 ####
color.set <- readRDS("input/ZhengLiangtao/panC.colSet.list.rds")
cd8.ref.sce <- readRDS("input/ZhengLiangtao/int.CD8.S35.sce.merged.rds")

cd8.ref.df <- reducedDim(cd8.ref.sce, "harmony.umap") %>% as.data.frame()
colnames(cd8.ref.df) <- paste0("UMAP_", 1:2)
cd8.ref.df$meta.cluster <- cd8.ref.sce$meta.cluster
cd8.ref.df$cluster.name <- plyr::mapvalues(cd8.ref.df$meta.cluster,
                                           from = names(color.set$meta.cluster)[1:41],
                                           to = names(color.set$cluster.name))

## a list of dgCMatrix object
expr.list <- readRDS("input/ZhengLiangtao/CD8.expr.list.rds")
sapply(expr.list, class)
sapply(expr.list, dim) # ~4500 genes kept in each dataset. To be noted, the genes are not same between datasets.

## 创建Seurat对象
seu.list <- lapply(expr.list, CreateSeuratObject)
seu.ref.cd8 <- base::Reduce(merge, seu.list)

## 写入cellmeta信息
cellmeta <- readRDS("input/ZhengLiangtao/CD8.cellmeta.rds")
colnames(cellmeta)[22:23] <- paste0("UMAP_", 1:2)
cellmeta <- cellmeta[rownames(seu.ref.cd8@meta.data), ]

umap.emb <- cellmeta[, paste0("UMAP_", 1:2)] %>% as.matrix()
cellmeta <- cellmeta[, -c(22,23)]

seu.ref.cd8@meta.data <- cellmeta
seu.ref.cd8[["umap"]] <- CreateDimReducObject(umap.emb, assay = "RNA")

seu.ref.cd8$meta.cluster.name <- plyr::mapvalues(
  x = seu.ref.cd8$meta.cluster,
  from = names(color.set$meta.cluster[1:41]),
  to = names(color.set$cluster.name))

## 写入metacell的信息
text.df <- data.frame(
  x = c(-6, -5, 2, 6,  3, -5, -4.5, 0),
  y = c( 1,  7, 6, 5, -7, -5, -1  , 0),
  label = c("ISG+", "Tex", "KIR+ NK-like", "Temra", "Tn", "Tc17", "Trm", "Tem")
)

seu.ref.cd8@misc$metacell <- list(
  meta.data = cd8.ref.df,
  colors = color.set$cluster.name[levels(cd8.ref.df$cluster.name)],
  text.pos = text.df
)

## 写入marker基因, 基因按照显著性排序
top.genes <- readRDS("input/ZhengLiangtao/CD8.markers.rds")
seu.ref.cd8@misc$markers <- top.genes

qs::qsave(seu.ref.cd8, "input/ZhengLiangtao.CD8.seurat.qs")
```