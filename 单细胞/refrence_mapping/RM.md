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
cellmeta <- cellmeta[rownames(seu.ref.cd4@meta.data),]

umap.emb <- cellmeta[, paste0("UMAP_", 1:2)] %>% as.matrix()
cellmeta <- cellmeta[, -c(22, 23)]

seu.ref.cd4@meta.data <- cellmeta
seu.ref.cd4[["umap"]] <- CreateDimReducObject(umap.emb, assay = "RNA")

seu.ref.cd4$meta.cluster.name <- plyr::mapvalues(
  x = seu.ref.cd4$meta.cluster,
  from = names(color.set$meta.cluster[1:41]),
  to = names(color.set$cluster.name))

## 写入metacell的信息
text.pos.df <- data.frame(
  x = c(-4.5, 4, 3.5, 2, -4.5, -7, -6),
  y = c(5, 4, -4, -6, -6, -2, 2),
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
cellmeta <- cellmeta[rownames(seu.ref.cd8@meta.data),]

umap.emb <- cellmeta[, paste0("UMAP_", 1:2)] %>% as.matrix()
cellmeta <- cellmeta[, -c(22, 23)]

seu.ref.cd8@meta.data <- cellmeta
seu.ref.cd8[["umap"]] <- CreateDimReducObject(umap.emb, assay = "RNA")

seu.ref.cd8$meta.cluster.name <- plyr::mapvalues(
  x = seu.ref.cd8$meta.cluster,
  from = names(color.set$meta.cluster[1:41]),
  to = names(color.set$cluster.name))

## 写入metacell的信息
text.df <- data.frame(
  x = c(-6, -5, 2, 6, 3, -5, -4.5, 0),
  y = c(1, 7, 6, 5, -7, -5, -1, 0),
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

#### 2. refference mapping 基于Seurat

```R
library(Seurat)
library(tidyverse)
#' 前后的Monocyte投影到DISCO Blood Atlas(v1.1)中的Monocyte参考图谱上。

#### 数据清洗(Reference) ####
#### https://www.immunesinglecell.org/atlasList
seu.ref <- readRDS("input/disco_blood_v1.1.rds")
seu.ref <- CreateSeuratObject(
  counts = seu.ref[["RNA"]]@counts,
  meta.data = seu.ref@meta.data)

## 修改cellmeta
seu.ref$disease[is.na(seu.ref$disease)] <- "control"
seu.ref$celltype <- seu.ref$ct
seu.ref$ct <- NULL
d <- c(
  "mild" = "covid-mild",
  "Mild COVID" = "covid-mild",
  "Moderate" = "covid-mild",
  "severe" = "covid-severe",
  "Severe" = "covid-severe",
  "Severe COVID" = "covid-severe"
)
seu.ref$disease_subtype <- d[seu.ref$disease_subtype]

## 选择Monocyte
seu.ref <- subset(seu.ref, celltype %in% c("CD14 monocyte", "CD14/CD16 monocyte", "CD16 monocyte"))

## 舍弃细胞数过少的研究(< 500 cells)
proj.stat <- table(seu.ref$project_id)
proj.used <- names(proj.stat)[proj.stat > 500]
seu.ref <- subset(seu.ref, project_id %in% proj.used)

qs::qsave(seu.ref, "input/disco_blood.mono.seurat.raw.qs")

#### 构建参考图谱 ####
seu.ref <- qs::qread("input/disco_blood.mono.seurat.raw.qs")
obj.list <- SplitObject(seu.ref, split.by = "project_id")

for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
}
## Anchored-CCA
system.time({
  anchors <- FindIntegrationAnchors(
    object.list = obj.list,
    dims = 1:20)
})

## Integration:
## 返回一个经过矫正的基因表达矩阵，存储在integrated assay中
## IntegrateData函数，看源代码
system.time({
  seu.ref <- IntegrateData(
    anchorset = anchors,
    dims = 1:20)
})

## Dimension reduction
DefaultAssay(seu.ref) <- "integrated"
seu.ref <- ScaleData(seu.ref, verbose = FALSE)
seu.ref <- RunPCA(seu.ref, npcs = 20, verbose = FALSE)
seu.ref <- RunUMAP(seu.ref, reduction = "pca",
                   dims = 1:20, return.model = TRUE,
                   verbose = FALSE)

p1 <- DimPlot(seu.ref, reduction = "umap",
              group.by = "celltype", label = T) +
  ggsci::scale_color_d3()
p2 <- DimPlot(seu.ref, reduction = "umap",
              group.by = "disease", label = F) +
  ggsci::scale_color_d3()
DimPlot(seu.ref, reduction = "umap",
        group.by = "disease_subtype", label = F)
DimPlot(seu.ref, reduction = "umap",
        cells.highlight = rownames(subset(seu.ref@meta.data, disease == "thrombocytopenia syndrome")),
        label = F)
p1 + p2


## 保存Reference data
qs::qsave(anchors, "input/disco_blood.mono.seurat.anchors.qs")
qs::qsave(seu.ref, "input/disco_blood.mono.seurat.reference.qs")

#### 参考映射 ####
seu.q <- readRDS("input/infb-pbmc.seurat.rds")
seu.q <- subset(seu.q, celltype %in% c("CD14 Mono", "CD16 Mono"))

## 参考映射分两步
## step1: Find anchors between query and ref
system.time({
  anchors <- FindTransferAnchors(
    reference = seu.ref,
    query = seu.q,
    dims = 1:20,
    reference.reduction = "pca",
    verbose = T)
})

## step2: map query to reference (fast!)
seu.q <- MapQuery(
  anchorset = anchors,
  reference = seu.ref,
  query = seu.q,
  reference.reduction = "pca",
  reduction.model = "umap")

p1 <- DimPlot(seu.q, reduction = "ref.umap", group.by = "group")
p2 <- DimPlot(seu.q, reduction = "ref.umap", group.by = "celltype")
p1 + p2

## 可视化
source("R/plot_projection.R")
seu.q1 <- subset(seu.q, group == "CTRL")
seu.q2 <- subset(seu.q, group == "STIM")

p1 <- PlotProjection(
  ref = seu.ref, query = seu.q1, query.reduction = "ref.umap",
  labels.col = "disease_subtype", ref.alpha = 1) +
  ggtitle("CTRL")
p2 <- PlotProjection(
  ref = seu.ref, query = seu.q2, query.reduction = "ref.umap",
  labels.col = "disease_subtype", ref.alpha = 1) +
  ggtitle("STIM")
p1 + p2

#### 补充说明 ####
#‘ 在整合效果不好的时候，可以考虑参考映射
# https://satijalab.org/seurat/articles/integration_mapping.html
```

#### 3. refference mapping 基于harmony

```R
library(symphony)
library(Seurat)
library(tidyverse)
#' phony(based on Harmony)流程提供的Reference Mapping接口，
#' 将INFB刺激前后的Monocyte投影到DISCO Blood Atlas(v1.1)中的Monocyte参考图谱上。
#### 3. 构建参考图谱 ####
seu.ref <- qs::qread("input/disco_blood.mono.seurat.raw.qs")
seu.ref <- NormalizeData(seu.ref)
## 3.1 鉴定高可变基因(HVGs)
obj.list <- SplitObject(seu.ref, split.by = "project_id")

for (i in 1:length(obj.list)) {
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
}
vfeatures <- SelectIntegrationFeatures(obj.list, nfeatures = 2000)
seu.ref[["RNA"]]@var.features <- vfeatures

## 3.2 Harmony整合
source("R/symphony_utils.R")
## 3.2a. 记录Scale data的相关参数(mean, sd)用于后续参考映射
## `CalMeanSds()` + `DoScale()` = `ScaleData()`
vfeatures_means_sds <- CalMeanSDs(seu.ref)
seu.ref <- DoScale(seu.ref, vfeatures_means_sds = vfeatures_means_sds)
## 3.2b. PCA
seu.ref <- RunPCA(seu.ref, npcs = 20, verbose = FALSE)

## 3.2c. Harmony (return with full harmony object): ~ 259s (4.3 min)
system.time({
  set.seed(0)
  Z_pca_ref <- seu.ref[["pca"]]@cell.embeddings
  ref_metadata <- seu.ref@meta.data
  ref_harmObj <- harmony::HarmonyMatrix(
    data_mat = Z_pca_ref,     ## PCA embedding matrix of cells
    meta_data = ref_metadata, ## data.frame with cell labels
    vars_use = c('project_id'),   ## variable to integrate out
    nclust = NULL,            ## number of clusters in Harmony model
    max.iter.harmony = 10,
    return_object = TRUE,     ## return the full Harmony model object
    do_pca = FALSE,           ## don't recompute PCs
  )
})

## 3.3 将harmony对象转变为symphony对象，用于参考映射 (~ 1 min)
system.time({
  vfeature_loadings <- seu.ref[["pca"]]@feature.loadings
  save_uwot_path <- file.path(getwd(), "input/disco_blood.mono.symphony.uwot_model")
  reference <- symphony::buildReferenceFromHarmonyObj(
    ref_harmObj,            # output object from HarmonyMatrix()
    ref_metadata,           # reference cell metadata
    vfeatures_means_sds,    # gene names, means, and std devs for scaling
    vfeature_loadings,      # genes x PCs matrix
    verbose = TRUE,         # verbose output
    do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
    save_uwot_path = save_uwot_path) # save_uwot_path should be full path.
})
## 保存参考模型
saveRDS(reference, 'input/disco_blood.mono.symphony.reference.rds')

## 将harmony的结果保存到Seurat对象中
seu.ref[["harmony"]] <- CreateDimReducObject(
  embeddings = t(reference$Z_corr), key = "harmony_", assay = "RNA")

seu.ref[["umap"]] <- CreateDimReducObject(
  embeddings = reference$umap$embedding, key = "UMAP_", assay = "RNA")

p1 <- DimPlot(seu.ref, reduction = "umap",
              group.by = "celltype", label = T) +
  ggsci::scale_color_d3()
p2 <- DimPlot(seu.ref, reduction = "umap",
              group.by = "disease", label = F) +
  ggsci::scale_color_d3()

p1 + p2

DimPlot(seu.ref, reduction = "umap",
        group.by = "disease_subtype", label = F)
DimPlot(seu.ref, reduction = "umap",
        cells.highlight = rownames(subset(seu.ref@meta.data, disease == "thrombocytopenia syndrome")),
        label = F)


#### 4. 参考映射 ####
seu.q <- readRDS("input/infb-pbmc.seurat.rds")
seu.q <- subset(seu.q, celltype %in% c("CD14 Mono", "CD16 Mono"))

source("R/symphony_utils.R")
seu.q <- MapQuery.Symphony(seu.q, reference = reference, assay.q = "RNA")

p1 <- DimPlot(seu.q, reduction = "ref.umap", group.by = "group")
p2 <- DimPlot(seu.q, reduction = "ref.umap", group.by = "celltype")
p1 + p2

## 可视化
source("R/plot_projection.R")
seu.q1 <- subset(seu.q, group == "CTRL")
seu.q2 <- subset(seu.q, group == "STIM")

p1 <- PlotProjection(
  ref = seu.ref, query = seu.q1, query.reduction = "ref.umap",
  labels.col = "disease_subtype", ref.alpha = 1) +
  ggtitle("CTRL")
p2 <- PlotProjection(
  ref = seu.ref, query = seu.q2, query.reduction = "ref.umap",
  labels.col = "disease_subtype", ref.alpha = 1) +
  ggtitle("STIM")
p1 + p2

# https://github.com/immunogenomics/symphony
# https://www.nature.com/articles/s41467-021-25957-x
```

#### 4.reference mapping 基于ProjecTILs

```R
### 下载参考图谱
# human CD8+ TIL atlas:
# https://doi.org/10.6084/m9.figshare.21931875.v2

# human CD4+ TIL atlas:
# https://doi.org/10.6084/m9.figshare.21981536.v1

# human blood and tumor-infiltrating DC atlas:
# https://doi.org/10.6084/m9.figshare.22040801.v1

# mouse TIL atlas:
# https://doi.org/10.6084/m9.figshare.12478571

# mouse acute and chronic viral infection CD8 T cell atlas:
# https://doi.org/10.6084/m9.figshare.12489518

# mouse acute and chronic viral infection CD4 T cell atlas:
# https://doi.org/10.6084/m9.figshare.16592693

rm(list = ls())
library(ProjecTILs)
library(STACAS)
library(Seurat)
library(SignatuR)
library(UCell)
library(scGate)
library(tidyverse)
library(patchwork)
setwd(here::here())

#### 参考图谱的数据结构 ####
seu.ref <- load.reference.map(ref = "input/Atlas-ProjectTILs/ref_TILAtlas_mouse_v1.rds")
misc <- seu.ref@misc

ref.cols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(seu.ref, label = T, cols = ref.cols)

#### 案例1: 示例数据 ####
### 1. Quick start
# 读入Query数据
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86031
seu.q <- ProjecTILs::query_example_seurat

# Reference mapping：注意query.projected的数据结构
# Run.ProjectTILs()
# => 1. make.projection()
# => 2. cellstate.predict()
query.projected <- Run.ProjecTILs(
  query = seu.q,
  ref = seu.ref,
  filter.cells = F) # Try to set it TRUE

T_model <- gating_model(name = "T", signature = c("Cd3d", "Cd3e", "Cd3g", "Cd2"))
query.projected <- Run.ProjecTILs(
  query = seu.q,
  ref = seu.ref,
  filter.cells = T, # Try to set it TRUE
  scGate_model = T_model)


# 可视化(密度图)
plot.projection(seu.ref, query.projected, linesize = 0.5, pointsize = 0.5)

# 可视化(柱状图)
plot.statepred.composition(seu.ref, query.projected, metric = "Percent")

# 可视化(marker基因的比较: query vs reference)
plot.states.radar(seu.ref, query = query.projected, min.cells = 30)

# 可视化(program的比较: query vs reference)
programs <- GetSignature(SignatuR$Mm$Programs)
seu.ref <- AddModuleScore_UCell(seu.ref, features = programs,
                                assay = "RNA", name = NULL)
query.projected <- AddModuleScore_UCell(query.projected, features = programs,
                                        assay = "RNA", name = NULL)
plot.states.radar(seu.ref, query = query.projected, meta4radar = names(programs))

# 细胞亚群注释(Label transfer)
seu.q <- seu.q %>%
  FindVariableFeatures(nfeatures = 500) %>%
  ScaleData() %>%
  RunPCA(npcs = 10) %>%
  RunUMAP(dims = 1:10)

DimPlot(seu.q)

# Label transfer in PCA space
seu.q <- ProjecTILs.classifier(
  query = seu.q,
  ref = seu.ref,
  filter.cells = F,
  reduction = "pca",
  ndim = 15,
  labels.col = "functional.cluster")
sel.cols <- grepl("functional.cluster", colnames(seu.q@meta.data))
colnames(seu.q@meta.data)[sel.cols] <- "celltype.pred.pca"

# Label transfer in UMAP space
seu.q <- ProjecTILs.classifier(
  query = seu.q,
  ref = seu.ref,
  filter.cells = F,
  reduction = "umap",
  ndim = 2,
  labels.col = "functional.cluster")
sel.cols <- grepl("functional.cluster", colnames(seu.q@meta.data))
colnames(seu.q@meta.data)[sel.cols] <- "celltype.pred.umap"

palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd", "#777777")
names(palette) <- c(levels(seu.ref$functional.cluster), "NA")

p1 <- DimPlot(seu.q, group.by = "celltype.pred.pca", cols = palette)
p2 <- DimPlot(seu.q, group.by = "celltype.pred.umap", cols = palette)
p1 + p2

#### 案例2: 将人的数据投影到小鼠的图谱上 ####
# https://pubmed.ncbi.nlm.nih.gov/30388456/
# Q: 那些T细胞和黑色素瘤病人免疫检查点阻断(ICB)响应程度相关？
# 注意这是小鼠的T细胞图谱
seu.ref <- load.reference.map(ref = "input/Atlas-ProjectTILs/ref_TILAtlas_mouse_v1.rds")
# Query data为人的T细胞数据
seu.q <- qs::qread("input/SadeFeldman.seurat.qs")
# View(seu.q@meta.data)
seu.q <- subset(seu.q, Time == "Pre")

# Reference mapping: ~ 382s (6.3 min)
system.time({
  seu.q1 <- make.projection(query = seu.q, ref = seu.ref)
})

# 并行计算
system.time({
  seu.q.list <- SplitObject(seu.q, split.by = "Sample")
  seu.q.list <- make.projection(query = seu.q.list, ref = seu.ref, ncores = 10)
  seu.q2 <- base::Reduce(f = merge.Seurat.embeddings, x = seu.q.list)
})

# 可视化（密度图）
p1 <- plot.projection(ref = seu.ref, query = seu.q1, ref.alpha = 1)
p2 <- plot.projection(ref = seu.ref, query = seu.q2, ref.alpha = 1)
p1 + p2

# Label transfer
seu.q1 <- cellstate.predict(ref = seu.ref, query = seu.q1, reduction = "pca")
sel.cols <- grepl("functional.cluster", colnames(seu.q1@meta.data))
colnames(seu.q1@meta.data)[sel.cols] <- c("celltype.pred.pca", "celltype.conf.pca")

seu.q1 <- cellstate.predict(ref = seu.ref, query = seu.q1, reduction = "umap")
sel.cols <- grepl("functional.cluster", colnames(seu.q1@meta.data))
colnames(seu.q1@meta.data)[sel.cols] <- c("celltype.pred.umap", "celltype.conf.umap")

table(seu.q1$celltype.pred.pca, seu.q1$celltype.pred.umap)

p1 <- DimPlot(seu.q1, reduction = "umap", group.by = "celltype.pred.pca", label = T)
p2 <- DimPlot(seu.q1, reduction = "umap", group.by = "celltype.pred.umap", label = T)
(p1 + p2) & ggsci::scale_color_d3("category20")

# 比较ICB响应组与不响应组
seu.q1$functional.cluster <- seu.q1$celltype.pred.pca
seu.q.list <- SplitObject(seu.q1, split.by = "Response")

genes.show <- c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7",
                "Gzmb", "Pdcd1", "Havcr2", "Tox", "Entpd1", "Cxcr5",
                "Ifng", "Cxcl13", "Xcl1", "Itgae")

plot.states.radar(ref = seu.ref,
                  query = seu.q.list,
                  min.cells = 20,
                  genes4radar = genes.show)

p1 <- plot.projection(seu.ref, seu.q.list$Responder, ref.alpha = 1) +
  ggtitle("Responder")
p2 <- plot.projection(seu.ref, seu.q.list$`Non-responder`, ref.alpha = 1) +
  ggtitle("Non-responder")
p1 + p2

p3 <- plot.statepred.composition(seu.ref, seu.q.list$Responder, metric = "Percent") +
  ggtitle("Responder")
p4 <- plot.statepred.composition(seu.ref, seu.q.list$`Non-responder`, metric = "Percent") +
  ggtitle("Non-responder")
(p3 + p4) & theme(plot.title = element_text(hjust = .5, face = "bold"))


#### 案例3: 用ProjectTILs创建自己的参考图谱 ####
## 读入Reference
seu.ref <- qs::qread("input/disco_blood.mono.seurat.raw.qs")
## 通过`STACAS`进行数据整合,以project_id作为批次
obj.list <- SplitObject(seu.ref, split.by = "project_id")
system.time({
  seu.ref <- Run.STACAS(
    object.list = obj.list,
    dims = 1:20,
    anchor.features = 2000)
})

seu.ref <- RunUMAP(seu.ref, dims = 1:20)

p1 <- DimPlot(seu.ref, reduction = "umap",
              group.by = "celltype", label = T) +
  ggsci::scale_color_d3()
p2 <- DimPlot(seu.ref, reduction = "umap",
              group.by = "disease", label = F) +
  ggsci::scale_color_d3()
DimPlot(seu.ref, reduction = "umap",
        group.by = "disease_subtype", label = F)
DimPlot(seu.ref, reduction = "umap",
        cells.highlight = rownames(subset(seu.ref@meta.data, disease == "thrombocytopenia syndrome")),
        label = F)
p1 + p2

qs::qsave(seu.ref, "input/disco_blood.mono.seurat.qs")

## make.reference()
## => 1. 重新计算PCA: `prcomp.seurat()`
## => 2. 计算UMAP投影相关的参数
## 比较seu.ref.1@misc和seu.ref.2@misc
system.time({
  seu.ref.1 <- make.reference(
    ref = seu.ref,
    recalculate.umap = FALSE,
    atlas.name = "mono-atlas",
    annotation.column = "celltype")
})

system.time({
  seu.ref.2 <- make.reference(
    ref = seu.ref,
    umap.method = 'umap',
    recalculate.umap = TRUE,
    atlas.name = "mono-altas",
    annotation.column = "celltype")
})

## 保存Reference data
qs::qsave(seu.ref.1, "input/disco_blood.mono.projecTILs.origUMAP.qs")
qs::qsave(seu.ref.2, "input/disco_blood.mono.projecTILs.recalUMAP.qs")

## Query data
seu.q <- readRDS("input/infb-pbmc.seurat.rds")
seu.q <- subset(seu.q, celltype %in% c("CD14 Mono", "CD16 Mono"))
## 保存de novo UMAP
seu.q[["umap.orig"]] <- seu.q[["umap"]]

## Project on original UMAP
## Save to seu.q1[["umap"]]
system.time({
  seu.q1 <- make.projection(
    query = seu.q,
    filter.cells = F,
    ref = seu.ref.1)
})

## Project on recalculated UMAP
system.time({
  seu.q2 <- make.projection(
    query = seu.q,
    filter.cells = F,
    ref = seu.ref.2)
})

p1 <- DimPlot(seu.q, reduction = "umap", group.by = "group")
p2 <- DimPlot(seu.q, reduction = "umap", group.by = "celltype")
p1 + p2

# 可视化（密度图）
seu.q1 <- subset(seu.q, group == "CTRL")
seu.q2 <- subset(seu.q, group == "STIM")

p1 <- plot.projection(
  ref = seu.ref, query = seu.q1,
  labels.col = "celltype", ref.alpha = 1) +
  ggtitle("CTRL")
p2 <- plot.projection(
  ref = seu.ref, query = seu.q2,
  labels.col = "celltype", ref.alpha = 1) +
  ggtitle("STIM")
p1 + p2

# https://github.com/carmonalab/ProjecTILs
# https://carmonalab.github.io/ProjecTILs.demo/tutorial.html
# https://carmonalab.github.io/ProjecTILs_CaseStudies/
# https://www.nature.com/articles/s41467-021-23324-4
```









