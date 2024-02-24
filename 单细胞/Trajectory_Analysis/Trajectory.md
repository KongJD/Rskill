## 单细胞轨迹分析

```properties
## 背景
https://www.nature.com/articles/s41587-019-0068-4
```

#### 1.seurat router

```R
library(Seurat)
library(tidyverse)
library(scutilsR)

samples <- list.dirs("data/10X_matrix/", full.names = F, recursive = F)

### pblapply 是给apply函数添加进度条
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/10X_matrix/", sn))
  ## replace str "_" 到"-"
  sn <- gsub("_", "-", sn) # 注意"_"在`CreateSeuratObject()`里有特殊的意义
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts)
  return(seu)
})
## 合并样本
seu <- base::Reduce(f = merge, x = seu.list)
seu$donor <- sub("HS-BM-(P[1-3])-cells-[1-4]", "\\1", seu$orig.ident)

#### Cell QC ####
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

VlnPlot(seu, features = c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = c(1000, 5000), linetype = "dashed", color = "blue")

VlnPlot(seu, features = c("percent.mt"), group.by = "orig.ident", pt.size = 0, log = T) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "blue")

## 这里做一个数据备份，可以让我们调整细胞过滤的指标
seu.backup <- seu
qs::qsave(seu.backup, "tmp/HS_BM.seurat.init.qs")

seu <- subset(seu.backup, nFeature_RNA >= 1000)
seu <- subset(seu, percent.mt <= 10)
seu.backup

#### Mark Doublets ####
seu <- scutilsR::MarkDoublets(seu, split.by = "orig.ident")
prop.table(table(seu$orig.ident, seu$DF.classifications), 1)
qs::qsave(seu, "tmp/HS_BM.seurat.filtered.mark_doublets.qs")

#### Seurat Routine ####
seu.list <- SplitObject(seu, split.by = "donor")

seu <- seu.list$P1
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

## 计算细胞周期
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)

## PCA降维
seu <- ScaleData(seu)
seu <- RunPCA(seu)

## 检查细胞周期基因的影响
## CD34: pan stem marker
## GATA1: Ery
## MPO: Mono
## CD79A: CLP
DimPlot(seu, reduction = "pca", group.by = "Phase")
FeaturePlot(seu, reduction = "pca", features = c("CD34", "GATA1", "MPO", "CD79A"))

## 去除细胞周期的影响
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"))
seu <- RunPCA(seu)
##
DimPlot(seu, reduction = "pca", group.by = "Phase")
FeaturePlot(seu, reduction = "pca", features = c("CD34", "GATA1", "MPO", "CD79A"))


## 查看数据的复杂程度
ElbowPlot(seu, reduction = "pca", ndims = 30)
xx <- cumsum(seu[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.9) # 33 PCs解释了90%的方差，假设10%的方差来自于噪声
ndim = 33

## tSNE/UMAP降维
## perplexity: 20 (defualt)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:ndim, perplexity = 20)
p1 <- DimPlot(seu, reduction = "tsne", group.by = "Phase") + NoLegend()

## try a larger perplexity!
seu <- RunTSNE(seu, reduction = "pca", dims = 1:ndim, perplexity = 150)
p2 <- DimPlot(seu, reduction = "tsne", group.by = "Phase") + NoLegend()

seu <- RunTSNE(seu, reduction = "pca", dims = 1:ndim, perplexity = 150)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:ndim, metric = "", n.neighbors = 100)

p1 <- DimPlot(seu, reduction = "pca", group.by = "Phase") + NoLegend()
p2 <- DimPlot(seu, reduction = "tsne", group.by = "Phase") + NoLegend()
p3 <- DimPlot(seu, reduction = "umap", group.by = "Phase")
p1 + p2 + p3

p1 <- DimPlot(seu, reduction = "pca", group.by = "DF.classifications") + NoLegend()
p2 <- DimPlot(seu, reduction = "tsne", group.by = "DF.classifications") + NoLegend()
p3 <- DimPlot(seu, reduction = "umap", group.by = "DF.classifications")
p1 + p2 + p3


#### Clustering ####
## On PCA
resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:ndim, k.param = 20)
seu <- FindClusters(seu, resolution = resolutions)
DimPlot(seu, reduction = "tsne", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T)
DimPlot(seu, reduction = "umap", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T)

#### Save results ####
qs::qsave(seu, "tmp/HS_BM_donor1.seurat.qs")
```

#### 2.注释

```R
library(Seurat)
library(tidyverse)
library(scutilsR)

#### Find all markers ####
seu$seurat_clusters <- seu$RNA_snn_res.0.6
Idents(seu) <- seu$seurat_clusters
table(seu$seurat_clusters)
prop.table(table(seu$seurat_clusters, seu$DF.classifications), 1)
DimPlot(seu, reduction = "tsne", label = T)
DimPlot(seu, reduction = "umap", label = T)

## 抽样，加快计算速度
seu.ds <- subset(seu, downsample = 200)
## 并行计算
all.markers <- scutilsR::mcFindAllMarkers(seu.ds, do.flatten = F, only.pos = T, n.cores = 1)
all.markers <- lapply(all.markers, function(xx) subset(xx, p_val_adj < 1e-6 & avg_log2FC > log2(1.5)))
all.markers.df <- do.call(rbind, all.markers)
table(all.markers.df$cluster)

top3.markers <- all.markers.df %>%
  group_by(cluster) %>%
  slice_head(n=3) %>%
  select(Gene.name.uniq)

DotPlot(seu, features = unique(top3.markers$Gene.name.uniq)) + RotatedAxis()

#### Enrich cell types ####
data("mca_hsa")
data("hcl_hsa")
t2g <- rbind(mca_hsa, hcl_hsa)

e.res <- pbapply::pblapply(all.markers, function(xx) {
  yy <- xx$Gene.name.uniq
  tmp <- clusterProfiler::enricher(yy, TERM2GENE = t2g)
  res <- tmp@result
  res$cluster <- xx$cluster[1]
  res
})
e.res.df <- do.call(rbind, e.res)

#### Assign cell type to cluster ####
FeaturePlot(seu, reduction = "tsne", features = "HBB")
FeaturePlot(seu, reduction = "tsne", features = "MPO")
FeaturePlot(seu, reduction = "tsne", features = "IRF8")
FeaturePlot(seu, reduction = "tsne", features = "CD79A")
FeaturePlot(seu, reduction = "tsne", features = "PRG2")
FeaturePlot(seu, reduction = "tsne", features = "IRF7")
FeaturePlot(seu, reduction = "tsne", features = "SELP")
FeaturePlot(seu, reduction = "tsne", features = "PPBP")

cluster.annot <- read_tsv("tmp/03-3.annotation.tsv")
seu$celltype <- plyr::mapvalues(x = as.character(seu$seurat_clusters),
                                from = as.character(cluster.annot$cluster),
                                to = cluster.annot$annotation)
DimPlot(seu, reduction = "tsne", group.by = "celltype", label = T)
DimPlot(seu, reduction = "umap", group.by = "celltype", label = T)

# pDC
FeaturePlot(seu, reduction = "tsne", features = "SCT")
# cDC
FeaturePlot(seu, reduction = "tsne", features = "FCER1A")

seu$celltype[seu$RNA_snn_res.0.8==13] <- "cDC"
```

#### 3.

```R
library(Seurat)
library(tidyverse)
source("R/DimReduction.R")

seu <- qs::qread("tmp/HS_BM_donor1.seurat.annotated.qs")

## data flow:
## raw count matrix -> normalized matrix -(feature selection: xx HVGs)-> scaled matrix
## scaled matrix -> PCA embeddings -> (KNN graph) -> SNN(k shared nearest neighbor) graph (cell x cell)
## cluster(Louvain/Leiden) based on SNN graph
## if integration (remove batch effect) -> corrected PCA embeddings -> SNN graph
## visualization: (corrected) PCA embeddings -(distance choice)-> (cell2cell graph) -> tSNE/UMAP
## trajectory: (corrected) PCA embeddings -(distance choice)-> (cell2cell graph) -> diffusion map(DM)/force directed graph(FDG)

## Why give a so many npcs?
## be sure that no batch effect in your data.
seu <- RunPCA(seu, npcs = 300)

## Force directed graph (FDG)
seu <- RunFDG(seu, reduction = "pca", dims = 1:300)

p1 <- DimPlot(seu, reduction = "fr", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p2 <- DimPlot(seu, reduction = "tsne", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p3 <- DimPlot(seu, reduction = "umap", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p1 + p2 + p3

DimPlot(seu, reduction = "fr", group.by = "celltype", label = T) + ggsci::scale_color_d3()
DimPlot(seu, reduction = "fr", group.by = "Phase", label = T) + ggsci::scale_color_d3()
FeaturePlot(seu, reduction = "fr", features = c("CD34", "MPO", "IRF8", "CD79A", "GATA1", "HBB"), ncol = 3)

## Try different dimensions
# p.list <- pbapply::pblapply(c(10,20,30,50,100,150,200,300), function(ndim) {
#   seu <- RunFDG(seu, reduction = "pca", dims = 1:ndim)
#   p <- DimPlot(seu, reduction = "fr", group.by = "celltype", label = T) + ggsci::scale_color_d3()
#   p <- p + ggtitle(glue::glue("ndim = {ndim}"))
#   return(p)
# })
# cowplot::plot_grid(plotlist = p.list, nrow = 2)

## Diffusion map
seu <- RunDiffusion(seu, reduction = "pca", dims = 1:300, verbose = T)
DimPlot(seu, reduction = "diffusion.map", dims = 3:4, group.by = "celltype", label = T) + ggsci::scale_color_d3()

## Run tSNE/UMAP on DCs
## 对于发育过程的数据，调大`perplexity`有惊喜
seu <- RunTSNE(seu, reduction = "diffusion.map", dims = 1:10, perplexity = 150,
               reduction.name = "dm.tsne", reduction.key = "DMtSNE_")
## 对于发育过程的数据，调大`n.neighbors`有惊喜
seu <- RunUMAP(seu, reduction = "diffusion.map", dims = 1:10, n.neighbors = 100,
               reduction.name = "dm.umap", reduction.key = "DMUMAP_")

p1 <- DimPlot(seu, reduction = "tsne", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p2 <- DimPlot(seu, reduction = "dm.tsne", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p1 + p2

p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p2 <- DimPlot(seu, reduction = "dm.umap", group.by = "celltype", label = T) + ggsci::scale_color_d3()
p1 + p2

resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1)
DimPlot(seu, reduction = "fr", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T) & NoLegend()

#### save data ####
qs::qsave(seu, "tmp/HS_BM_donor1.seurat.trajectory.qs")
## 注意`sceasy::convertFormat`默认只保留Seurat对象中的`data`矩阵(也就是log1p normalized matrix).
## 如何想要改变这一行为，可以使用`main_layer`参数，例如
sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "data", outFile = "tmp/04.HS_BM_donor1.log1pnorm.h5ad")
sceasy::convertFormat(seu, from = "seurat", to = "anndata", main_layer = "counts", outFile = "tmp/04.HS_BM_donor1.counts.h5ad")
seu2 <- seu[rownames(seu[["RNA"]]@scale.data), ]
sceasy::convertFormat(seu2, from = "seurat", to = "anndata", main_layer = "scale.data", outFile = "tmp/04.HS_BM_donor1.scaled.h5ad")

#### Reference ####
# Force directed graphs compared to other dimensionality reduction methods for scRNA-seq
# https://www.biostars.org/p/415885/


#### MISC ####
seu <- qs::qread("tmp/HS_BM_donor1.seurat.trajectory.qs")

root.cells <- CellSelector(FeaturePlot(seu, reduction = "fr", features = "CD34"))
root.cells <- "HS-BM-P1-cells-1_GTACTCCGTCCAAGTT-1"
DimPlot(seu, reduction = "fr", cells.highlight = root.cells)
DimPlot(seu, reduction = "dm.tsne", cells.highlight = root.cells)
DimPlot(seu, reduction = "tsne", cells.highlight = root.cells)

terminal.p <- c(
  "HS-BM-P1-cells-1_CGTCACTTCAGGCGAA-1",
  "HS-BM-P1-cells-1_CGTCACTTCGCTTAGA-1",
  "HS-BM-P1-cells-1_CTTACCGCATCTACGA-1",
  "HS-BM-P1-cells-2_ACAGCTAGTAAGGATT-1"
)
DimPlot(seu, reduction = "fr", cells.highlight = terminal.p)

terminal.cells <- c(
  "pDC" = "HS-BM-P1-cells-1_CGTCACTTCGCTTAGA-1",
  "cDC" = "HS-BM-P1-cells-2_CCTTCCCCAGGGATTG-1",
  "CLP" = "HS-BM-P1-cells-1_CTTACCGCATCTACGA-1",
  "Ery" = "HS-BM-P1-cells-2_ACAGCTAGTAAGGATT-1",
  "Mega" = "HS-BM-P1-cells-1_CGTCACTTCAGGCGAA-1",
  "Mono" = "HS-BM-P1-cells-1_TTAGTTCTCGGACAAG-1"
)
# terminal.cells = "HS-BM-P1-cells-1_CTTACCGCATCTACGA-1"
DimPlot(seu, reduction = "fr", cells.highlight = terminal.cells)
# DimPlot(seu, reduction = "dm.tsne", cells.highlight = terminal.cells)
# DimPlot(seu, reduction = "tsne", cells.highlight = terminal.cells)

CellSelector(FeaturePlot(seu, reduction = "fr", features = "IRF8"))
CellSelector(DimPlot(seu, reduction = "fr", group.by = "celltype"))
```