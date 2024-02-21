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

seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/10X_matrix/", sn))
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



