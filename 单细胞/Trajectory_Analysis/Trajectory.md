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
  slice_head(n = 3) %>%
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

seu$celltype[seu$RNA_snn_res.0.8 == 13] <- "cDC"
```

#### 3.轨迹展示

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
seu2 <- seu[rownames(seu[["RNA"]]@scale.data),]
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

#### 4. 定义路径

```R
#### Load packages ####
library(Seurat)
library(tidyverse)
library(igraph)
setwd(here::here())

#' 定义lineage

#### Load data ####
seu <- qs::qread("tmp/HS_BM_donor1.seurat.trajectory.qs")
seu <- subset(seu, celltype != "Mast cell")

## Load the palantir results
pr.res <- read.table("tmp/palantir_results.csv", sep = ",", header = T, row.names = 1)

## Write to Seurat object
seu@meta.data <- pr.res

DimPlot(seu, reduction = "fr", group.by = "seurat_clusters", label = T)

#### Abstract Lineage: 最小生成树算法(Monocle V1) ####

## 过度聚类:
# k.clusters <- kmeans(x = Embeddings(seu, reduction = "diffusion.map")[,1:10], centers = 100)
# seu$seurat_clusters <- k.clusters$cluster
seu <- FindNeighbors(seu, reduction = "fr", dims = 1:2, k.param = 20)
seu <- FindClusters(seu, resolution = 8)

DimPlot(seu, reduction = "fr", group.by = "seurat_clusters", label = T) + NoLegend()

## the median cells for each cluster
data.use <- FetchData(seu, vars = c("seurat_clusters", "FDG_1", "FDG_2"))
data.use <- data.use %>%
  mutate(cellID = rownames(.)) %>%
  group_by(seurat_clusters) %>%
  mutate(x = abs(FDG_1 - median(FDG_1)),
         y = abs(FDG_2 - median(FDG_2)),
         d = x + y) %>%
  arrange(seurat_clusters, d) %>%
  slice_head(n = 1)

DimPlot(seu, reduction = "fr", cells.highlight = data.use$cellID)


##
cell.embeddings <- Embeddings(seu, reduction = "fr")[data.use$cellID,]
cell.pair.dis <- dist(cell.embeddings)
g <- igraph::graph_from_adjacency_matrix(adjmatrix = as.matrix(cell.pair.dis),
                                         mode = "undirected",
                                         weighted = TRUE,
                                         add.colnames = "cellID")

nodes2cellID <- get.vertex.attribute(g)$cellID
cellID2nodes <- setNames(1:length(nodes2cellID), nodes2cellID)

## 最小生成树
mst <- igraph::minimum.spanning.tree(g)
edges <- igraph::ends(mst, igraph::E(mst))
## Other spanning tree
# cell.pair.dis <- dist(cell.embeddings)
# dis.threshold <- 1
# cell.pair.dis <- as.matrix(cell.pair.dis)
# cell.pair.dis[cell.pair.dis > dis.threshold] <- 0
# g2 <- igraph::graph_from_adjacency_matrix(adjmatrix = cell.pair.dis,
#                                           mode = "undirected",
#                                           weighted = TRUE,
#                                           add.colnames = "cellID")
# edges <- igraph::ends(g2, igraph::E(g2))

edges <- as.data.frame(edges)
edges$from <- nodes2cellID[edges$V1]
edges$to <- nodes2cellID[edges$V2]
edges$from.x <- cell.embeddings[edges$from, 1]
edges$from.y <- cell.embeddings[edges$from, 2]
edges$to.x <- cell.embeddings[edges$to, 1]
edges$to.y <- cell.embeddings[edges$to, 2]

DimPlot(seu, reduction = "fr", cells.highlight = data.use$cellID) +
  geom_segment(inherit.aes = F, data = edges,
               mapping = aes(x = from.x, y = from.y, xend = to.x, yend = to.y), alpha = 1)


#### Lineage choice: 最短路径算法 ####
root.cells <- "HS-BM-P1-cells-1_GTACTCCGTCCAAGTT-1"
terminal.cells <- c(
  "pDC" = "HS-BM-P1-cells-1_CGTCACTTCGCTTAGA-1",
  "cDC" = "HS-BM-P1-cells-2_CCTTCCCCAGGGATTG-1",
  "CLP" = "HS-BM-P1-cells-1_CTTACCGCATCTACGA-1",
  "Ery" = "HS-BM-P1-cells-2_ACAGCTAGTAAGGATT-1",
  "Mega" = "HS-BM-P1-cells-1_CGTCACTTCAGGCGAA-1",
  "Mono" = "HS-BM-P1-cells-1_TTAGTTCTCGGACAAG-1"
)

nodes2cluster <- seu$seurat_clusters[nodes2cellID]
cluster2nodes <- setNames(names(nodes2cluster), nodes2cluster)

terminal.nodes <- cellID2nodes[cluster2nodes[seu$seurat_clusters[terminal.cells]]]
names(terminal.nodes) <- names(terminal.cells)

root.nodes <- cellID2nodes[cluster2nodes[seu$seurat_clusters[root.cells]]]

result <- igraph::shortest_paths(mst, from = root.nodes, to = terminal.nodes) ## modify me
lineages <- result$vpath
names(lineages) <- names(terminal.nodes)

DimPlot(seu, reduction = "fr", cells.highlight = nodes2cellID[lineages$Mega])

## cellID2lineage
seu$lineage.Mega <- seu$seurat_clusters %in% nodes2cluster[lineages$Mega]
DimPlot(seu, reduction = "fr", group.by = "lineage.Mega")

## for loop
for (i in seq_along(lineages)) {
  ln <- names(lineages)[i]
  meta.name <- glue::glue("lineage.{ln}")
  seu[[meta.name]] <- seu$seurat_clusters %in% nodes2cluster[lineages[[i]]]
}

DimPlot(seu, reduction = "fr", group.by = paste0("lineage.", names(lineages)), ncol = 3) &
  scale_color_manual(values = c("grey", "red")) &
  NoLegend()


#### Fine tune cells along the lineage: Principle curve ####
source("R/lineage.R")

## fit the principle curve
p.curve.Mega <- fit_pc(seu, lineage = "lineage.Mega", reduction = "fr", sample.n = 100) # 4
p.curve.pDC <- fit_pc(seu, lineage = "lineage.pDC", reduction = "fr") # 7
p.curve.cDC <- fit_pc(seu, lineage = "lineage.cDC", reduction = "fr") # 7
p.curve.CLP <- fit_pc(seu, lineage = "lineage.CLP", reduction = "fr", sample.n = 100) # 5
p.curve.Ery <- fit_pc(seu, lineage = "lineage.Ery", reduction = "fr") # 5
p.curve.Mono <- fit_pc(seu, lineage = "lineage.Mono", reduction = "fr") # 2


## plot on FR
data.point <- FetchData(seu, vars = c(paste0("FDG_", 1:2), "celltype"))
data.path.Mega <- get_path(p.curve.Mega, df = 5)
data.arrow.Mega <- get_arrow(data.path.Mega, reverse = T)

data.path.pDC <- get_path(p.curve.pDC, df = 7)
data.arrow.pDC <- get_arrow(data.path.pDC, reverse = F)

data.path.cDC <- get_path(p.curve.cDC, df = 7)
data.arrow.cDC <- get_arrow(data.path.cDC, reverse = F)

data.path.CLP <- get_path(p.curve.CLP, df = 5)
data.arrow.CLP <- get_arrow(data.path.CLP, reverse = F)

data.path.Ery <- get_path(p.curve.Ery, df = 5)
data.arrow.Ery <- get_arrow(data.path.Ery, reverse = T)

data.path.Mono <- get_path(p.curve.Mono, df = 2)
data.arrow.Mono <- get_arrow(data.path.Mono, reverse = T)

ggplot() +
  geom_point(data = data.point, aes(FDG_1, FDG_2, color = celltype), size = .2) +

  geom_path(data = data.path.Mega, aes(X, Y), size = .8) +
  geom_segment(data = data.arrow.Mega, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +

  geom_path(data = data.path.pDC, aes(X, Y), size = .8) +
  geom_segment(data = data.arrow.pDC, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +

  geom_path(data = data.path.cDC, aes(X, Y), size = .8) +
  geom_segment(data = data.arrow.cDC, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +

  geom_path(data = data.path.CLP, aes(X, Y), size = .8) +
  geom_segment(data = data.arrow.CLP, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +

  geom_path(data = data.path.Ery, aes(X, Y), size = .8) +
  geom_segment(data = data.arrow.Ery, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +

  geom_path(data = data.path.Mono, aes(X, Y), size = .8) +
  geom_segment(data = data.arrow.Mono, aes(x = X, xend = Xend, y = Y, yend = Yend),
               arrow = arrow(length = unit(0.1, "in"), angle = 30, type = "closed"), size = .5) +

  theme_classic(base_size = 15)


## cells around the curve
seu <- cells_on_lineage(seu, lineage = "Mega", reduction = "fr", data.path = data.path.Mega, delta = 0.1)
DimPlot(seu, reduction = "fr", group.by = "lineage.finetune.Mega") +
  scale_color_manual(values = c("grey", "red"))

seu <- cells_on_lineage(seu, lineage = "pDC", reduction = "fr", data.path = data.path.pDC)
seu <- cells_on_lineage(seu, lineage = "cDC", reduction = "fr", data.path = data.path.cDC)
seu <- cells_on_lineage(seu, lineage = "CLP", reduction = "fr", data.path = data.path.CLP)
seu <- cells_on_lineage(seu, lineage = "Ery", reduction = "fr", data.path = data.path.Ery)
seu <- cells_on_lineage(seu, lineage = "Mono", reduction = "fr", data.path = data.path.Mono)

## 结合lineage probability信息
seu$lineage.finetune.Mega <- seu$lineage.finetune.Mega | seu$Mega > 0.5
seu$lineage.finetune.cDC <- seu$lineage.finetune.cDC | seu$cDC > 0.5
seu$lineage.finetune.pDC <- seu$lineage.finetune.pDC | seu$pDC > 0.5
seu$lineage.finetune.Ery <- seu$lineage.finetune.Ery | seu$Ery > 0.5
seu$lineage.finetune.CLP <- seu$lineage.finetune.CLP | seu$CLP > 0.5
seu$lineage.finetune.Mono <- seu$lineage.finetune.Mono | seu$Mono > 0.5

DimPlot(seu, reduction = "fr", group.by = paste0("lineage.finetune.", names(terminal.cells)), cols = 3) &
  scale_color_manual(values = c("grey", "red")) &
  NoLegend()

qs::qsave(seu, "tmp/HS_BM_donor1.seurat.lineage.qs")
```

#### 5.differential potential

```R
library(Seurat)
library(tidyverse)
library(splines)
library(mgcv)

seu <- qs::qread("tmp/HS_BM_donor1.seurat.lineage.qs")

## branch plot
data.use <- FetchData(seu, vars = c("pseudotime", "celltype", "Ery"))
ggplot(data.use, aes(pseudotime, Ery, color = celltype)) +
  geom_point(size = .3) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggsci::scale_color_d3()


## Key concept: DP (differential potential)
data.use <- FetchData(seu, vars = c("DP", "pseudotime", "celltype"))

ggplot(data.use, aes(pseudotime, DP, color = celltype, group = celltype)) +
  geom_point(size = .5) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

## Early differentiation
seu.early <- subset(seu, celltype %in% c("HSC", "precursor"))
data.use <- FetchData(seu.early, vars = c("DP", "pseudotime", "celltype"))
data.use <- arrange(data.use, pseudotime)

ggplot(data.use, aes(pseudotime, DP, color = celltype, group = celltype)) +
  geom_point(size = .5) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

### fit the curve 线性可加模型(LAM)
model <- gam(DP ~ s(pseudotime), data = data.use) ## via `mgcv` package
data.use$fitted <- model$fitted.values

## 一阶微分:找到DP变化最大的位置
data.use$diff <- abs(c(NA, diff(data.use$fitted) / diff(data.use$pseudotime)))
## Max-min标准化到0-1区间
data.use$diff <- (data.use$diff - min(data.use$diff, na.rm = T)) / (max(data.use$diff, na.rm = T) - min(data.use$diff, na.rm = T))

## 二阶微分：确定极值点
data.use$diff.diff <- abs(c(NA, diff(data.use$diff) / diff(data.use$pseudotime)))

ggplot(data.use, aes(pseudotime, DP, color = celltype, group = celltype)) +
  geom_point(size = .5) +
  geom_vline(xintercept = c(0.069, 0.156)) +
  geom_point(aes(pseudotime, fitted), color = "red", size = .1) +
  geom_point(aes(pseudotime, diff), color = "blue", size = .1) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

data.DP <- data.use

## Key molecular events
source("R/compute_module_score.R")
pathway <- readLines("resource/h.all.v2022.1.Hs.symbols.gmt")
gene.list <- lapply(strsplit(pathway, split = "\t"), function(xx) xx[-c(1, 2)])
names(gene.list) <- sapply(strsplit(pathway, split = "\t"), function(xx) xx[1])

## 17s
system.time({
  seu.early <- ComputeModuleScore(seu.early, gene.sets = gene.list, min.size = 20, cores = 10)
})
DefaultAssay(seu.early) <- "AUCell"

data.use <- FetchData(seu.early, vars = c(rownames(seu.early), "pseudotime", "DP", "celltype"))
dim(data.use)

colnames(data.use) <- gsub("-", ".", colnames(data.use))
colnames(data.use) <- sub("^HALLMARK.", "", colnames(data.use))
dvars <- setdiff(colnames(data.use), c("pseudotime", "DP", "celltype"))
dvars <- tolower(dvars) %>% Hmisc::capitalize()
colnames(data.use) <- c(dvars, c("pseudotime", "DP", "celltype"))

cor.test <- cor(data.use[, dvars], data.use$DP)

cor.res.df <- data.frame(
  term = rownames(cor.test),
  corr = cor.test[, 1]
)

cor.res.df <- arrange(cor.res.df, desc(corr))
cor.res.df$rank <- 1:nrow(cor.res.df)

ggplot(cor.res.df, aes(rank, corr)) +
  geom_point()

ggplot(data.use, aes(pseudotime, `Oxidative.phosphorylation`)) +
  geom_point(size = .5, aes(color = celltype)) +
  geom_point(inherit.aes = F, data = data.DP, aes(pseudotime, diff / 3), size = .2, color = "red") +
  geom_smooth() +
  # scale_color_viridis_c() +
  geom_vline(xintercept = c(0.069, 0.156)) +
  theme_classic(base_size = 15)

## Metabolism remodel (有氧呼吸 上调)
```

#### 6.DEGS within lineage

```R
library(Seurat)
library(tidyverse)
library(mgcv)

setwd(here::here())

#### Load data ####
seu <- qs::qread("tmp/HS_BM_donor1.seurat.lineage.qs")

## Imputation by `MAGIC`(https://github.com/KrishnaswamyLab/MAGIC)
expr.in.cells <- Matrix::rowSums(seu[["RNA"]]@counts > 0)
select.features <- names(expr.in.cells[expr.in.cells >= 50])
seu.magic <- Rmagic::magic(seu, genes = select.features, seed = 1024, npcs = 30)
DefaultAssay(seu.magic) <- "MAGIC_RNA"

## Ery lineage
seu.Ery <- subset(seu.magic, lineage.finetune.Ery)

## check the results
DefaultAssay(seu.Ery) <- "RNA"
p1 <- FeaturePlot(seu.Ery, reduction = "fr", features = "HBB")
DefaultAssay(seu.Ery) <- "MAGIC_RNA"
p2 <- FeaturePlot(seu.Ery, reduction = "fr", features = "HBB")
(p1 | p2) & scale_color_viridis_c()

### DP
data.use <- FetchData(seu.Ery, vars = c("DP", "pseudotime", "celltype"))
data.use <- arrange(data.use, pseudotime)

ggplot(data.use, aes(pseudotime, DP, color = celltype, group = celltype)) +
  geom_point(size = .5) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

### fit the curve 线性可加模型(LAM)
model <- gam(DP ~ s(pseudotime), data = data.use) ## via `mgcv` package
data.use$fitted <- model$fitted.values

## 一阶微分:找到DP变化最大的位置
data.use$diff <- abs(c(NA, diff(data.use$fitted) / diff(data.use$pseudotime)))
## Max-min标准化到0-1区间
data.use$diff <- (data.use$diff - min(data.use$diff, na.rm = T)) / (max(data.use$diff, na.rm = T) - min(data.use$diff, na.rm = T))

## 二阶微分：确定极值点
data.use$diff.diff <- abs(c(NA, diff(data.use$diff) / diff(data.use$pseudotime)))

ggplot(data.use, aes(pseudotime, DP, color = celltype, group = celltype)) +
  geom_point(size = .5) +
  geom_vline(xintercept = c(0.152, 0.382)) +
  geom_point(aes(pseudotime, fitted), color = "red", size = .1) +
  geom_point(aes(pseudotime, diff), color = "blue", size = .1) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

data.DP <- data.use

### related genes
data.use <- FetchData(seu.Ery, vars = c(rownames(seu.Ery), "pseudotime", "Ery"))
cor.test <- cor(data.use[, 1:(length(data.use) - 2)], data.use$Ery)

cor.res.df <- data.frame(
  term = rownames(cor.test),
  corr = cor.test[, 1]
)

FeaturePlot(seu, reduction = "fr", features = c("Ery", "KLF1", "GATA1", "MPST"), ncol = 2) &
  scale_color_viridis_c()

FeaturePlot(seu.magic, reduction = "fr", features = c("Ery", "KLF1", "GATA1", "MPST"), ncol = 2) &
  scale_color_viridis_c()

## LMM model
seu.Ery$pseudotime.bin <- infotheo::discretize(seu.Ery$pseudotime, disc = "equalwidth", 100)[[1]]
expr.in.cells <- Matrix::rowSums(seu.Ery[["RNA"]]@counts > 0)
select.features <- names(expr.in.cells[expr.in.cells >= 50])
imputed.data <- FetchData(seu.Ery, vars = select.features)
vd.vars <- c("pseudotime.bin")
meta.data <- seu.Ery@meta.data[, vd.vars, drop = F]

dim(imputed.data) ## 10750 genes

## 48s
system.time({
  vd.res <- Gadest::VarDecompose(data = imputed.data, meta.data = meta.data, vd.vars = vd.vars, cores = 20)
})

vd.res$Ery.corr <- cor.res.df[vd.res$gene, "corr"]

ggplot(vd.res, aes(pseudotime.bin, Ery.corr)) +
  geom_point(size = .1)

## 感受一下方差分解的好处
head(subset(vd.res, abs(Ery.corr) < 0.1 & pseudotime.bin > 0.75))

FeaturePlot(seu.magic, reduction = "fr", features = c("Ery", "THRB", "RGS18", "GATA2"), ncol = 2) &
  scale_color_viridis_c()

##
all.DEGs <- subset(vd.res, pseudotime.bin > 0.5)
dim(all.DEGs)

## gene cluster:

## visualization
data.use <- FetchData(seu.Ery, vars = c("pseudotime", "DP", "KLF1", "GATA1", "celltype"))

ggplot(data.use, aes(pseudotime, GATA1)) +
  geom_point(aes(color = celltype), size = 1, alpha = .5) +
  geom_vline(xintercept = c(0.152, 0.382)) +
  geom_smooth() +
  geom_point(inherit.aes = F, data = data.DP, aes(pseudotime, fitted), size = .1, color = "red") +
  geom_point(inherit.aes = F, data = data.DP, aes(pseudotime, diff), size = .1, color = "orange") +
  theme_classic(base_size = 15)

```

#### 7.DEGS across lineage

```R
library(Seurat)
library(tidyverse)
library(patchwork)

seu <- qs::qread("tmp/HS_BM_donor1.seurat.lineage.qs")

## Imputation by `MAGIC`(https://github.com/KrishnaswamyLab/MAGIC)
expr.in.cells <- Matrix::rowSums(seu[["RNA"]]@counts > 0)
select.features <- names(expr.in.cells[expr.in.cells >= 50])
seu.magic <- Rmagic::magic(seu, genes = select.features, seed = 1024, npcs = 30)
DefaultAssay(seu.magic) <- "MAGIC_RNA"

lineages <- c("Mono", "CLP")
p1 <- FeaturePlot(seu, reduction = "fr", features = c(lineages), ncol = 3) &
  scale_color_viridis_c()

p2 <- DimPlot(seu, reduction = "fr", group.by = paste0("lineage.finetune.", lineages), ncol = 3) &
  scale_color_manual(values = c("grey", "red")) &
  NoLegend()

p1 / p2


seu2 <- subset(seu, lineage.finetune.Mono | lineage.finetune.CLP)
seu2 <- subset(seu2, pseudotime <= 0.42)
VlnPlot(seu2, group.by = "celltype", features = "pseudotime") +
  geom_hline(yintercept = 0.42)
DimPlot(seu2, reduction = "fr", group.by = "celltype")

## LMM model
seu2$pseudotime.bin <- infotheo::discretize(seu2$pseudotime, disc = "equalwidth", 100)[[1]]
seu2$is.Mono <- seu2$Mono > 0.5
seu2$is.CLP <- seu2$CLP > 0.5

DimPlot(seu2, reduction = "fr", group.by = c("is.Mono", "is.CLP"), ncol = 2)

expr.in.cells <- Matrix::rowSums(seu2[["RNA"]]@counts > 0)
select.features <- names(expr.in.cells[expr.in.cells >= 50])
imputed.data <- FetchData(seu2, vars = select.features)
vd.vars <- c("pseudotime.bin", "is.Mono", "is.CLP")
meta.data <- seu2@meta.data[, vd.vars, drop = F]

dim(imputed.data) ## 12131 genes

## 1.5 min
system.time({
  vd.res <- Gadest::VarDecompose(data = imputed.data, meta.data = meta.data, vd.vars = vd.vars, cores = 20)
})

ggplot(vd.res, aes(is.Mono, is.CLP)) +
  geom_point(size = .2) +
  theme_bw(base_size = 15)

top5.mono <- head(arrange(vd.res, desc(is.Mono))$gene)
top5.CLP <- head(arrange(vd.res, desc(is.CLP))$gene)
top5.seudo <- head(arrange(vd.res, desc(pseudotime.bin))$gene)

FeaturePlot(seu.magic, reduction = "fr", features = c(top5.mono, top5.CLP), ncol = 6) &
  scale_color_viridis_c()
```
