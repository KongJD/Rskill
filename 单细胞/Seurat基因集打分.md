## 单细胞基因集打分

```properties
## 背景
基因集打分
```

#### 1.数据准备

```R
## 数据准备
rm(list = ls())

library(Seurat)
scdata <- Read10X(data.dir = "./data/10x/filtered_gene_bc_matrices/hg19/")
scobj <- CreateSeuratObject(counts = scdata,
                            project = "pbmc3k",
                            min.cells = 3,
                            min.features = 200)
library(dplyr)
scobj <- scobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:10) %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.5)
DimPlot(scobj, reduction = "umap", label = TRUE)
Idents(scobj) <- "seurat_clusters"
scobj <- RenameIdents(scobj,
                      "0" = "Naive CD4+ T",
                      "1" = "CD14+ Mono",
                      "2" = "Memory CD4+",
                      "3" = "B cell",
                      "4" = "CD8+ T",
                      "5" = "FCGR3A+ Mono",
                      "6" = "NK",
                      "7" = "DC",
                      "8" = "Platelet"
)
scobj@meta.data$celltype = Idents(scobj)
DimPlot(scobj, reduction = "umap", label = T, repel = T) + NoLegend()
```

#### 2.Seuarat基因集打分

```R
## AddModuleScore 函数
## 需要是一群基因
## 以NK的marker举例
library(tibble)
nk_enriched <- FindMarkers(scobj, ident.1 = "NK") %>%
  arrange(-avg_log2FC) %>%
  rownames_to_column(var = "gene")
nk_enriched_top <- nk_enriched$gene[1:50]
## AddModuleScore给的是基因名称就可打分
## 结果保存在 scobj@meta_data里面，行名为NK_enriched1
scobj <- AddModuleScore(scobj, features = list(nk_enriched_top),
                        name = "NK_enriched")
FeaturePlot(scobj, features = "NK_enriched1", label = TRUE, repel = TRUE)

library(ggplot2)
library(RColorBrewer)
FeaturePlot(scobj, features = "NK_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


### 多个特征
library(Seurat)
library(clusterProfiler)
genesets <- read.gmt("data/h.all.v2022.1.Hs.symbols.gmt")
signatures <- split(genesets$gene, genesets$term)
scobj <- AddModuleScore(scobj,
                        features = signatures,
                        name = "seurat",
                        assay = "RNA")

names(scobj@meta.data)[grep("seurat\\d", names(scobj@meta.data))] <- names(signatures)
FeaturePlot(scobj, features = "HALLMARK_INTERFERON_ALPHA_RESPONSE", label = TRUE, repel = TRUE)
```

#### 3. Ucell打分

```R
library(UCell)
library(Seurat)
library(clusterProfiler)
genesets <- read.gmt("data/h.all.v2022.1.Hs.symbols.gmt")
signatures <- split(genesets$gene, genesets$term)

scobj <- AddModuleScore_UCell(scobj, features = signatures, name = "_UCell")

signature.names <- paste0(names(signatures), "_UCell")
DimPlot(scobj, reduction = "umap", split.by = "group")
marker_genes <- "HALLMARK_INTERFERON_ALPHA_RESPONSE_UCell"
FeaturePlot(scobj, features = marker_genes, order = T)

library(ggplot2)
FeaturePlot(scobj, features = marker_genes, order = T) +
  labs(title = "Interferon alpha")
FeaturePlot(scobj, features = marker_genes, split.by = "group", order = T) ##右边处理组偏高

## 修改其中一列的名称
original_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE_UCell"
expected_name = "Interferon alpha"
scobj_new = scobj
location = grep(original_name, colnames(scobj_new@meta.data))
colnames(scobj_new@meta.data)[location] = expected_name
FeaturePlot(scobj_new, features = expected_name, split.by = "group", order = T)

library(patchwork)
plots <- VlnPlot(scobj, features = marker_genes, split.by = "group", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1) + labs(title = "Interferon alpha")
```

#### 4.aucell 打分

```R
library(AUCell)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)

celltype.markers <- FindAllMarkers(scobj, assay = "RNA", only.pos = T)
signatures <- celltype.markers %>%
  filter(p_val_adj < 0.01) %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(seq(n() * .1))

signatures <- split(signatures$gene, signatures$cluster)


### AUCell评分并分配活性阈值
exprMatrix <- scobj@assays$RNA@data
# 表达矩阵转rank矩阵
cells_rankings <- AUCell_buildRankings(exprMatrix, useNames = FALSE)
# min      1%      5%     10%     50%    100% 
# 212.00  325.00  434.95  539.90  816.00 3400.00
#以上数据意味着排名3400以后的基因没有价值

# 参数aucMaxRank的设置最好参考上面的结果，此处先不调整默认参数
cells_AUC <- AUCell_calcAUC(signatures, cells_rankings,
                            nCores = 1,
                            aucMaxRank = nrow(cells_rankings) * 0.05)
# 给每个基因集分配活性阈值
#cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
## 提取AUC, 可视化作图
sigmat <- getAUC(cells_AUC)
sigmat <- t(as.data.frame(sigmat))
sigmat <- sigmat[rownames(scobj@meta.data),]
colnames(sigmat) <- paste0(colnames(sigmat), "_AUCell")
scobj@meta.data <- cbind(scobj@meta.data, sigmat)

marker_genes <- colnames(scobj@meta.data)[grep("_AUCell", colnames(scobj@meta.data))]
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 3)
DimPlot(scobj, label = T)
## ssGSEA(慢)
library(GSVA)
exprMatrix <- scobj@assays$RNA@data
sigmat <- gsva(exprMatrix, signatures, method = "ssgsea", parallel.sz = 6)

sigmat <- t(as.data.frame(sigmat))
sigmat <- sigmat[rownames(scobj@meta.data),]
colnames(sigmat) <- paste0(colnames(sigmat), "_ssGSEA")
scobj@meta.data <- cbind(scobj@meta.data, sigmat)

marker_genes <- colnames(scobj@meta.data)[grep("_ssGSEA", colnames(scobj@meta.data))]
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 3)

## GSVA(十分慢)
exprMatrix <- scobj@assays$RNA@data
sigmat <- gsva(exprMatrix, signatures, method = "gsva", kcdf = "Gaussian", parallel.sz = 6)
sigmat <- t(as.data.frame(sigmat))
sigmat <- sigmat[rownames(scobj@meta.data),]
colnames(sigmat) <- paste0(colnames(sigmat), "_gsva")
scobj@meta.data <- cbind(scobj@meta.data, sigmat)
marker_genes <- colnames(scobj@meta.data)[grep("_gsva", colnames(scobj@meta.data))]
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 3)

## AddModuleScore
scobj <- AddModuleScore(scobj, features = signatures, name = "seurat", assay = "RNA")
names(scobj@meta.data)[grep("seurat\\d", names(scobj@meta.data))] <- paste0(names(signatures), "_seurat")

marker_genes <- colnames(scobj@meta.data)[grep("_seurat", colnames(scobj@meta.data))]
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 3)

## UCell
## 使用make.name 对features 的名称进行确认, 特殊符号会变成点
library(UCell)
scobj <- AddModuleScore_UCell(scobj, features = signatures, name = "_UCell")

marker_genes <- colnames(scobj@meta.data)[grep("_UCell", colnames(scobj@meta.data))]
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 3)


### 相关性分析
metaNames = colnames(scobj@meta.data)
myMetaNames <- metaNames[grep("FCGR3A", metaNames)]
FeaturePlot(scobj, features = myMetaNames, order = T, ncol = 3)

### 两两相关性
M <- cor(scobj@meta.data[, myMetaNames])
library(corrplot)
corrplot(M, method = "circle")
corrplot(M, method = "pie")

### 两个因素相关性
FeatureScatter(scobj, feature1 = 'CD9', feature2 = 'CD3E')
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(scobj, feature1 = "Naive.CD4.T_UCell", feature2 = "Memory.CD4.T_UCell")
FeatureScatter(scobj, feature1 = "FCGR3A+ Mono_AUCell", feature2 = "FCGR3A..Mono_UCell")

## 修改名称
original_name = "FCGR3A+ Mono_AUCell"
expected_name = "FCGR3A_Mono_AUCell"
scobj_new = scobj
location = grep(original_name, colnames(scobj_new@meta.data))
location = grep(original_name, colnames(scobj_new@meta.data), fixed = T)
colnames(scobj_new@meta.data)[location] = expected_name
FeatureScatter(scobj_new, feature1 = "FCGR3A_Mono_AUCell", feature2 = "FCGR3A..Mono_UCell")

```