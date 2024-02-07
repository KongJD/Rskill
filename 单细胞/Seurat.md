## 基于Seurat包单细胞流程

```properties
## 背景：
https://satijalab.org/seurat/articles/get_started_v5_new
https://github.com/gao-lab/GLUE/
https://github.com/carmonalab/ProjecTILs
```

#### 1.基本步骤

```R
rm(list = ls())

library(Seurat)
scdata <- Read10X(data.dir = "./data/10x/filtered_gene_bc_matrices/hg19")
scobj <- CreateSeuratObject(counts = scdata, project = "pbmc",
                            min.cells = 3, min.features = 200)
library(dplyr)
scobj <- scobj %>%
  NormalizeData() %>%  #归一化
  FindVariableFeatures() %>%  #筛选特征,找高变基因
  ScaleData() %>%  #标准化
  RunPCA() %>% #降维
  RunUMAP(dims = 1:10) %>%  #非线性降维
  FindNeighbors() %>% #建立SNN
  FindClusters(resolution = 0.5)  #分群5

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
DimPlot(scobj, reduction = "umap", label = T, repel = T) + NoLegend()
```

#### 1. 质控、筛选、找高变基因、

```R
scobj <- CreateSeuratObject(counts = scdata, project = "pbmc",
                            min.cells = 3, min.features = 200)

### 数据质控
### 主要PercentageFeatureSet函数计算线粒体含量
### 人类使用pattern = "^MT-"，小鼠使用pattern = "^mt-"
scobj@meta.data$percent.mt = PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")

scobj@meta.data

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### 正式筛选，筛选的是细胞，最终细胞减少
### nFeature_RNA > 200
### nFeature_RNA < 2500
### percent.mt < 5
scobj <- subset(scobj, subset = nFeature_RNA > 200 &
  nFeature_RNA < 2500 &
  percent.mt < 5)

## 数据标准化: NormalizeData
### 先除以总数，再乘以系数，然后取log
scobj <- NormalizeData(scobj, normalization.method = "LogNormalize", scale.factor = 10000)

## 特征筛选(高变基因)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)

### 高变基因可视化
### 使用VariableFeatures 函数提取高变基因
top10 <- head(VariableFeatures(scobj), 10)
### 使用VariableFeaturePlot 画图
plot1 <- VariableFeaturePlot(scobj) ## 横坐标基因平均表达量，纵坐标为变异度
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


t1 <- cbind(as.data.frame(rownames(scobj@assays$RNA@features)), scobj@assays$RNA@meta.data)

### 自己画
t1 <- cbind(as.data.frame(rownames(scobj@assays$RNA@features)), scobj@assays$RNA@meta.data)
colnames(t1)[1] <- "gene_name"
rownames(t1) <- t1$gene_name
library(ggplot2)

ggplot(t1, aes(x = vf_vst_counts_mean, y = vf_vst_counts_variance.standardized, color = vf_vst_counts_variable)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Average Expression", y = "Standardized Variance")

##加标签
library(ggrepel)
ggplot(t1, aes(x = vf_vst_counts_mean, y = vf_vst_counts_variance.standardized, color = vf_vst_counts_variable)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Average Expression", y = "Standardized Variance") +
  geom_text_repel(data = t1[top10,], label = top10, color = "black")
```

#### 2.数据缩放、降维

```R
### 降维之前的必备操作
### 缩放的效果是，基因的平均值是0，方差是1
## 保存时把它删掉
scobj <- ScaleData(scobj, features = rownames(scobj))
scale_data <- scobj@assays$RNA@layers$scale.data
hist(apply(scale_data, 1, mean))
hist(apply(scale_data, 1, var))
hist(matrixStats::rowMeans2(scale_data))
hist(matrixStats::rowVars(scale_data))

## PCA线性降维

```

