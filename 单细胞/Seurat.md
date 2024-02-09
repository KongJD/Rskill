## 基于Seurat包单细胞流程

```properties
## 背景：
https://satijalab.org/seurat/articles/get_started_v5_new
https://github.com/gao-lab/GLUE/
https://github.com/carmonalab/ProjecTILs
https://pair-code.github.io/understanding-umap/
https://www.youtube.com/watch?v=FgakZw6K1QQ&t=684s&ab_channel=StatQuestwithJoshStarmer
```

#### 1.基本步骤

```R
rm(list = ls())

options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
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

#### 1. 质控、筛选、找高变基因

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
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj))
DimPlot(scobj, reduction = "pca")

##pca data
pca_data <- as.data.frame(scobj@reductions$pca@cell.embeddings)
ggplot(pca_data, aes(x = PC_1, y = PC_2, color = "red")) + geom_point()

### 选择合适的PCA维度
ElbowPlot(scobj) ##这里图片显示10之后pca就没啥变化了

##umap非线性降维，依赖pca结果
scobj <- RunUMAP(scobj, dims = 1:10)
umap_data <- as.data.frame(scobj@reductions$umap@cell.embeddings)
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = "red")) + geom_point()
```

#### 3.聚类分群

```R
### 聚类分群
### 找紧邻,dims = 1:10 跟UMAP相同
scobj <- FindNeighbors(scobj, dims = 1:10)
### 分群(需要在找近邻后)
scobj <- FindClusters(scobj, resolution = 0.5) ### 会在metadata中增加两列数据"RNA_snn_res.0.5" "seurat_clusters"
DimPlot(scobj, reduction = "umap", label = T) ## 将umap结果分群

metdaat <- scobj@meta.data
scobj <- FindClusters(scobj, resolution = seq(0.2, 1.2, 0.2)) ## 分辨率改为0.2到1.2 ，分群结果细致
DimPlot(scobj, reduction = "umap", label = T) ## 将umap结果分群

library(clustree)
clustree(scobj) ## 看不同分辨率下哪些群被分的 更细致了

### 选择特定分辨率得到的分群此处为RNA_snn_res.0.5
scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.5
DimPlot(scobj, reduction = "umap", label = T)
```

#### 4.分群注释

```R
### 分群注释
### 先分大群, marker哪里来
### http://bio-bigdata.hrbmu.edu.cn/CellMarker/ (此网站找细胞的marker基因)
### B: "MS4A1", "CD79A"
### NK: "GNLY", "NKG7"
### T: "CD3E","CD8A","CD4","IL7R", 
### Monocyte: "CD14", "FCGR3A", "LYZ"
### DC "FCER1A"
### Platelet: "PPBP"

"MS4A1" %in% rownames(scobj)
"CD79A" %in% rownames(scobj)
### 使用VlnPlot画marker小提琴图
VlnPlot(scobj, features = c("MS4A1", "CD79A"))
### 使用FeaturePlot画出特征分布图
FeaturePlot(scobj, features = c("MS4A1", "CD79A"), order = TRUE, ncol = 2)


VlnPlot(scobj, features = c("GNLY", "NKG7"))
FeaturePlot(scobj, features = c("GNLY", "NKG7"), order = TRUE, ncol = 2)

VlnPlot(scobj, features = c("CD3E", "CD8A", "CD4", "IL7R"))
FeaturePlot(scobj, features = c("CD3E", "CD8A", "CD4", "IL7R"), order = TRUE, ncol = 2)

VlnPlot(scobj, features = c("CD14", "FCGR3A", "LYZ"))
FeaturePlot(scobj, features = c("CD14", "FCGR3A", "LYZ"), order = TRUE, ncol = 2)

VlnPlot(scobj, features = c("FCER1A", "PPBP"))
FeaturePlot(scobj, features = c("FCER1A", "PPBP"), order = TRUE, ncol = 2)

marker_genes <- c("MS4A1", "CD79A",
                  "GNLY", "NKG7",
                  "CD3E", "CD8A", "CD4", "IL7R",
                  "CD14", "FCGR3A", "LYZ",
                  "FCER1A",
                  "PPBP"
)
FeaturePlot(scobj, features = marker_genes, order = TRUE, ncol = 5)

### 汇总画图
marker_genes <- c("MS4A1", "CD79A",
                  "GNLY", "NKG7",
                  "CD3E", "CD8A", "CD4", "IL7R",
                  "CD14", "FCGR3A", "LYZ",
                  "FCER1A",
                  "PPBP"
)
FeaturePlot(scobj, features = marker_genes, order = TRUE, ncol = 5)

### 再确定细节
marker_genes <- c("CD3E", "CD4", "CD8A",
                  "CCR7", "SELL", "CREM", "TCF7", "S100A4")
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 5)

marker_genes <- c("CCR7", "SELL", "CREM", "TCF7")
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 2)
marker_genes <- c("RPS12", "FASLG")
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 2)

### 不好确定的时候,换一种作图方式
library(Nebulosa)
marker_genes <- c("CCR7", "SELL", "TCF7", "S100A4")
plot_density(scobj, features = marker_genes) + plot_layout(ncol = 2)

### 确定不了的时候自己分析maker
all_markers <- FindAllMarkers(object = scobj)
all_markers$rank <- all_markers$pct.1 / all_markers$pct.2
## 这里得到了差异表达矩阵，后续会写一下对它单细胞的富集分析(见下一篇)
## all_markers 这里自己找的时候找特异性强的，尽量其他群不表达
## 这里找avg_logfc差异倍数大的,pct1大一点，pct2小点，这样特异性强
marker_genes <- c("RPS12", "FASLG")
FeaturePlot(scobj, features = marker_genes, order = T, ncol = 2)

library(dplyr)

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:10) %>%
  ungroup()
```

#### 5. 确定注释结果

```R
### 确认群的个数
Idents(scobj) <- "seurat_clusters"
head(Idents(scobj))
## 给每个群加注释
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
DimPlot(scobj, label = T)
head(Idents(scobj))
scobj@meta.data$celltype = Idents(scobj)
metadata <- scobj@meta.data
### 保存结果
saveRDS(scobj, file = "output/Seurat_single_sample_scobj.rds")
```

#### 6.注释结果可视化

```R
DimPlot(scobj, reduction = "umap", label = T)
DimPlot(scobj, reduction = "umap", label = T) + NoLegend()
scCustomize::DimPlot_scCustom(scobj, figure_plot = TRUE)

### 各个细胞群，画图展示
### 出现在scobj@meta.data 中的 列都可以FeaturePlot 展示
scobj@meta.data$cd14 <- ifelse(scobj@meta.data$celltype == "CD14+ Mono", 1, 0)
FeaturePlot(scobj, features = "cd14")

desgin <- model.matrix(~0 + metadata$celltype)
desgin <- as.data.frame(desgin)
colnames(desgin) <- levels(metadata$celltype)
scobj@meta.data <- cbind(scobj@meta.data, desgin)
FeaturePlot(scobj, features = levels(metadata$celltype))

### 可视化展示DotPlot
### 因子来限定细胞的顺序
Idents(scobj) <- factor(Idents(scobj), levels = c("B cell",
                                                  "NK",
                                                  "Naive CD4+ T",
                                                  "Memory CD4+",
                                                  "CD8+ T",
                                                  "CD14+ Mono",
                                                  "FCGR3A+ Mono",
                                                  "DC",
                                                  "Platelet"))

markers.to.plot <- c("MS4A1", "CD79A",
                     "GNLY", "NKG7",
                     "CD3E", "CD8A", "CD4", "IL7R",
                     "CD14", "FCGR3A", "LYZ",
                     "FCER1A",
                     "PPBP")
DotPlot(scobj, features = markers.to.plot, dot.scale = 8)
DotPlot(scobj, features = markers.to.plot, dot.scale = 8) + RotatedAxis()
DotPlot(scobj, features = markers.to.plot, dot.scale = 8) +
  coord_flip() +
  RotatedAxis()

### Clustered_DotPlot
### 可以自己挑选marker来呈现
### 也可以使用FindAllMarkers批量提取marker

all_markers <- FindAllMarkers(object = scobj)
library(dplyr)
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

DotPlot(scobj, features = top5_markers) +
  coord_flip() +
  RotatedAxis()
scCustomize::Clustered_DotPlot(scobj, features = top5_markers)

### marker 热图
DoHeatmap(scobj, features = top5_markers)

### Identity的大小修改，通过size参数
DoHeatmap(scobj, features = top5_markers, size = 3)

### 基因的大小修改theme(axis.text.y = element_text(size = 8))
DoHeatmap(scobj, features = top5_markers, size = 3) +
  theme(axis.text.y = element_text(size = 8))

### subset和downsample 可以随机取每个群的细胞数
DoHeatmap(subset(scobj, downsample = 50), features = top5_markers, size = 3) +
  theme(axis.text.y = element_text(size = 8))
```

#### geo读取

```R
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188987
### 1.分步读入
count_matrix <- Matrix::readMM("data/GSE188987/GSE188987_YJ_Lib9_1_matrix.mtx.gz")
gene_id <- read.table("data/GSE188987/GSE188987_YJ_Lib9_1_features.tsv.gz")
barcode <- read.table("data/GSE188987/GSE188987_YJ_Lib9_1_barcodes.tsv.gz")
rownames(count_matrix) <- gene_id$V2
colnames(count_matrix) <- barcode$V1

### 2.修改后Read10X读入
scdata1 <- Read10X(data.dir = "data/GSE188987_2/")
```

#### 画更好看的图

```R
### https://enblacar.github.io/SCpubr-book/index.html
### https://samuel-marsh.github.io/scCustomize/
### https://github.com/junjunlab/scRNAtoolVis
### https://github.com/zhanghao-njmu/SCP
### https://samuel-marsh.github.io/scCustomize

if (!requireNamespace("scRNAtoolVis", quietly = TRUE)) {
  install.packages("./resource/scRNAtoolVis-master/", type = "source", repos = NULL)
}

library(scRNAtoolVis)

scobj <- readRDS(file = "output/Seurat_scobj.rds")
all_markers <- readRDS(file = "output/Seurat_all_markers.rds")
library(dplyr)
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:5) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

AverageHeatmap(object = scobj,
               markerGene = top5_markers)


library(Seurat)
DimPlot(scobj, label = T)
metadata = scobj@meta.data
design = model.matrix(~0 + metadata$celltype)
design <- as.data.frame(design)
colnames(design) = levels(metadata$celltype)

scobj@meta.data <- cbind(scobj@meta.data, design)
FeaturePlot(scobj, features = levels(metadata$celltype))

### 上色
color_use <- scCustomize::DiscretePalette_scCustomize(num_colors = 9, palette = "varibow")

DimPlot(scobj, label = T, cols = color_use)
FeaturePlot(scobj, features = "Naive CD4+ T", cols = c("lightgrey", "#FF7373"))

### 批量操作
plotlist = list()
features = levels(metadata$celltype)
for (i in 1:9) {
  plotlist[[features[i]]] = FeaturePlot(scobj, features = features[i], cols = c("lightgrey", color_use[i])) + NoLegend()
}

library(patchwork)
wrap_plots(plotlist)
```


