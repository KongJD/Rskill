## Seurat多样本

```properties
## 背景
多样本分析流程
https://satijalab.org/seurat/archive/v3.0/immune_alignment.html
https://github.com/satijalab/seurat-wrappers
```

#### 1.读取多个样本数据并合并

```R
library(Seurat)
### 样本1
scdata1 <- data.table::fread('data/immune_control_expression_matrix.txt.gz', data.table = F)
rownames(scdata1) <- scdata1$V1
scdata1 <- scdata1[, -1]
scobj1 <- CreateSeuratObject(counts = scdata1,
                             project = "pbmc_stim",
                             min.cells = 3,
                             min.features = 200)
scobj1@meta.data$group = "STIM"

### 样本2
scdata2 <- data.table::fread("data/immune_stimulated_expression_matrix.txt.gz", data.table = F)
rownames(scdata2) <- scdata2$V1
scdata2 <- scdata2[, -1]
scobj2 <- CreateSeuratObject(counts = scdata2,
                             project = "pbmc_ctrl",
                             min.cells = 3,
                             min.features = 200)
scobj2@meta.data$group = "CTRL"

### 合并数据
data <- list(scobj1, scobj2)
scobj <- merge(x = data[[1]], y = data[-1])

### 也可以在GEO下载数据后自己调整
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
library(Matrix)
scdata1 <- Matrix::readMM("data/GSE96583/GSM2560248_2.1.mtx.gz")
gene_id1 <- read.table("data/GSE96583/GSE96583_batch2.genes.tsv.gz")
barcode1 <- read.table("data/GSE96583/GSM2560248_barcodes.tsv.gz")
rownames(scdata1) <- gene_id1$V2
colnames(scdata1) <- barcode1$V1

### 读取另外一个数据
scdata2 <- Matrix::readMM("data/GSE96583/GSM2560249_2.2.mtx.gz")
gene_id2 <- read.table("data/GSE96583/GSE96583_batch2.genes.tsv.gz")
barcode2 <- read.table("data/GSE96583/GSM2560249_barcodes.tsv.gz")
rownames(scdata2) <- gene_id2$V2
colnames(scdata2) <- barcode2$V1

### 如果需要Read10X 就需要三个文件,
### 改造一下,barcodes.tsv.gz, features.tsv.gz,matrix.mtx.gz
scdata <- Read10X(data.dir = "data/GSE96583/stim/")
scobj <- CreateSeuratObject(counts = scdata,
                            project = "pbmc_stim",
                            min.cells = 3,
                            min.features = 200)
scobj@meta.data$group = "STIM"


scdata <- Read10X(data.dir = "data/GSE96583/ctrl/")
scobj <- CreateSeuratObject(counts = scdata,
                            project = "pbmc_ctrl",
                            min.cells = 3,
                            min.features = 200)
scobj@meta.data$group = "CTRL"
```

#### 2.多样本之间进行批次矫正

```R
### 不批次矫正的结果 两个样本 "umap_native"和“umap”有明显的区分，矫正后样本融合在一起了
library(harmony)
library(SeuratWrappers)

scobj@meta.data$percent.mt <- PercentageFeatureSet(scobj, pattern = "^MT-")
VlnPlot(scobj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))


## 筛选
scobj <- subset(scobj, subset = nFeature_RNA > 200 &
  nFeature_RNA < 1500
  &
  percent.mt < 5)

## 标准化
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj), reduction.name = "pca")
ElbowPlot(scobj)
scobj <- RunUMAP(scobj, reduction = "pca", dims = 1:15, reduction.name = "umap_native")
## 批次矫正
scobj <- RunHarmony(scobj, reduction = "pca", group.by.vars = "group", reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:15, reduction.name = "umap")

p1 <- DimPlot(scobj, reduction = "umap_native", group.by = "group")
p2 <- DimPlot(scobj, reduction = "umap", group.by = "group")
p1 + p2


cobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30)
### 设置多个resolution选择合适的resolution
scobj <- FindClusters(scobj, resolution = seq(0.2, 1.2, 0.1))

library(clustree)
clustree(scobj)

## 上述看图后 resloution 结果为0.5 差不多 后面的结果跟0.5 没有区别
scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.5
Idents(scobj) <- "seurat_clusters"
DimPlot(scobj, reduction = "umap", label = T)
DimPlot(scobj, reduction = "umap", group.by = "group")
DimPlot(scobj, reduction = "umap", split.by = "group")

scobj@assays$RNA$scale.data <- as.matrix(scobj@assays$RNA$scale.data)
scobj@assays$RNA$scale.data <- matrix()

FeaturePlot(scobj, features = "MS4A1", order = TRUE, reduction = "umap")
scobj@reductions$umap_naive <- NULL
saveRDS(scobj, file = "out/harmoy.rds")
```

#### 3.分群

```R
### 先分大群, 
### http://bio-bigdata.hrbmu.edu.cn/CellMarker/
### B: "MS4A1", "CD79A"
### NK: "GNLY", "NKG7"
### T: "CD3E","CD8A","CD4","IL7R", 
### Monocyte: "CD14", "FCGR3A", "LYZ"
### DC "FCER1A"
### Megakaryocytes/Platelet: "PPBP"
### Erythrocytes: "HBB","HBA2"

marker_genes <- c("MS4A1", "CD79A", "CD19")
VlnPlot(scobj, features = marker_genes)
FeaturePlot(scobj, features = marker_genes, order = TRUE, ncol = 2)

marker_genes <- c("CD14", "FCGR3A", "LYZ")
VlnPlot(scobj, features = marker_genes)
FeaturePlot(scobj, features = marker_genes, order = TRUE, ncol = 2)

### 汇总画图
marker_genes <- c("MS4A1",
                  "GNLY", "NKG7",
                  "CD3E", "CD8A", "CD4", "IL7R",
                  "CD14", "FCGR3A", "LYZ",
                  "FCER1A",
                  "PPBP"
)
FeaturePlot(scobj, features = marker_genes, order = TRUE, ncol = 5)

## cDC
marker_genes <- c("CLEC9A", "ITGAM", "ITGAE", "FCER1A")
## pDC
marker_genes <- c("IL3RA", "HLA-DRA")

## T细胞以及B细胞激活
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4110661/
marker_genes <- c("CCR7", "SELL", "CREM", "CD69")
VlnPlot(scobj, features = marker_genes)
FeaturePlot(scobj, features = marker_genes, order = TRUE, ncol = 3)

### 不好确定的时候,换一种作图方式(上述的群不好分辨的时候)
library(Nebulosa)
marker_genes <- c("CCR7", "SELL", "CREM", "CD69")
plot_density(scobj, features = marker_genes) + plot_layout(ncol = 2)

### 找出所有的marker
### 还有哪些群不确定的找它排名靠前的marker ,查marker的文献，能查到其图谱就可以
scobj <- JoinLayers(scobj)
all_markers <- FindAllMarkers(scobj)

library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:15) %>%
  ungroup()

#### 拿到群的注释结果，进行marker对应的群注释
Idents(scobj) <- "seurat_clusters"
### 给每个群添加注释
scobj <- RenameIdents(scobj,
                      "0" = "CD14 Mono",
                      "1" = "CD4 Naive T",
                      "2" = "CD4 Memory T",
                      "3" = "CD16 Mono",
                      "4" = "B cell",
                      "5" = "CD8 T",
                      "6" = "NK",
                      "7" = "T activated",
                      "8" = "DC",
                      "9" = "B Activated",
                      "10" = "Mk",
                      "11" = "pDC",
                      "12" = "Mono/Mk Doublets",
                      "13" = "Eryth"
)
DimPlot(scobj, reduction = "umap", label = T)
scobj@meta.data$celltype = Idents(scobj)

### 加入筛选细胞
Idents(scobj) <- "seurat_clusters"
scobj_subset <- subset(scobj, subset = seurat_clusters != 12)
DimPlot(scobj_subset, reduction = "umap", label = T)

scobj_subset1 <- subset(scobj, subset = celltype != "Mono/Mk Doublets")
DimPlot(scobj_subset1, reduction = "umap", label = T)

## 不要小群
scobj_subset2 <- subset(scobj, subset = seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
DimPlot(scobj_subset2, reduction = "umap", label = T)

## 反选
scobj_subset3 <- subset(scobj, subset = celltype %in% c("Mono/Mk Doublets", "Mk", "Eryth"), invert = TRUE)
DimPlot(scobj_subset3, reduction = "umap", label = T)

## 自己确定 反选
`%notin%` <- Negate(`%in%`)
scobj_subset4 <- subset(scobj, celltype %notin% c("Mono/Mk Doublets", "pDC", "Mk"))
DimPlot(scobj_subset4, reduction = "umap", label = T)

## 保存注释结果
Idents(scobj) <- "celltype"
DimPlot(scobj, reduction = "umap", label = T)
#saveRDS(scobj,file = "output/hamony_annotaion.rds")

## 对分好的1群想再分
Idents(scobj) <- "seurat_clusters"
## snn是shared nearest neighbors
## nn是nearest neighbors
scobj1 <- FindSubCluster(
  scobj,
  cluster = "1",
  graph.name = "RNA_snn",
  subcluster.name = "RNA_snn_res.0.5_c1_sub",
  resolution = 0.5
)

Idents(scobj1) <- "RNA_snn_res.0.5_c1_sub"
scobj1 <- RenameIdents(scobj1,
                       "0" = "CD14 Mono",
                       "1_0" = "CD4 Naive T",
                       "1_1" = "CD4 Naive T",
                       "1_2" = "CD8 Naive T",
                       "2" = "CD4 Memory T",
                       "3" = "CD16 Mono",
                       "4" = "B cell",
                       "5" = "CD8 T",
                       "6" = "NK",
                       "7" = "T activated",
                       "8" = "DC",
                       "9" = "B Activated",
                       "10" = "Mk",
                       "11" = "pDC",
                       "12" = "Mono/Mk Doublets",
                       "13" = "Eryth"
)

DimPlot(scobj1, reduction = "umap", label = T)
scobj1@meta.data$celltype_sub = Idents(scobj1)
plot_density(scobj1, c("CD8A", "CCR7"), joint = TRUE) + plot_layout(nrow = 1)
### 要注意reduction指定的问题
### https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html
### https://satijalab.org/seurat/articles/integration_introduction.html
### https://satijalab.org/seurat/archive/v3.0/immune_alignment.html
```

