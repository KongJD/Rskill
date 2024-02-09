## Seurat多样本

```properties
## 背景
多样本分析流程
https://satijalab.org/seurat/archive/v3.0/immune_alignment.html
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
data <- list(scobj1,scobj2)
scobj <- merge(x=data[[1]], y = data[-1])

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

#### 2.样本进行批次矫正

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
ElbowPlot(scobj)
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30, reduction.name = "umap")

p1 <- DimPlot(scobj, reduction = "umap_native", group.by = "group")
p2 <- DimPlot(scobj, reduction = "umap", group.by = "group")
p1 + p2
```

