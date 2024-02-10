## 单细胞基因集打分

```properties
## 背景
基因集打分
```

#### 1.

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
