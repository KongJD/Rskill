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
