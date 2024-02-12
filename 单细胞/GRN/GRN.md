## GRN 基因调控网络

```properties
## 背景
http://www.seqyuan.com/GRN.html
https://lishensuo.github.io/posts/bioinfo/010%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E8%BD%AC%E5%BD%95%E5%9B%A0%E5%AD%90pyscenic/
## 数据securt多样本数据
外源刺激 -> TF(GRN) -> 转录组的变化 -> 通路/生物学过程的变化
```

#### 1.pySCENIC

```R
## 准备pySCENIC的输入文件:
## 1. 基因表达矩阵, meta cell matrix (For grnboost, step1)
## 2. TF list, download from
###   (1) 自定义
###   (2) SCENIC: https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
###   (3) AnimalTFDB: http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/TF_list_final/Homo_sapiens_TF
## 3. 参考文件:
###   (1) motif ranking tracks (.feather, 二进制文件, rows are motifs, columns are genes, values are ranks)
###   (2) motif annotation     (.tbl, txt文件, motif -> TFs)
rm(list = ls())

library(Seurat)
library(tidyverse)

source("R/makeMetaCells.R")
if (!dir.exists("output")) {
  dir.create("output")
}

### 读入数据
## clusters were annotated after batch effect removal via harmony.
seu <- readRDS("data/infb-pbmc.seurat.rds")
DimPlot(seu, group.by = "celltype", split.by = "group") & ggsci::scale_color_d3("category20")

### 创建meta cell
## 按样本处理
seu.list <- SplitObject(seu, split.by = "group")
seu.list[[1]]@project.name <- "pbmc_stim"
seu.list[[2]]@project.name <- "pbmc_ctrl"

metacells.list <- lapply(seq_along(seu.list), function(ii) {
  makeMetaCells(
    seu = seu.list[[ii]],
    min.cells = 10,
    reduction = "umap",
    dims = 1:2,
    k.param = 10,
    cores = 1)
})

mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)

## meta cell
seu2 <- CreateSeuratObject(mc.mat)
seu2 <- NormalizeData(seu2)
FeatureScatter(seu, feature1 = "CD8A", feature2 = "CD8B") +
  FeatureScatter(seu2, feature1 = "CD8A", feature2 = "CD8B")


## 准备pySCENIC的输入文件
## (1) TF list文件(optional)，可以使用预定义的TF list，例如pySCENIC官方提供的，或者animalTFDB提供的。
motif2tfs <- data.table::fread("cisTarget_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "cisTarget_db/hsa_hgnc_tfs.motifs-v10.txt")

## (2) meta cell matrix (for step1): *.csv or *.loom

## (2.1) 过滤低表达基因
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5,]
## (2.2) 过滤不在cisTargetDB中的基因
cisdb <- arrow::read_feather("cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
mc.mat <- mc.mat[genes.use,]


## 保存为loom文件
loom <- SCopeLoomR::build_loom(
  file.name = "output/00-2.mc_mat_for_step1.loom",
  dgem = mc.mat,
  default.embedding = NULL
)
loom$close()

## 释放loom对象，允许其它文件占用loom文件
rm(loom)
gc()
```