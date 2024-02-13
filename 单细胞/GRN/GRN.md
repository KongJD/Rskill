## GRN 基因调控网络

```properties
## 背景
http://www.seqyuan.com/GRN.html
https://lishensuo.github.io/posts/bioinfo/010%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--%E8%BD%AC%E5%BD%95%E5%9B%A0%E5%AD%90pyscenic/
## 数据securt多样本数据
外源刺激 -> TF(GRN) -> 转录组的变化 -> 通路/生物学过程的变化
```

#### 1.pySCENIC数据准备

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

makeMetaCells <- function(seu, min.cells = 10, reduction = "umap", dims = 1:2, k.param = 10, cores = 1) {
  seu <- seu %>%
    FindNeighbors(reduction = reduction, dims = dims, k.param = k.param) %>%
    FindClusters(res = 50)
  metadata <- seu@meta.data
  metadata$METACELL_ID <- factor(metadata$seurat_clusters)
  dge_mat <- seu[["RNA"]]@counts

  dge_mat_mc <- parallel::mclapply(levels(metadata$METACELL_ID), function(xx) {
    cells <- rownames(subset(metadata, METACELL_ID == xx))
    Matrix::rowSums(dge_mat[, cells])
  }, mc.cores = cores)
  dge_mat_mc <- do.call(cbind, dge_mat_mc)

  metacell_metadata <-
    metadata[["METACELL_ID"]] %>%
      table() %>%
      as.data.frame()
  colnames(metacell_metadata) <- c("METACELL_ID", "CELL_COUNT")
  rownames(metacell_metadata) <- metacell_metadata[["METACELL_ID"]]

  kept.cells <- subset(metacell_metadata, CELL_COUNT >= min.cells)[["METACELL_ID"]]
  metacells <- list(
    mat = dge_mat_mc[, kept.cells],
    metadata = metacell_metadata[kept.cells,]
  )
  colnames(metacells$mat) <- paste0(seu@project.name, ".METACELL_", kept.cells)
  rownames(metacells$metadata) <- colnames(metacells$mat)
  metacells
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

### 在这里做metacell的好处就是在干扰素刺激下 STAT1和STAT2两个转录因子应该有很强的相关性，但原始数据没有，但是做了metacell就有了相关性
## 上述的代码暂时没能理解为啥做了，两个转录因子有这么强的相关性，后续过来理解下
FeatureScatter(seu, feature1 = "STAT1", feature2 = "STAT2") +
  FeatureScatter(seu2, feature1 = "STAT1", feature2 = "STAT2")

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

#### 2.pySCENIC流程

```R
## inputs
f_loom_grn = .. / output / 00 - 2.mc_mat_for_step1.loom

## outputs
grn_output=../output/step1_adj.tsv
ctx_output=../output/step2_reg.tsv

## reference
f_tfs=../cisTarget_db/hsa_hgnc_tfs.motifs-v10.txt
f_motif_path=../cisTarget_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
f_db_500bp=../cisTarget_db/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
f_db_10kb=../cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


#### 1. build GRN
pyscenic grn \
--seed 777 \
--num_workers 10 \
--method grnboost2 \
-o $grn_output \
$f_loom_grn \
/cisdb/$f_tfs
# time arboreto_with_multiprocessing.py \
# $f_loom_grn \
# $f_tfs \
# --method grnboost2 \
# --output $grn_output \
# --num_workers 10 \
# --seed 777


#### 2. cisTarget
time pyscenic ctx \
$grn_output \
$f_db_500bp $f_db_10kb \
--annotations_fname$f_motif_path\
--expression_mtx_fname $f_loom_grn\
--output$ctx_output\
--num_workers10

###### 这里推荐用 singularity ，后续学一下这个容器用法
```

```python
#### 上述结果转换格式
import sys
import pandas as pd
from pyscenic.cli.utils import load_signatures

# pyscenic ctx的输出文件
regulon_file = "/data/step2_reg.tsv"
project = "/data/ifnb_pbmc"
# 最小的regulon基因数，用于过滤regulon
min_regulon_size = 10


def get_motif_logo(regulon):
    base_url = "http://motifcollections.aertslab.org/v10nr_clust/logos/"
    for elem in regulon.context:
        if elem.endswith('.png'):
            return (base_url + elem)


sys.stdout.flush()
regulons = load_signatures(regulon_file)
select_cols = [i.name for i in regulons if len(i.genes) >= min_regulon_size]
gmt_file = project + ".regulons.gmt"
txt_file = project + ".regulons.txt"
fo1 = open(gmt_file, 'w')
fo2 = open(txt_file, 'w')
for i in regulons:
    if i.name in select_cols:
        motif = get_motif_logo(i)
        genes = "\t".join(i.genes)
        tf = "%s(%sg)" % (i.transcription_factor, len(i.genes))
        fo1.write("%s\t%s\t%s\n" % (tf, motif, genes))
        fo2.write("%s\t%s\t%s\n" % (tf, motif, genes.replace("\t", ",")))
```

#### 3.pySCENIC的结果探索:用各细胞的转录因子分析

#### （1）刺激后影响最大的细胞类型是什么

```R
library(Seurat)
library(tidyverse)
library(patchwork)

seu <- subset(seu, celltype %in% c("Mono/Mk Doublets", "Eryth"), invert = TRUE)
DimPlot(seu, group.by = "celltype", reduction = "umap", label = T)
celltype.levels <- c("CD14 Mono", "CD16 Mono", "DC", "pDC", "B cell", "B Activated",
                     "NK", "CD8 T", "CD4 Memory T", "T activated", "CD4 Naive T", "Mk")
seu$celltype <- factor(seu$celltype, levels = celltype.levels)

## 导入regulon (gene list)
regulons <- clusterProfiler::read.gmt("output/02-ifnb_pbmc.regulons.gmt")

## data.frame -> list, list中的每个元素为一个gene set
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
print(regulon.list[1])

ComputeModuleScore <- function(x, ...) UseMethod('ComputeModuleScore')

ComputeModuleScore.default <- function(x, gene.sets, min.size = 20, batch.size = 500, cores = 1, ...) {
  if (!is.list(gene.sets)) {
    stop("'gene.sets' should be a list or data.frame!")
  }
  gene.sets <- gene.sets[sapply(gene.sets, length) >= min.size]
  n.cells <- ncol(x)
  batches <- floor((1:n.cells - 1) / batch.size)
  batch.levels <- unique(batches)

  aucell <- function(i) {
    dge.tmp <- x[, batches == i]
    cr <- AUCell::AUCell_buildRankings(dge.tmp, nCores = 1, plotStats = F, verbose = F)
    auc <- AUCell::AUCell_calcAUC(gene.sets, cr, nCores = 1, verbose = F)
    AUCell::getAUC(auc)
  }

  auc_scores <- parallel::mclapply(batch.levels, aucell, mc.cores = cores)
  do.call(cbind, auc_scores)
}

ComputeModuleScore.Seurat <- function(x, gene.sets, min.size = 20, batch.size = 500, cores = 1, assay = Seurat::DefaultAssay(x), ...) {
  dge <- x[[assay]]@counts
  ras_mat <- ComputeModuleScore.default(x = dge, gene.sets, min.size, batch.size, cores)
  x[["AUCell"]] <- Seurat::CreateAssayObject(data = ras_mat)
  return(x)
}


## 用AUCell计算RAS matrix
## RAS = regulon activity score
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, min.size = 10, cores = 1)
seu
DefaultAssay(seu) <- "AUCell"

p1 <- FeaturePlot(seu, features = "STAT2(+)", split.by = "group")
p2 <- FeaturePlot(seu, features = "STAT2", split.by = "group")
(p1 / p2) & scale_color_viridis_c()

VlnPlot(seu, group.by = "celltype", features = "STAT2(+)", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = "STAT2", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red"))


## 用RAS matrix计算UMAP
seu <- RunUMAP(object = seu,
               features = rownames(seu),
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRAS",
               reduction.key = "umapRAS_")

## 可视化：UMAP on harmony
p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype") +
  ggsci::scale_color_d3("category20") +
  NoLegend()
p2 <- DimPlot(seu, reduction = "umap", group.by = "group") + NoLegend()

## 可视化：UMAP on RAS 
##(批次矫正后一定程度掩盖了组间差别，此方法在批次矫正后也一定程度上显示了组间的差别)
## 这里确保组间的差异为什么不是批次，因为基因打分是天然去批次，给的是基因名，而不是表达量
p3 <- DimPlot(seu, reduction = "umapRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p4 <- DimPlot(seu, reduction = "umapRAS", group.by = "group")

(p1 + p3) / (p2 + p4)

## 推测：INFB对髓系细胞的影响更大(cd14 和 cd16)

## 换一种方式：PCA
DefaultAssay(seu) <- "AUCell"
seu <- ScaleData(seu)
seu <- RunPCA(object = seu,
              features = rownames(seu),
              reduction.name = "pcaRAS",
              reduction.key = "pcaRAS_")

## 可视化：PCA on RAS
p3 <- DimPlot(seu, reduction = "pcaRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p4 <- DimPlot(seu, reduction = "pcaRAS", group.by = "group")
p3 + p4

## PC1 encoding the regulons related to cell type
## PC2 encoding the regulons affected by INFB treatment
## The INFB induced transcriptome shift is orthogonal to the cell identity transcriptional programs.

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_1", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red"))

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_2", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red"))
```

#### （2）哪些转录因子（TF）响应了外源刺激后导致样本的的转录组变化?

##### A:鉴定出差异的调控的转录因子，其在组间和在哪种细胞类型中变化最大

```R
library(Seurat)
library(tidyverse)
library(patchwork)

#### Method 1: feature loadings of PCA ####
p1 <- DimPlot(seu, reduction = "pcaRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p2 <- DimPlot(seu, reduction = "pcaRAS", group.by = "group")
p1 + p2

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_2", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red"))

## PC2-high-loading regulons
## gene expression matrix (cells x genes) = cell embeddings (cells x nPCs) %*% feature loadings (nPCs x genes)
## Here, only the high variable genes were involved.
## If you want to fetch all genes loadings, you should perform ProjectDim first
# seu <- ProjectDim(seu, reduction = "pcaRAS", do.center = T)
# We didn't need do this because we use all regulons to perform PCA.
Loadings(seu, reduction = "pcaRAS")
VizDimLoadings(seu, reduction = "pcaRAS", dims = 2, balanced = T, projected = F)
VlnPlot(seu, group.by = "celltype", features = "FOS(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red")) + ylab("TF activity")

#### Method 2: Variance Decomposition ####
VarDecompose <- function(data, meta.data, vd.vars, genes = "all", cores = -1) {
  ## check params
  if (missing(data) ||
    missing(meta.data) ||
    missing(vd.vars)) {
    stop("Must provide 'data', 'meta.data', and 'vd.vars'.")
  }
  if (is.null(colnames(meta.data)) || is.null(rownames(meta.data))) {
    stop("The row and column names of 'meta.data' should be provided.")
  }
  if (is.null(colnames(data)) || is.null(rownames(data))) {
    stop("The row and column names of 'data' should be provided.")
  }
  if (!all(rownames(data) == rownames(meta.data))) {
    stop("The row names of 'data' and 'meta.data' should be matched.")
  }
  if (!all(vd.vars %in% colnames(meta.data))) {
    vd.vars.404 <- setdiff(vd.vars, colnames(meta.data))
    stop(paste("vd.vars:", vd.vars.404, "is(are) not found in 'meta.data'"))
  }
  if (length(genes) == 1) {
    if (genes == "all") {
      genes <- colnames(data)
    } else {
      stop("'genes' should be 'all' or a vector.")
    }
  } else {
    genes.404 <- setdiff(genes, colnames(data))
    if (length(genes.404) > 0) {
      warning(paste(length(genes.404), "gene(s) are not found in 'data'."))
      genes <- setdiff(genes, genes.404)
    }
  }
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  ## prepare data for VD
  vd.vars.str <- sapply(vd.vars, function(xx) sprintf("(1|%s)", xx))
  modelFormulaStr <- paste("expression ~", paste(vd.vars.str, collapse = "+"))
  data.use <- cbind(data[, genes], meta.data)
  ## exe VD
  vd.res <- do.call(rbind, parallel::mclapply(genes, function(genename) {
    data.model <- data.use[, c(vd.vars, genename)]
    colnames(data.model) <- c(vd.vars, "expression")
    tryCatch({
      model <- suppressWarnings(lme4::lmer(stats::as.formula(modelFormulaStr), data = data.model, REML = TRUE, verbose = FALSE))
      results <- as.data.frame(lme4::VarCorr(model))
      rownames(results) <- results$grp
      results <- results[c(vd.vars, "Residual"),]
      frac.var <- results$vcov / sum(results$vcov)

      res.tmp <- c("OK", frac.var)
    },
      error = function(e) {
        print(e)
        res.tmp <- c("FAIL", rep(-1, length(vd.vars) + 1))
      })
    names(res.tmp) <- c("status", vd.vars, "residual")
    as.data.frame(as.list(res.tmp)) # return
  }, mc.cores = cores)
  )
  rownames(vd.res) <- genes
  vd.res %<>% as.data.frame()
  vd.res <- vd.res %>% dplyr::mutate(gene = rownames(vd.res), .before = 1)
  # vd.res <- vd.res %>% as.data.frame() %>% dplyr::mutate(gene = rownames(.), .before=1)
  for (i in 3:ncol(vd.res)) {
    vd.res[[i]] %<>% as.numeric()
  }
  return(vd.res)
}

DefaultAssay(seu) <- "AUCell"
vd.vars <- c("celltype", "group")
meta.data <- seu@meta.data[, vd.vars]
ras.data <- FetchData(seu, vars = rownames(seu))
vd.res <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)


## PCA和VD结果的一致性
vd.res$PC1.loadings <- Loadings(seu, reduction = "pcaRAS")[vd.res$gene, 1]
vd.res$PC2.loadings <- Loadings(seu, reduction = "pcaRAS")[vd.res$gene, 2]

to.label <- subset(vd.res, celltype > 0.25 & group > 0.15)

ggplot(vd.res, aes(celltype, group, color = abs(PC1.loadings))) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.15, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0.25, color = "blue", linetype = "dashed") +
  ggrepel::geom_text_repel(inherit.aes = F, data = to.label, aes(celltype, group, label = gene)) +
  xlab("Fraction of variance by cell type") +
  ylab("Fraction of variance by group") +
  scale_color_viridis_c() +
  theme_bw(base_size = 15)

ggplot(vd.res, aes(celltype, group, color = abs(PC2.loadings))) +
  geom_point(size = 2) +
  xlab("Fraction of variance by cell type") +
  ylab("Fraction of variance by group") +
  scale_color_viridis_c() +
  theme_bw(base_size = 15)

ggplot(vd.res, aes(group)) +
  geom_density()

to.label <- subset(vd.res, group > 0.4)
ggplot(vd.res, aes(PC2.loadings, group)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(inherit.aes = F, data = to.label, aes(PC2.loadings, group, label = gene)) +
  ylab("Fraction of variance by group") +
  theme_bw(base_size = 15)

VlnPlot(seu, group.by = "celltype", features = "HESX1(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red")) + ylab("TF activity")
FeaturePlot(seu, reduction = "umap", features = "HESX1(+)", split.by = "group")
VlnPlot(seu, group.by = "celltype", features = "ETV6(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red")) + ylab("TF activity")
FeaturePlot(seu, reduction = "umap", features = "ETV6(+)", split.by = "group")
VlnPlot(seu, group.by = "celltype", features = "IRF7(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red")) + ylab("TF activity")
FeaturePlot(seu, reduction = "umap", features = "IRF7(+)", split.by = "group")


#### Method 3: DE test ####
DERegulon <- function(seu, celltype, group, test.use = "wilcox") {
  seu$new.group <- paste(seu[[celltype, drop = T]], seu[[group, drop = T]], sep = "_")
  Idents(seu) <- factor(seu$new.group)

  celltypes <- levels(seu[[celltype, drop = T]])
  groups <- levels(seu[[group, drop = T]])
  seu2 <- subset(seu, group == groups[1])
  Idents(seu2) <- seu2[[celltype, drop = T]]
  baseline.levels <- AverageExpression(seu2, assays = "AUCell")$AUCell

  de.list <- lapply(celltypes, function(ct) {
    message(glue::glue("processing {ct} ..."))
    ct1 <- paste(ct, groups[1], sep = "_")
    ct2 <- paste(ct, groups[2], sep = "_")
    de <- FindMarkers(seu, ident.1 = ct1, ident.2 = ct2, test.use = "wilcox", fc.name = "avg_diff", logfc.threshold = 0)
    de$change <- ifelse(de$avg_diff > 0,
                        paste("up in", groups[1]),
                        paste("up in", groups[2]))
    de$avg_diff <- -de$avg_diff
    de$diff_rate <- de$avg_diff / baseline.levels[rownames(de), ct]
    de$group <- ct
    return(de)
  })
  names(de.list) <- celltypes
  de.list
}

#' regulon activity ~ (0,1)
format_change.df <- . %>%
  mutate(change = case_when(p_val_adj < 1e-6 & avg_diff > 0.3 ~ '+++',
                            p_val_adj < 1e-6 &
                              avg_diff > 0.1 &
                              avg_diff <= 0.3 ~ '++',
                            p_val_adj < 1e-6 &
                              avg_diff > 0.05 &
                              avg_diff <= 0.1 ~ '+',
                            p_val_adj < 1e-6 & avg_diff < -0.3 ~ '---',
                            p_val_adj < 1e-6 &
                              avg_diff < -0.1 &
                              avg_diff >= -0.3 ~ '--',
                            p_val_adj < 1e-6 &
                              avg_diff < -0.05 &
                              avg_diff >= -0.1 ~ '-',
                            TRUE ~ "")) %>%
  mutate(gene = rownames(.)) %>%
  dplyr::select(gene, change, group)

#'
format_change <- function(mylist, celltype.levels) {
  f.change <- lapply(mylist, format_change.df) %>%
    Reduce(rbind, .) %>%
    pivot_wider(names_from = "group", values_from = "change")
  f.change <- f.change %>% dplyr::select(c("gene", all_of(celltype.levels)))
  f.change[is.na(f.change)] <- ""
  f.change
}

DefaultAssay(seu) <- "AUCell"
# should be factors
# change = STIM / CTRL
seu$group <- factor(seu$group, levels = c("CTRL", "STIM"))
de.list <- DERegulon(seu, celltype = "celltype", group = "group", test.use = "wilcox")

# avg_diff = STIM - CTRL
View(de.list[[1]])
summary(de.list[[1]]$avg_diff)

ggplot(de.list[[1]], aes(diff_rate, -log10(p_val_adj))) +
  geom_point()

de.regulons <- format_change(mylist = de.list, celltype.levels = levels(seu$celltype))

data.use <- full_join(vd.res, de.regulons, by = "gene")

## check results
VlnPlot(seu, group.by = "celltype", features = "THRB(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F)

FeaturePlot(seu, reduction = "umap", split.by = "group", features = "THRB(+)") &
  scale_color_viridis_c()

FeaturePlot(seu, reduction = "umap", split.by = "group", features = "THRB") &
  scale_color_viridis_c()
```

##### B:

```R


```







