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

#### 2.Seuarat打分

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

#### 5.corr

```R
library(UCell)
library(Seurat)
library(clusterProfiler)
## 读入基因集
genesets <- read.gmt("data/ENCODE_TF_ChIP-seq_2015.txt")
signatures <- split(genesets$gene, genesets$term)
max(sapply(signatures, length))
#signatures <-lapply(signatures, function(i) sample(i,size = length(i)*0.1))
### 数目太大，选择自己感兴趣的基因集
### 来自于前期GSEA分析的结果
signatures <- signatures["STAT2 K562 hg19"]
scobj <- AddModuleScore_UCell(scobj, features = signatures, name = "_UCell")

## UCell 对名称敏感
marker_genes <- "STAT2.K562.hg19_UCell"
FeaturePlot(scobj, features = marker_genes, split.by = "group", order = T)

### 如果有多个基因集要查看
genesets <- read.gmt("data/ENCODE_TF_ChIP-seq_2015.txt")
signatures <- split(genesets$gene, genesets$term)
terms <- c("STAT2 K562 hg19", "STAT1 K562 hg19", "STAT1 HeLa-S3 hg19",
           "IKZF1 GM12878 hg19", "STAT3 HeLa-S3 hg19", "EP300 CH12.LX mm9",
           "BATF GM12878 hg19")
signatures <- signatures[terms]
max(sapply(signatures, length))
names(signatures) <- gsub(" ", "_", names(signatures))
#signatures <-lapply(signatures, function(i) sample(i,size = length(i)*0.1))
scobj <- AddModuleScore_UCell(scobj, features = signatures, name = "_UCell", maxRank = 2000)
test <- scobj@meta.data

## 画图
term <- c("STAT2 K562 hg19")
marker_genes <- gsub(" ", "_", term)
marker_genes <- paste0(marker_genes, "_UCell")
FeaturePlot(scobj, features = marker_genes, split.by = "group", order = T)
```

#### 6.genset

```R
library(Seurat)

ComputeModuleScore <- function(x, gene.sets, min.size = 20, batch.size = 500, cores = 1, ...) {
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


regulons <- clusterProfiler::read.gmt("data/s4_post_scenic.regulons.gmt")
regulons.list <- split(regulons$gene, regulons$term)
names(regulons.list) <- sub("[0-9]+g", "\\+", names(regulons.list))

scobj <- ComputeModuleScore(scobj, gene.sets = regulons.list, min.size = 10, cores = 1)


### 使用新的降维数据进行umap
DefaultAssay(scobj) <- "AUCell"

scobj <- RunUMAP(scobj,
                 features = rownames(scobj),
                 metric = "correlation",
                 reduction.name = "umap_ras",
                 reduction.key = "umap_ras")

p1 <- DimPlot(scobj, reduction = "umap", group.by = "celltype") + ggsci::scale_color_d3("category20")
p2 <- DimPlot(scobj, reduction = "umap", group.by = "group")
p1 + p2

## ras-umap
p3 <- DimPlot(scobj, reduction = "umap_ras", group.by = "celltype") + ggsci::scale_color_d3("category20")
p4 <- DimPlot(scobj, reduction = "umap_ras", group.by = "group")

p3 + p4

library(patchwork)
(p1 + p2) / (p3 + p4)

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


vd.vars <- c("celltype", "group")
meta.data <- scobj@meta.data[, vd.vars]
ras.data = t(scobj@assays$AUCell@data)

vd.res <- VarDecompose(data = ras.data,
                       meta.data = meta.data,
                       vd.vars = vd.vars, cores = 1)

## 作图检测
DefaultAssay(scobj) <- "AUCell"
FeaturePlot(scobj, features = "IRF9(+)", split.by = "group", reduction = "umap")
FeaturePlot(scobj, features = "STAT1(+)", split.by = "group", reduction = "umap")
FeaturePlot(scobj, features = "KLF4(+)", split.by = "group", reduction = "umap")
FeaturePlot(scobj, features = "IRF7(+)", split.by = "group", reduction = "umap")
FeaturePlot(scobj, features = "BACH1(+)", reduction = "umap")
FeaturePlot(scobj, features = "IRF7(+)", reduction = "umap")


## 探索
library(dplyr)
dput(vd.res %>%
       arrange(desc(group)) %>%
       pull(gene) %>%
       head(60))

module <- c("IRF9(+)", "ELF1(+)", "IRF1(+)", "STAT1(+)", "ETV7(+)", "NFE2L3(+)",
            "STAT2(+)", "IRF2(+)", "IRF7(+)", "IRF8(+)", "BATF3(+)", "NFYC(+)",
            "MAX(+)", "ELK3(+)", "ZNF580(+)", "JUND(+)", "ZSCAN16(+)", "CEBPG(+)",
            "JUN(+)", "HOXB2(+)", "THAP11(+)", "TCF7L2(+)", "NR3C1(+)", "CEBPD(+)",
            "ZNF76(+)", "CREB1(+)", "TAF1(+)", "HDX(+)", "ETV6(+)", "FOSB(+)",
            "ETS1(+)", "ZFHX3(+)", "ATF6B(+)", "PHF20(+)", "TFEC(+)", "ZBTB25(+)",
            "SREBF2(+)", "ZBTB7A(+)", "BCLAF1(+)", "ATF4(+)", "TFEB(+)",
            "CLOCK(+)", "FLI1(+)", "ATF5(+)", "FOSL2(+)", "YY1(+)", "ATF3(+)",
            "POU2F2(+)", "KLF13(+)", "STAT6(+)", "TFE3(+)", "NFE2L2(+)",
            "HMGA1(+)", "ELF2(+)", "ING4(+)", "RUNX1(+)", "NFIC(+)", "ELK1(+)",
            "MAFB(+)", "NFKB2(+)")


scobj <- RunUMAP(scobj,
                 features = setdiff(rownames(scobj), module),
                 metric = "correlation",
                 reduction.name = "umap_ras2",
                 reduction.key = "umap_ras2")

p5 <- DimPlot(scobj, reduction = "umap_ras2", group.by = "celltype") + ggsci::scale_color_d3("category20")
p6 <- DimPlot(scobj, reduction = "umap_ras2", group.by = "group")

p5 + p6


### one more time

dput(vd.res %>%
       arrange(desc(celltype)) %>%
       pull(gene) %>%
       head(80))

module <- c("KLF4(+)", "CEBPB(+)", "RXRA(+)", "BACH1(+)", "JUNB(+)", "VDR(+)",
            "SPI1(+)", "FOS(+)", "ETS2(+)", "NFKB1(+)", "MAFB(+)", "IRF5(+)",
            "ATF5(+)", "CREM(+)", "ATF3(+)", "MITF(+)", "NFE2L2(+)", "SPIB(+)",
            "TCF4(+)", "TFEC(+)", "CEBPD(+)", "ZNF467(+)", "CUX1(+)", "ETV6(+)",
            "GFI1(+)", "ETS1(+)", "TCF7L2(+)", "MYC(+)", "HMGA1(+)", "EOMES(+)",
            "REL(+)", "BCLAF1(+)", "TBX21(+)", "HES1(+)", "NFATC2(+)", "NR3C1(+)",
            "ATF4(+)", "ATF1(+)", "GATA3(+)", "LEF1(+)", "GABPB1(+)", "ELF4(+)",
            "YY1(+)", "FLI1(+)", "RFX5(+)", "IKZF1(+)", "HOXB2(+)", "PPARG(+)",
            "NFKB2(+)", "KLF2(+)", "RUNX3(+)", "ATF2(+)", "RELB(+)", "THAP11(+)",
            "ARID3A(+)", "JUN(+)", "ZBTB25(+)", "TFE3(+)", "ELK4(+)", "ELK3(+)",
            "ELF2(+)", "XBP1(+)", "KLF13(+)", "TFEB(+)", "HOXB4(+)", "FOSL2(+)",
            "MAFF(+)", "CEBPG(+)", "THAP1(+)", "RELA(+)", "IRF8(+)", "ETV3(+)",
            "IRF7(+)", "FOSB(+)", "ZBTB7A(+)", "PHF20(+)", "BATF3(+)", "MAX(+)",
            "ZBTB7B(+)", "TFAP4(+)")


scobj <- RunUMAP(scobj,
                 features = setdiff(rownames(scobj), module),
                 metric = "correlation",
                 reduction.name = "umap_ras3",
                 reduction.key = "umap_ras3")

p7 <- DimPlot(scobj, reduction = "umap_ras3", group.by = "celltype") + ggsci::scale_color_d3("category20")
p8 <- DimPlot(scobj, reduction = "umap_ras3", group.by = "group")

p7 + p8


library(patchwork)
(p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8)
```
