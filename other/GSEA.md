## GSEA基因集富集分析

```properties
## 背景：
基因集富集分析
```

#### 1.fgsea

```R
### https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
rm(list = ls())

library(fgsea)
library(data.table)
library(ggplot2)
data(examplePathways)
data(exampleRanks)
set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways,
                  stats = exampleRanks,
                  minSize = 15,
                  maxSize = 500,
                  nperm = 10000)
psw <- examplePathways[["5991130_Programmed_Cell_Death"]]

plotEnrichment(psw, exampleRanks) +
  labs(title = "Programmed Cell Death")

index <- grep("5991130", fgseaRes$pathway)
anno <- fgseaRes[index, c("NES", "pval", "padj")]
lab <- paste0(names(anno), "=", round(anno, 3), collapse = "\n")

plotEnrichment(psw, exampleRanks) +
  labs(title = "Programmed Cell Death") +
  annotate("text", 0, as.numeric(fgseaRes[index, "ES"]) * .9, label = lab, hjust = 0, vjust = 0)

#### 多个
topPathwaysUp <- fgseaRes[ES > 0,][head(order(pval), n = 10), pathway]
topPathwaysDown <- fgseaRes[ES < 0,][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes,
              gseaParam = 0.5)

### 精简
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                      examplePathways, exampleRanks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes,
              gseaParam = 0.5)

## 自带数据
rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package = "fgsea")
ranks <- read.table(rnk.file, header = TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)

gmt.file <- system.file("extdata", "mouse.reactome.gmt", package = "fgsea")
pathways <- gmtPathways(gmt.file)
fgseaRes <- fgsea(pathways, ranks, minSize = 15, maxSize = 500, nperm = 10000)
```

#### 2.clusterprofiler

```R
rm(list = ls())
library(clusterProfiler)

# 差异分析得到的表达矩阵
gene <- rownames(Diff_data)
gene = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gene <- dplyr::distinct(gene, SYMBOL, .keep_all = TRUE)
gene_df <- data.frame(logFC = allDiff$logFC, SYMBOL = rownames(allDiff))
gene_df <- merge(gene_df, gene, by = "SYMBOL")

### genelist获取
geneList <- gene_df$logFC
names(geneList) = gene_df$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

## 读入hallmarks gene set，broad网站
hallmarks <- read.gmt("h.all.v7.1.entrez.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)

### 作图
cnetplot(y, foldChange = geneList)
y2 <- setReadable(y, "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(y2, showCategory = 4,
         foldChange = geneList,
         colorEdge = T)

### 看整体分布
dotplot(y, showCategory = 12, split = ".sign") + facet_grid(~.sign)

### 选择需要呈现的来作图
yd <- data.frame(y)
library(enrichplot)
gseaplot2(y, "HALLMARK_MYC_TARGETS_V2", color = "red", pvalue_table = T)
gseaplot2(y, 11, color = "red", pvalue_table = T)
ridgeplot(y)
gseaplot2(y, geneSetID = 1:3)

### 自己添加文字加文字
index <- "HALLMARK_MYC_TARGETS_V2"
gseaplot2(y, index, color = "green")

anno <- yd[index, c("enrichmentScore", "NES", "pvalue", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=", round(anno, 3), collapse = "\n")

p <- gseaplot2(y, index, color = "green")
p[[1]] <- p[[1]] + annotate("text", 15000, -0.6, label = lab, hjust = 0, vjust = 0, size = 5)
p
```

#### 3.misgdb

```R


```
