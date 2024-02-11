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


```
