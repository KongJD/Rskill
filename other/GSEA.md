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
library(msigdbr)
msigdbr_species()
m_df = msigdbr(species = "Homo sapiens")
library(dplyr)
m_df %>%
  distinct(gs_cat, gs_subcat) %>%
  arrange(gs_cat, gs_subcat)

### 人
h_df = msigdbr(species = "Homo sapiens", category = "H")

### 鼠
m_df = msigdbr(species = "Mus musculus", category = "H")

### 选择的方式
m_df1 = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
dd <- msigdbr(species = "Mus musculus")
colnames(dd)

###选取鼠的C2，CGP，entrizID
library(dplyr)
m_df2 <- dd %>%
  filter(species_name == "Mus musculus") %>%
  filter(gs_cat == "C2", gs_subcat == "CGP") %>%
  select(gs_name, entrez_gene)

### 选人，hallmarks，两列
dd <- msigdbr(species = "Homo sapiens")
h_df2 <- dd %>%
  filter(gs_cat == "H") %>%
  select(gs_name, entrez_gene)

### 比较
h_df3 <- read.gmt("h.all.v7.1.entrez.gmt")
m_list = split(h_df2$entrez_gene, h_df2$gs_name)

### fgsea的输入
pathways <- gmtPathways("h.all.v7.1.entrez.gmt")
### https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
```

#### 4.experdata

```R
y <- as.numeric(exprSet[1,])
rownames <- rownames(exprSet)
cor_data_df <- data.frame(rownames)
### 批量相关性分析
for (i in 1:length(rownames)) {
  test <- cor.test(as.numeric(exprSet[i,]), y, method = "spearman")
  cor_data_df[i, 2] <- test$estimate
  cor_data_df[i, 3] <- test$p.value
}
names(cor_data_df) <- c("symbol", "correlation", "pvalue")

## geneList 
geneList <- cor_data_df$correlation
names(geneList) = cor_data_df$symbol
geneList = sort(geneList, decreasing = TRUE)


### 基因集
library(msigdbr)
dd <- msigdbr(species = "Homo sapiens")
hallmarks <- dd %>%
  filter(gs_cat == "H") %>%
  select(gs_name, gene_symbol)

### GSEA分析
library(clusterProfiler)
y <- GSEA(geneList, TERM2GENE = hallmarks)

### 看整体分布
dotplot(y, showCategory = 12, split = ".sign") + facet_grid(~.sign)
```

#### 5.gsva+limma

```R
library(GSVA)
kegggmt <- read.gmt("c2.cp.kegg.v7.1.symbols.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
## 下面有个报错错误: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
## 把包降级 remotes::install_version("matrixStats", version="1.1.0")
kegg1 <- gsva(expr = as.matrix(exprSet), kegg_list, kcdf = "Gaussian", method = "gsva", parallel.sz = 1)

keggSet <- getGmt("c2.cp.kegg.v7.1.symbols.gmt")
kegg2 <- gsva(expr = as.matrix(exprSet), keggSet, kcdf = "Gaussian", method = "gsva", parallel.sz = 1)

library(pheatmap)
##用名称提取部分数据用作热图绘制
group <- c(rep("con", 3), rep("treat", 3))
annotation_col <- data.frame(group)
rownames(annotation_col) <- colnames(kegg2)

pheatmap(kegg2,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         annotation_legend = TRUE,
         show_rownames = F,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cellwidth = 80, cellheight = 2,
         fontsize = 10)


library(limma)
exprSet <- kegg2
### 差异分析

group <- c(rep("con", 3), rep("treat", 3))
## 分组变成向量，并且限定leves的顺序
## levels里面，把对照组放在前面
group <- factor(group, levels = c("con", "treat"), ordered = F)
## 1.构建比较矩阵
design <- model.matrix(~group)
## 比较矩阵命名
colnames(design) <- levels(group)

##2.线性模型拟合
fit <- lmFit(exprSet, design)
##3.贝叶斯检验
fit2 <- eBayes(fit)
#4.输出差异分析结果,其中coef的数目不能操过design的列数
# 此处的2代表的是第二列和第一列的比较
allDiff = topTable(fit2, adjust = 'fdr', coef = 2, number = Inf) 
```

#### 6.ssgsva

```R
### 每个细胞对应表达的mark基因（类似通路对应基因）
expr <- data.table::fread("exprMat.txt", data.table = F)
rownames(expr) <- expr[, 1]
expr <- expr[, -1]
expr <- as.matrix(expr)

### 3.使用ssGSEA量化免疫浸润
library(GSVA)
gsva_data <- gsva(expr, cellMarker, method = "ssgsea")
```


