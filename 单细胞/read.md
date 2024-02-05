## 单细胞

```properties
## 背景
https://satijalab.org/seurat/articles/get_started_v5_new
```

#### 基本流程

```properties
## 1. 上游分析：略
见上述的转录组的流程到 featurecounts 得到 原始的表达矩阵
## 数据直接下载
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111229
```

```R
## 2. 表达矩阵
a = read.table('./GSE111229_Mammary_Tumor_fibroblasts_768samples_rawCounts.txt.gz',
               header = T, sep = '\t')
dat <- a[apply(a, 1, function(x) sum(x > 1) > floor(ncol(a) / 50)),]
dat = log2(edgeR::cpm(dat) + 1)
##归一化的一种选择，这里是CPM(count-per-million，每百万碱基中每个转录本的count值)
###CPM只对read count相对总reads数做了数量的均一化，去除文库大小差异。

hc = hclust(dist(t(dat)))

clus = cutree(hc, 4) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
group_list = as.factor(clus)
table(group_list)

#提取批次信息
library(stringr)
plate = str_split(colnames(dat), '_', simplify = T)[, 3]
table(plate)


n_g = apply(a, 2, function(x) sum(x > 1))
df = data.frame(g = group_list, plate = plate, n_g = n_g)
df$all = 'all'
metadata = df


## rpkm
a = read.table('./GSE111229_Mammary_Tumor_fibroblasts_768samples_rpkmNormalized.txt.gz',
               header = T, sep = '\t')
dat = a[apply(a, 1, function(x) sum(x > 0) > floor(ncol(a) / 50)),] #筛选表达量合格的行,列数不变
#层次聚类
hc = hclust(dist(t(log(dat + 0.1)))) ##样本间层次聚类
# 如果是基因聚类，可以选择 wgcna 等算法
plot(hc, labels = F)
clus = cutree(hc, 4) #对hclust()函数的聚类结果进行剪枝，即选择输出指定类别数的系谱聚类结果。
group_list = as.factor(clus)
table(group_list)

#提取批次信息
colnames(dat)
library(stringr)
plate = str_split(colnames(dat), '_', simplify = T)[, 3]
table(plate)

n_g = apply(a, 2, function(x) sum(x > 0)) #统计每个样本有表达的有多少行（基因）
df = data.frame(g = group_list, plate = plate, n_g = n_g)
##(样本为行名，列分别为：样本分类信息，样本分组，样本表达的基因数【注意：不是表达量的和，而是种类数或者说个数】)
df$all = 'all'
metadata = df
```

```R
## 2.表达矩阵画图探索
group_list = metadata$g
table(group_list)
cg = names(tail(sort(apply(dat, 1, sd)), 100)) ##取表达量标准差最大的100行的行名

library(pheatmap)
mat = log2(dat[cg,] + 0.01)
pheatmap(mat, show_colnames = F, show_rownames = F,)
n = t(scale(t(mat)))
n[n > 2] = 2
n[n < -2] = -2
ac = data.frame(g = group_list)
rownames(ac) = colnames(n)
pheatmap(n, show_colnames = F, show_rownames = F,
         annotation_col = ac,
)

## pca图
dat_back = dat
dat = dat_back
dat = t(dat)
dat = as.data.frame(dat)
dat = cbind(dat, group_list)
dat[, ncol(dat)]
library("FactoMineR")
library("factoextra")
dat.pca <- PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca, repel = T,
             geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
             col.ind = dat$group_list, # color by groups 颜色组
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses 集中成椭圆
             legend.title = "Groups"
)

## 
exprSet = dat
mean_per_gene <- apply(exprSet, 1, mean, na.rm = TRUE) #对表达矩阵每行求均值
sd_per_gene <- apply(exprSet, 1, sd, na.rm = TRUE) #对表达矩阵每行求标准差
mad_per_gene <- apply(exprSet, 1, mad, na.rm = TRUE) #对表达矩阵每行求绝对中位差

cv_per_gene <- data.frame(mean = mean_per_gene,
                          sd = sd_per_gene,
                          mad = mad_per_gene,
                          cv = sd_per_gene / mean_per_gene)
rownames(cv_per_gene) <- rownames(exprSet)
with(cv_per_gene, plot(log10(mean), log10(cv)))
with(cv_per_gene, plot(log10(mean), log10(cv^2)))
cv_per_gene$log10cv2 = log10(cv_per_gene$cv^2)
cv_per_gene$log10mean = log10(cv_per_gene$mean)
library(ggpubr)
cv_per_gene = cv_per_gene[cv_per_gene$log10mean < 4,]
cv_per_gene = cv_per_gene[cv_per_gene$log10mean > 0,]
ggscatter(cv_per_gene, x = 'log10mean', y = 'log10cv2',
          color = "black", shape = 16, size = 1, # Points color, shape and size
          xlab = 'log10(mean)RPKM', ylab = "log10(cv^2)",
          add = "loess", #添加拟合曲线
          add.params = list(color = "red", fill = "lightgray"),
          cor.coeff.args = list(method = "spearman"),
          label.x = 3, label.sep = "\n",
          dot.size = 2,
          ylim = c(-0.5, 3),
          xlim = c(0, 4)
)
```

```R
## 3.批次效应
dat_back = dat
dat = dat_back
dat = t(dat)
dat = as.data.frame(dat)
dat = cbind(dat, plate)
table(dat$plate)
library("FactoMineR")
library("factoextra")

dat.pca <- PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca, #repel =T,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$plate, # color by groups
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_cells_PCA_by_plate.png')

library(Rtsne)
dat_matrix <- t(dat_back)


dat_matrix = log2(dat_matrix + 0.01)
set.seed(42)
tsne_out <- Rtsne(dat_matrix, pca = FALSE, perplexity = 30, theta = 0.0) # Run TSNE
plot(tsne_out$Y, col = plate)
# https://distill.pub/nc16/misread-tsne/


library(ggpubr)
df = as.data.frame(tsne_out$Y)
colnames(df) = c("X", 'Y')
df$plate = plate
df$g = metadata$g
ggscatter(df, x = "X", y = "Y", color = "g"
          # palette = c("#00AFBB", "#E7B800" ) 
)
```

```R
## 4.基因数量
library(ggpubr)
ggviolin(df, x = "all", y = "n_g", fill = "all",
         add = "boxplot", add.params = list(fill = "white"))
ggviolin(df, x = "plate", y = "n_g", fill = "plate",
         #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))
ggviolin(df, x = "g", y = "n_g", fill = "g",
         add = "boxplot", add.params = list(fill = "white")) + stat_compare_means()
```