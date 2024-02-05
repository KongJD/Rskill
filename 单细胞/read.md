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
## 2.

```