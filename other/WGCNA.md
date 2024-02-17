## WGCNA分析

#### 背景

```properties
## https://omicverse.readthedocs.io/en/latest/Tutorials-bulk/t_wgcna/
```

#### 1.清洗数据

```R
### Network analysis of liver expression data from female mice: finding modules related to body weight
### 分析雌性鼠肝脏数据，找出跟体重相关的模块
rm(list = ls())
library(WGCNA)
femData = data.table::fread("LiverFemale3600.csv", data.table = F)
### 提取正确的表达量数据
### 样本在行，基因在列，本次前八个信息不要
datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
colnames(datExpr0) = femData$substanceBXH
rownames(datExpr0) = names(femData)[-c(1:8)]

### 第一是样本，在行
### 第二是基因，在列
### 达标的标准是什么啊？goodSamplesGenes
gsg = goodSamplesGenes(datExpr0, verbose = 3)

### 如果没有达标就需要筛选
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

### 聚类看异常样本
dist(mtcars[1:5, 1:4])
dist(mtcars[1:5, 1:5])
sampleTree = hclust(dist(datExpr0), method = "average")
### 画图
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

### 去除异常样本
abline(h = 15, col = "red")
### 使用WGCNA中的cutreeStatic函数来切割
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)

# clust 1 contains the samples we want to keep.
keepSamples = (clust == 1)
### 选取样本
datExpr = datExpr0[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

### 读取表型信息
traitData = data.table::fread("ClinicalTraits.csv", data.table = F)
### 选取自己需要的
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36)]
femaleSamples = rownames(datExpr)
### 获取第一个数据在第二个数据中的位置，返回的是位置
traitRows = match(femaleSamples, allTraits$Mice)

### 提取交集
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

### 重新聚类
sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2)

# Convert traits to a color representation: white means low, red means high, grey means missing entry
## 颜色,白色是low，红色是高，灰色是缺失值
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
## traitColors该部分可以是单独一列，也可以用cbind传入多列，此处传入的是矩阵,26列
## Each column will be plotted as a horizontal row of colors under the dendrogram.
plot(sampleTree2)
plotDendroAndColors(sampleTree2,
                    traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
```

#### 2

```R


```

