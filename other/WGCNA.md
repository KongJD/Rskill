## WGCNA分析

#### 背景

```properties
用来描述不同样品之间基因关联模式的系统生物学方法，可以用来鉴定高度协同变化的基因集
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
library(WGCNA)
### 设置多线程
enableWGCNAThreads()
### 软阈值的确定，
### 使用函数pickSoftThreshold
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
length(powers)
### 这个power的长度，决定了最终返回结果的行数目

# Call the network topology analysis function
### 模块这里都是作用于表达量数据，没有性状数据的事情
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

test <- sft$fitIndices

sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"));
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

### 根据图，选择6
## 一步法网络构建以及模块发现
### power = 6，来自于上一步
##  maxBlockSize = 5000, 根据自己的电脑可以调整，16GB内存20000，32GB内存30000
### TOMType 
### 这一步是所有分析的基石
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
### 返回的net是个列表，Module Epigene 表达数据也在里面
table(net$colors)
### 总共找到18个模块
### 0表示的是不在模块内的基因
# open a graphics window
sizeGrWindow(12, 9)
### 标签转为颜色
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plot(net$dendrograms[[1]])
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
```

