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

#### 2 网络构建 

```R
library(WGCNA)
### 设置多线程
enableWGCNAThreads()
### 软阈值的确定，
### 使用函数pickSoftThreshold
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

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

### 总共找到18个模块
### 0表示的是不在模块内的基因
# open a graphics window
sizeGrWindow(12, 9)
### 标签转为颜色
mergedColors = labels2colors(net$colors)
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

#### 3.模块与性状相连

```R
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
test <- moduleEigengenes(datExpr, moduleColors)
### 提取结果
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
### 重新排序，改变的是列的位置 ，没有实际意义，主要用于绘图
MEs = orderMEs(MEs0)
#############################################################################
### 很重要的步骤，很简单的操作
### 求相关性，datTraits 
moduleTraitCor = cor(MEs, datTraits, use = "p")
### 算p值
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10, 6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g)
names(weight) = "weight"
# names (colors) of the modules,cong
## 从第三位开始选取
modNames = substring(names(MEs), 3)

### 计算基因的相关性
### module membership MM
### Gene Significance GS
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

### 单独一个模块
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(weight), sep = "")
names(GSPvalue) = paste("p.GS.", names(weight), sep = "");

# 举例子：brown 和weight
module = "brown"
column = match(module, modNames)
table(moduleColors)
### 提取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1, 1))
### 本质上就是一个散点图
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])

verboseScatterplot(MM, GS,
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


names(datExpr)[moduleColors == "brown"]

### 模块注释
annot = data.table::fread(file = "GeneAnnotation.csv", data.table = F)

probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
### 排序
modOrder = order(-abs(cor(MEs, weight, use = "p")))
# Add module membership information in the chosen order
### 按列操作
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""), paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
geneInfo = geneInfo0[geneOrder,]
```

