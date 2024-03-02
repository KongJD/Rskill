load(file = "data/EGG.Rdata")
### B加载R包
library(ggplot2)
library(clusterProfiler)
### C画图
barplot(EGG)
dotplot(EGG)


load(file = "data/go.Rdata")
### 画直方图
barplot(go, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
barplot(go, label_format = 60, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
### 画点图
dotplot(go, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")
dotplot(go, label_format = 60, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free")


if (T) {
  x = EGG
  ## 提取富集分析的数据
  dd = x@result
  ## 计算富集分数
  dd$richFactor = dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))
  ## 提取p值小于0.05 的数据
  dd <- dd[dd$p.adjust < 0.05,]

  library(ggplot2)
  ## 正式画图
  ggplot(dd, aes(richFactor, forcats::fct_reorder(Description, richFactor))) +
    ## 画横线
    geom_segment(aes(xend = 0, yend = Description)) +
    ## 画点
    geom_point(aes(color = p.adjust, size = Count)) +
    ## 调整颜色的区间,begin越大，整体颜色越明艳
    scale_color_viridis_c(begin = 0.3, end = 1) +
    ## 调整泡泡的大小
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    xlab("Rich factor") +
    ylab(NULL) +
    ggtitle("")
}

################################################

rm(list = ls())
## 判断癌和癌旁组织
## 1.加载数据
load(file = "data/TCGA_exprSet_plot.Rdata")
metadata <- data.frame(TCGA_id = rownames(exprSet))

## for循环批量处理
for (i in 1:nrow(metadata)) {
  print(i)
  num <- as.numeric(substring(metadata[i, 1], 14, 15))
  if (num %in% seq(1, 9)) { metadata[i, 2] = "Tumor" }
  if (num %in% seq(10, 29)) { metadata[i, 2] = "Normal" }
}
## 修改列名
colnames(metadata) <- c("TCGA_id", "sample")

#################################################

load(file = "data/plotdf.Rdata")
if (T) {
  library(ggplot2)
  options(scipen = 2)
  ggplot(plotdf, aes(-log10(p.value), cor)) +
    geom_point(size = 6, fill = "purple", shape = 21, colour = "black", stroke = 2) +
    scale_y_continuous(expand = c(0, 0), limits = c(-1.1, 1.1), breaks = seq(-1, 1, 0.2)) +
    scale_x_log10(limits = c(0.01, 1000), breaks = c(0.01, 0.1, 10, 1000)) +
    geom_hline(yintercept = 0, size = 1.5) +
    geom_vline(xintercept = -log10(0.05), size = 1.5) +
    labs(x = bquote(-log[10] ~ italic("P")), y = "Pearson correlation (r)") +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(face = "bold", size = 16),
          axis.ticks.length = unit(.4, "cm"),
          axis.ticks = element_line(colour = "black", size = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
}

pancorplot <- function(data = plotdf, anotate = "none") {
  require(ggplot2)
  require(ggrepel)
  options(scipen = 2)
  p = ggplot(plotdf, aes(-log10(p.value), cor)) +
    geom_point(data = subset(plotdf, plotdf$p.value >= 0.05), size = 6, fill = "grey", alpha = 0.6, shape = 21, colour = "black", stroke = 1.5) +
    geom_point(data = subset(plotdf, plotdf$p.value < 0.05 & plotdf$cor >= 0), size = 6, fill = "red", alpha = 0.6, shape = 21, colour = "black", stroke = 1.5) +
    geom_point(data = subset(plotdf, plotdf$p.value < 0.05 & plotdf$cor < 0), size = 6, fill = "blue", alpha = 0.6, shape = 21, colour = "black", stroke = 1.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(-1.1, 1.1), breaks = seq(-1, 1, 0.2)) +
    scale_x_log10(limits = c(0.01, 1000), breaks = c(0.01, 0.1, 10, 1000)) +
    geom_hline(yintercept = 0, size = 1.5) +
    geom_vline(xintercept = -log10(0.05), size = 1.5) +
    labs(x = bquote(-log[10] ~ italic("P")), y = "Pearson correlation (r)") +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(face = "bold", size = 16),
          axis.ticks.length = unit(.4, "cm"),
          axis.ticks = element_line(colour = "black", size = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  index = anotate == "none"
  index = index[1]
  if (index) {
    p
  }else {
    p +
      geom_point(data = subset(plotdf, plotdf$type %in% anotate),
                 size = 6, shape = 21, colour = "green", stroke = 2) +
      geom_label_repel(data = subset(plotdf, plotdf$type %in% anotate),
                       aes(label = type), col = "black", size = 6,
                       box.padding = 10,
                       arrow = arrow(angle = 30, length = unit(0.25, "inches"),
                                     ends = "first", type = "closed"),
                       segment.size = 1,
                       segment.color = "red")
  }
}


pancorplot(plotdf)
pancorplot(plotdf, anotate = c("BRCA"))
pancorplot(plotdf, anotate = c("BRCA", "GBM"))
