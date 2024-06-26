## GEO数据

```properties
背景：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872
```

#### 基本流程：

```R
## 1. 下载数据
library(GEOquery)
f = 'GSE42872_eSet.Rdata'
if (!file.exists(f)) {
  gset <- getGEO('GSE42872', destdir = ".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset, file = f)
}
a <- read.table("./GSE42872_series_matrix.txt.gz", header = T, sep = "\t", fill = T, quote = "", comment.char = "!")
```

```R
## 2. 探针id转换
pd <- gset$GSE42872_series_matrix.txt.gz@phenoData@data
# dat1 <- exprs(gset[[1]])
dat <- exprs(gset[[1]])

library(hugene10sttranscriptcluster.db)
ids <- toTable(hugene10sttranscriptclusterSYMBOL)

dat_final <- dat[rownames(dat) %in% ids$probe_id,]
ids <- ids[match(rownames(dat_final), ids$probe_id),]

# 探针对应的基因一样时选平均数最大的
tmp <- by(dat_final, ids$symbol, function(x) {
  rownames(x)[which.max(rowMeans(x))]
})
dat_final_res <- dat_final[tmp,]
# x = dat_final[ids[,2]=="IGKC",]
# dat_final[rownames(x)[which.max(rowMeans(x))],]
```

```R
## 3. 差异分析
group_list <- unlist(lapply(as.character(pd$title), function(x) {
  strsplit(x, " ")[[1]][4] }
))

library(limma)

design <- model.matrix(~0 + factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(dat_final_res)

contrast.matrix <- makeContrasts(paste0(unique(group_list), collapse = "-"),
                                 levels = design)

deg <- function(exprSet, design, contrast.matrix) {
  fit <- lmFit(exprSet, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef = 1, n = Inf)
  nrDEG = na.omit(tempOutput)
  return(nrDEG)
}

res <- deg(dat_final_res, design, contrast.matrix)

## volcano 
library(ggpubr)
plot(res$logFC, -log10(res$P.Value))
res$v = -log10(res$P.Value)
ggscatter(res, x = "logFC", y = "v", size = 0.5)
res$g = ifelse(res$P.Value > 0.01, 'stable',
               ifelse(res$logFC > 2, 'up',
                      ifelse(res$logFC < -2, 'down', 'stable')) #接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
res$name = rownames(res)
ggscatter(res, x = "logFC", y = "v", size = 0.5, color = 'g')
ggscatter(res, x = "logFC", y = "v", color = "g", size = 0.5,
          label = "name", repel = T,
          #label.select = rownames(df)[df$g != 'stable'] ,
          label.select = head(rownames(res)), #挑选一些基因在图中显示出来
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))

ggscatter(res, x = "AveExpr", y = "logFC", size = 0.2)
res$p_c = ifelse(res$P.Value < 0.001, 'p<0.001',
                 ifelse(res$P.Value < 0.01, '0.001<p<0.01', 'p>0.01'))
ggscatter(res, x = "AveExpr", y = "logFC", color = "p_c", size = 0.2,
          palette = c("green", "red", "black"))


### pheatmap
library(pheatmap)

x = res$logFC
names(x) = rownames(res)
cg = c(names(head(sort(x), 100)),
       names(tail(sort(x), 100)))
pheatmap(dat[cg,], show_colnames = F, show_rownames = F)
n = t(scale(t(dat[cg,]))) #通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
pheatmap(n, show_colnames = F, show_rownames = F)
n[n > 2] = 2
n[n < -2] = -2
pheatmap(n, show_colnames = F, show_rownames = F)
ac = data.frame(g = group_list)
rownames(ac) = colnames(n)
pheatmap(n, show_colnames = F,
         show_rownames = F,
         cluster_cols = F,
         annotation_col = ac,) #列名注释信息为ac即分组信息
```

```R
## 4.结果注释（富集分析）
ibrary(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

res$probe_id <- rownames(res)
res <- merge(res, ids, by = "probe_id")


df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db)

final = merge(res, df, by.y = 'SYMBOL', by.x = 'symbol')

gene_up = final[final$g == 'UP', 'ENTREZID']
gene_down = final[final$g == 'DOWN', 'ENTREZID']
gene_diff = c(gene_up, gene_down)
gene_all = as.character(final[, 'ENTREZID'])

## GO kegg database analysis
if (T) {
  ###   over-representation test
  install.packages("R.utils")
  R.utils::setOption("clusterProfiler.download.method", 'auto')
  ## 这里因为clusterProfiler需要 r 4.2 以上，后面的 kegg 和 go的暂时没跑出来 ，后面创建个新环境跑
  kk.up <- enrichKEGG(gene = gene_up,
                      organism = 'hsa',
                      universe = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  head(kk.up)[, 1:6]
  dotplot(kk.up); ggsave('kk.up.dotplot.png')
  kk.down <- enrichKEGG(gene = gene_down,
                        organism = 'hsa',
                        universe = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff = 0.9)
  head(kk.down)[, 1:6]
  dotplot(kk.down); ggsave('kk.down.dotplot.png')
  kk.diff <- enrichKEGG(gene = gene_diff,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[, 1:6]
  dotplot(kk.diff); ggsave('kk.diff.dotplot.png')

  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg <- kegg_down_dt[kegg_down_dt$pvalue < 0.05,]; down_kegg$group = -1
  up_kegg <- kegg_up_dt[kegg_up_dt$pvalue < 0.05,]; up_kegg$group = 1

  kegg_plot <- function(up_kegg, down_kegg) {
    dat = rbind(up_kegg, down_kegg)
    colnames(dat)
    dat$pvalue = -log10(dat$pvalue)
    dat$pvalue = dat$pvalue * dat$group

    dat = dat[order(dat$pvalue, decreasing = F),]

    g_kegg <- ggplot(dat, aes(x = reorder(Description, order(pvalue, decreasing = F)), y = pvalue, fill = group)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "blue", high = "red", guide = FALSE) +
      scale_x_discrete(name = "Pathway names") +
      scale_y_continuous(name = "log10P-value") +
      coord_flip() +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("Pathway Enrichment")
  }

  g_kegg = kegg_plot(up_kegg, down_kegg)


  ggsave(g_kegg, filename = 'kegg_up_down.png')

  ###  GSEA 
  data(geneList, package = "DOSE")
  geneList = final$logFC
  names(geneList) = final$ENTREZID
  geneList = sort(geneList, decreasing = T)

  kk_gse <- gseKEGG(geneList = geneList,
                    organism = 'hsa',
                    nPerm = 1000,
                    minGSSize = 120,
                    pvalueCutoff = 0.9,
                    verbose = FALSE)
  head(kk_gse)[, 1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))

  down_kegg <- kk_gse[kk_gse$pvalue < 0.05 & kk_gse$enrichmentScore < 0,]; down_kegg$group = -1
  up_kegg <- kk_gse[kk_gse$pvalue < 0.05 & kk_gse$enrichmentScore > 0,]; up_kegg$group = 1

  g_kegg = kegg_plot(up_kegg, down_kegg)
  print(g_kegg)
  ggsave(g_kegg, filename = 'kegg_up_down_gsea.png')


}

{

  g_list = list(gene_up = gene_up,
                gene_down = gene_down,
                gene_diff = gene_diff)

  if (F) {
    go_enrich_results <- lapply(g_list, function(gene) {
      lapply(c('BP', 'MF', 'CC'), function(ont) {
        cat(paste('Now process ', ont))
        ego <- enrichGO(gene = gene,
                        universe = gene_all,
                        OrgDb = org.Hs.eg.db,
                        ont = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.99,
                        qvalueCutoff = 0.99,
                        readable = TRUE)

        return(ego)
      })
    })
  }

  n1 = c('gene_up', 'gene_down', 'gene_diff')
  n2 = c('BP', 'MF', 'CC')
  for (i in 1:3) {
    for (j in 1:3) {
      fn = paste0('dotplot_', n1[i], '_', n2[j], '.png')
      cat(paste0(fn, '\n'))
      png(fn, res = 150, width = 1080)
      print(dotplot(go_enrich_results[[i]][[j]]))
      dev.off()
    }
  }

}
```

```R
## 5.生存分析 
# 用my.surv <- surv(OS_MONTHS,OS_STATUS=='DECEASED')构建生存曲线。
# 用kmfit2 <- survfit(my.surv~TUMOR_STAGE_2009)来做某一个因子的KM生存曲线。
# 用 survdiff(my.surv~type, data=dat)来看看这个因子的不同水平是否有显著差异，其中默认用是的logrank test 方法。
# 用coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung) 来检测自己感兴趣的因子是否受其它因子(age,gender等等)的影响。
load(file = './GSE11121_survival/step1-output.Rdata')
exprSet = dat
dim(exprSet)
colnames(phe) = c('event', 'grade', 'node', 'size', 'time')
phe = as.data.frame(apply(phe, 2, as.numeric))
boxplot(phe$size)
phe$size = ifelse(phe$size > median(phe$size), 'big', 'small')

library(survival)
library(survminer)
# 利用ggsurvplot快速绘制漂亮的生存曲线图
sfit <- survfit(Surv(time, event) ~ size, data = phe)
summary(sfit)
ggsurvplot(sfit, conf.int = F, pval = TRUE)
ggsurvplot(sfit, palette = c("#E7B800", "#2E9FDF"),
           risk.table = TRUE, pval = TRUE,
           conf.int = TRUE, xlab = "Time in months",
           ggtheme = theme_light(),
           ncensor.plot = TRUE)
## 多个 ggsurvplots作图生存曲线代码合并代码公布
sfit1 = survfit(Surv(time, event) ~ size, data = phe)
sfit2 = survfit(Surv(time, event) ~ grade, data = phe)
splots <- list()
splots[[1]] <- ggsurvplot(sfit1, pval = TRUE, data = phe, risk.table = TRUE)
splots[[2]] <- ggsurvplot(sfit2, pval = TRUE, data = phe, risk.table = TRUE)
# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 1, risk.table.height = 0.4)
# 可以看到grade跟生存显著相关，而size跟病人生存的关系并不显著。

## 挑选感兴趣的基因做差异分析
phe$CBX4 = ifelse(exprSet['CBX4',] > median(exprSet['CBX4',]), 'high', 'low')
table(phe$CBX4)
ggsurvplot(survfit(Surv(time, event) ~ CBX4, data = phe), conf.int = F, pval = TRUE)
phe$CBX6 = ifelse(exprSet['CBX6',] > median(exprSet['CBX6',]), 'high', 'low')
table(phe$CBX6)
ggsurvplot(survfit(Surv(time, event) ~ CBX6, data = phe), conf.int = F, pval = TRUE)
## 批量生存分析 使用  logrank test 方法
mySurv = with(phe, Surv(time, event))
log_rank_p <- apply(exprSet, 1, function(gene) {
  #gene=exprSet[1,]
  phe$group = ifelse(gene > median(gene), 'high', 'low')
  data.survdiff = survdiff(mySurv ~ group, data = phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p = sort(log_rank_p)
head(log_rank_p)
boxplot(log_rank_p)
phe$H6PD = ifelse(exprSet['H6PD',] > median(exprSet['H6PD',]), 'high', 'low')
table(phe$H6PD)
ggsurvplot(survfit(Surv(time, event) ~ H6PD, data = phe), conf.int = F, pval = TRUE)

# 前面如果我们使用了WGCNA找到了跟grade相关的基因模块，然后在这里就可以跟生存分析的显著性基因做交集
# 这样就可以得到既能跟grade相关，又有临床预后意义的基因啦。

## 批量生存分析 使用 coxph 回归方法
mySurv = with(phe, Surv(time, event))
cox_results <- apply(exprSet, 1, function(gene) {
  group = ifelse(gene > median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group, grade = phe$grade, size = phe$size, stringsAsFactors = F)
  m = coxph(mySurv ~ grade + size + group, data = survival_dat)

  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se

  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta / se, p = 1 - pchisq((beta / se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1) / HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])

})
cox_results = t(cox_results)
table(cox_results[, 4] < 0.05)
cox_results[cox_results[, 4] < 0.05,]

length(setdiff(rownames(cox_results[cox_results[, 4] < 0.05,]),
               names(log_rank_p[log_rank_p < 0.05])
))
length(setdiff(names(log_rank_p[log_rank_p < 0.05]),
               rownames(cox_results[cox_results[, 4] < 0.05,])
))
length(unique(names(log_rank_p[log_rank_p < 0.05]),
              rownames(cox_results[cox_results[, 4] < 0.05,])
))
save(log_rank_p,cox_results,exprSet,phe,file = 'surviva.Rdata')
```
