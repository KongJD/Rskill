## TCGA数据

```properties
## 背景:
https://portal.gdc.cancer.gov/

```

#### 基本流程

```properties
后续代码后续测试跑通，代码这里不是关键，思想、流程学会
```

```R
## 1.TCGA数据下载，比较繁琐，后续需要的话再查一下
library(RTCGA.miRNASeq)
s = rownames(KIRC.miRNASeq)[seq(1, nrow(KIRC.miRNASeq), by = 3)]
expr <- expressionsTCGA(KIRC.miRNASeq)
expr = as.data.frame(expr[seq(1, nrow(expr), by = 3), 3:ncol(expr)])
mi = colnames(expr)
expr = apply(expr, 1, as.numeric)
colnames(expr) = s
rownames(expr) = mi
expr = na.omit(expr)
expr = expr[apply(expr, 1, function(x) { sum(x > 1) > 10 }),]

library(RTCGA.clinical)
meta <- KIRC.clinical
tmp = as.data.frame(colnames(meta))
meta[(grepl('patient.bcr_patient_barcode', colnames(meta)))]
meta[(grepl('patient.days_to_last_followup', colnames(meta)))]
meta[(grepl('patient.days_to_death', colnames(meta)))]
meta[(grepl('patient.vital_status', colnames(meta)))]
meta = as.data.frame(meta[c('patient.bcr_patient_barcode', 'patient.vital_status',
                            'patient.days_to_death', 'patient.days_to_last_followup',
                            'patient.race',
                            'patient.age_at_initial_pathologic_diagnosis',
                            'patient.gender',
                            'patient.stage_event.pathologic_stage')])
```

```R
## 2.拿到表达矩阵，差异分析
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
exprSet = na.omit(expr)

draw_h_v <- function(exprSet, need_DEG, n = 'DEseq2', group_list, logFC_cutoff) {
  library(pheatmap)
  choose_gene = head(rownames(need_DEG), 50)
  choose_matrix = exprSet[choose_gene,]
  choose_matrix[1:4, 1:4]
  choose_matrix = t(scale(t(log2(choose_matrix + 1))))
  annotation_col = data.frame(group_list = group_list)
  rownames(annotation_col) = colnames(exprSet)
  pheatmap(choose_matrix, show_colnames = F, annotation_col = annotation_col,
           filename = paste0(n, '_need_DEG_top50_heatmap.png'))
  library(ggfortify)
  df = as.data.frame(t(choose_matrix))
  df$group = group_list
  png(paste0(n, '_DEG_top50_pca.png'), res = 120)
  p = autoplot(prcomp(df[, 1:(ncol(df) - 1)]), data = df, colour = 'group') + theme_bw()
  dev.off()
  if (!logFC_cutoff) {
    logFC_cutoff <- with(need_DEG, mean(abs(log2FoldChange)) + 2 * sd(abs(log2FoldChange)))
  }
  need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff, 'UP', 'DOWN'), 'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ', round(logFC_cutoff, 3),
                      '\nThe number of up gene is ', nrow(need_DEG[need_DEG$change == 'UP',]),
                      '\nThe number of down gene is ', nrow(need_DEG[need_DEG$change == 'DOWN',])
  )
  library(ggplot2)
  g = ggplot(data = need_DEG,
             aes(x = log2FoldChange, y = -log10(pvalue),
                 color = change)) +
    geom_point(alpha = 0.4, size = 1.75) +
    theme_set(theme_set(theme_bw(base_size = 20))) +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    ggtitle(this_tile) +
    theme(plot.title = element_text(size = 15, hjust = 0.5)) +
    scale_colour_manual(values = c('blue', 'black', 'red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g, filename = paste0(n, '_volcano.png'))
  dev.off()
}

### DESeq2
library(DESeq2)
(colData <- data.frame(row.names = colnames(exprSet),
                       group_list = group_list))
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~group_list)
dds <- DESeq(dds)
res <- results(dds,
               contrast = c("group_list", "tumor", "normal"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG = as.data.frame(resOrdered)
DESeq2_DEG = na.omit(DEG)
nrDEG = DESeq2_DEG[, c(2, 6)]
colnames(nrDEG) = c('log2FoldChange', 'pvalue')
draw_h_v(exprSet, nrDEG, 'DEseq2', group_list, 1)

## edgeR
library(edgeR)
d <- DGEList(counts = exprSet, group = factor(group_list))
keep <- rowSums(cpm(d) > 1) >= 2
d <- d[keep, , keep.lib.sizes = FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
dge = d
design <- model.matrix(~0 + factor(group_list))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group_list))
dge = d
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(-1, 1))
nrDEG = topTags(lrt, n = nrow(dge))
nrDEG = as.data.frame(nrDEG)
edgeR_DEG = nrDEG
nrDEG = edgeR_DEG[, c(1, 5)]
colnames(nrDEG) = c('log2FoldChange', 'pvalue')
draw_h_v(exprSet, nrDEG, 'edgeR', group_list, 1)

## limma
suppressMessages(library(limma))
design <- model.matrix(~0 + factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(exprSet)
dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)
v <- voom(dge, design, plot = TRUE, normalize = "quantile")
fit <- lmFit(v, design)
cont.matrix = makeContrasts(contrasts = c('tumor-normal'), levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)
tempOutput = topTable(fit2, coef = 'tumor-normal', n = Inf)
DEG_limma_voom = na.omit(tempOutput)
nrDEG = DEG_limma_voom[, c(1, 4)]
colnames(nrDEG) = c('log2FoldChange', 'pvalue')
draw_h_v(exprSet, nrDEG, 'limma', group_list, 1)
```

```R
## 3.生存分析
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
exprSet = na.omit(expr)
library(survival)
library(survminer)
# 这里做生存分析，已经不需要正常样本的表达矩阵了，所以需要过滤。
# 而且临床信息，有需要进行整理。
exprSet = na.omit(expr)
exprSet = exprSet[, group_list == 'tumor']
meta[, 3][is.na(meta[, 3])] = 0
meta[, 4][is.na(meta[, 4])] = 0
meta$days = as.numeric(meta[, 3]) + as.numeric(meta[, 4])
meta = meta[, c(1:2, 5:9)]
colnames(meta)
colnames(meta) = c('ID', 'event', 'race', 'age', 'gender', 'stage', "days")


library(survival)
library(survminer)
meta$event = ifelse(meta$event == 'alive', 0, 1)
meta$age = as.numeric(meta$age)
library(stringr)
meta$stage = str_split(meta$stage, ' ', simplify = T)[, 2]
boxplot(meta$age)
meta$age_group = ifelse(meta$age > median(meta$age), 'older', 'younger')

meta$time = meta$days / 30
phe = meta
phe$ID = toupper(phe$ID)
phe = phe[match(substr(colnames(exprSet), 1, 12), phe$ID),]

save(exprSet, phe,
     file =
       file.path(Rdata_dir, 'TCGA-KIRC-miRNA-survival_input.Rdata')
)

## 批量生存分析
library(survival)
library(survminer)
## 批量生存分析 使用 coxph 回归方法 http://www.sthda.com/english/wiki/cox-proportional-hazards-model
colnames(phe)
mySurv = with(phe, Surv(time, event))
## 有些基因恒定的或者表达都为0在所有样本中的直接删除
cox_results <- apply(exprSet, 1, function(gene) {
  group = ifelse(gene > median(gene), 'high', 'low')
  survival_dat <- data.frame(group = group, stage = phe$stage, age = phe$age,
                             gender = phe$gender,
                             stringsAsFactors = F)
  m = coxph(mySurv ~ gender + age + stage + group, data = survival_dat)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
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

## 批量生存分析 使用  logrank test 方法
mySurv = with(phe, Surv(time, event))
log_rank_p <- apply(exprSet, 1, function(gene) {
  phe$group = ifelse(gene > median(gene), 'high', 'low')
  data.survdiff = survdiff(mySurv ~ group, data = phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

require("VennDiagram")
VENN.LIST = list(cox = rownames(cox_results[cox_results[, 4] < 0.05,]),
                 log = names(log_rank_p[log_rank_p < 0.05]))
venn.plot <- venn.diagram(VENN.LIST, NULL,
                          fill = c("darkmagenta", "darkblue"),
                          alpha = c(0.5, 0.5), cex = 2,
                          cat.fontface = 4,
                          main = "overlap of coxph and log-rank test")
grid.draw(venn.plot)

library(pheatmap)
choose_gene = rownames(cox_results[cox_results[, 4] < 0.05,])
choose_matrix = expr[choose_gene,]
n = t(scale(t(log2(choose_matrix + 1))))
n[n > 2] = 2
n[n < -2] = -2

annotation_col = data.frame(group_list = group_list)
rownames(annotation_col) = colnames(expr)

pheatmap(n, show_colnames = F, annotation_col = annotation_col,
         filename = 'cox_genes.heatmap.png')

library(ggfortify)
df = as.data.frame(t(choose_matrix))
df$group = group_list
png('cox_genes.pca.png', res = 120)
autoplot(prcomp(df[, 1:(ncol(df) - 1)]), data = df, colour = 'group') + theme_bw()
dev.off()


library("FactoMineR")
library("factoextra")
dat.pca <- PCA(t(choose_matrix), graph = FALSE)
fviz_pca_ind(dat.pca, repel = T,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
```

```R
## 5.lasso 回归问题
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
exprSet = na.omit(expr)
exprSet = na.omit(expr)
all(substring(colnames(exprSet), 1, 12) == phe$ID)
library(lars)
library(glmnet)
x = t(log2(exprSet + 1))
y = phe$event
model_lasso <- glmnet(x, y, family = "binomial", nlambda = 50, alpha = 1)
plot.glmnet(model_lasso, xvar = "norm", label = TRUE)
plot(model_lasso, xvar = "lambda", label = TRUE)
cv_fit <- cv.glmnet(x = x, y = y, alpha = 1, nlambda = 1000)
plot.cv.glmnet(cv_fit)
c(cv_fit$lambda.min, cv_fit$lambda.1se)

model_lasso <- glmnet(x = x, y = y, alpha = 1, lambda = cv_fit$lambda.1se)
lasso.prob <- predict(cv_fit, newx = x, s = c(cv_fit$lambda.min, cv_fit$lambda.1se))
re = cbind(y, lasso.prob)
dat = as.data.frame(re[, 1:2])
colnames(dat) = c('event', 'prob')
dat$event = as.factor(dat$event)
library(ggpubr)
p <- ggboxplot(dat, x = "event", y = "prob",
               color = "event", palette = "jco",
               add = "jitter")
p + stat_compare_means()

library(ROCR)
library(glmnet)
library(caret)

pred <- prediction(re[, 2], re[, 1])
perf <- performance(pred, "tpr", "fpr")
performance(pred, "auc")
plot(perf, colorize = FALSE, col = "black")
lines(c(0, 1), c(0, 1), col = "gray", lty = 4)


fit <- glmnet(x = x, y = y, alpha = 1, lambda = cv_fit$lambda.1se)
choose_gene = rownames(fit$beta)[as.numeric(fit$beta) != 0]
length(choose_gene)
myexpr = x[, choose_gene]
mysurv = phe[, c("days", "event")]
mysurv$days[mysurv$days < 1] = 1
fit <- glmnet(myexpr, Surv(mysurv$days, mysurv$event),
              family = "cox")
plot(fit, xvar = "lambda", label = TRUE)
plot(fit, label = TRUE)

library(pheatmap)
choose_matrix = expr[choose_gene,]
choose_matrix[1:4, 1:4]
n = t(scale(t(log2(choose_matrix + 1))))
n[n > 2] = 2
n[n < -2] = -2
annotation_col = data.frame(group_list = group_list)
rownames(annotation_col) = colnames(expr)
pheatmap(n, show_colnames = F, annotation_col = annotation_col,
         filename = 'lasso_genes.heatmap.png')

library(ggfortify)
df = as.data.frame(t(choose_matrix))
df$group = group_list
png('lasso_genes.pca.png', res = 120)
autoplot(prcomp(df[, 1:(ncol(df) - 1)]), data = df, colour = 'group') + theme_bw()
dev.off()

library("FactoMineR")
library("factoextra")
dat.pca <- PCA(t(choose_matrix), graph = FALSE)
fviz_pca_ind(dat.pca, repel = T,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)
```

```R
## 6.coxpt-fortest
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
exprSet = na.omit(expr)
exprSet = na.omit(expr)
all(substring(colnames(exprSet), 1, 12) == phe$ID)
## 挑选感兴趣的基因构建coxph模型
e = t(exprSet[c('hsa-mir-21', 'hsa-mir-143', 'hsa-mir-10b', 'hsa-mir-192', 'hsa-mir-183'),])
e = log2(e)
colnames(e) = c('miR21', 'miR143', 'miR10b', 'miR192', 'miR183')
dat = cbind(phe, e)
dat$gender = factor(dat$gender)
dat$stage = factor(dat$stage)
s = Surv(time, event) ~ gender +
  stage +
  age +
  miR21 +
  miR143 +
  miR10b +
  miR192 +
  miR183
s = Surv(time, event) ~ miR21 + miR143 + miR10b + miR192 + miR183
model <- coxph(s, data = dat)
summary(model, data = dat)
options(scipen = 1)
ggforest(model, data = dat,
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0,
         refLabel = "1", noDigits = 4)

fp <- predict(model)
summary(model, data = dat)
library(Hmisc)
options(scipen = 200)
with(dat, rcorr.cens(fp, Surv(time, event)))
```

```R
## 7. risk-score-dtribution
library(survival)
library(survminer)


library(RTCGA.miRNASeq)
s = rownames(LUAD.miRNASeq)[seq(1, nrow(LUAD.miRNASeq), by = 3)]
expr <- expressionsTCGA(LUAD.miRNASeq)
expr = as.data.frame(expr[seq(1, nrow(expr), by = 3), 3:ncol(expr)])
mi = colnames(expr)
expr = apply(expr, 1, as.numeric)
colnames(expr) = s
rownames(expr) = mi
expr = na.omit(expr)
expr = expr[apply(expr, 1, function(x) { sum(x > 1) > 10 }),]

library(RTCGA.clinical)
meta <- LUAD.clinical
tmp = as.data.frame(colnames(meta))
meta[(grepl('patient.bcr_patient_barcode', colnames(meta)))]
meta[(grepl('patient.days_to_last_followup', colnames(meta)))]
meta[(grepl('patient.days_to_death', colnames(meta)))]
meta[(grepl('patient.vital_status', colnames(meta)))]
meta = as.data.frame(meta[c('patient.bcr_patient_barcode', 'patient.vital_status',
                            'patient.days_to_death', 'patient.days_to_last_followup',
                            'patient.race',
                            'patient.age_at_initial_pathologic_diagnosis',
                            'patient.gender',
                            'patient.stage_event.pathologic_stage')])

group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')


exprSet = na.omit(expr)
exprSet = exprSet[, group_list == 'tumor']
meta[, 3][is.na(meta[, 3])] = 0
meta[, 4][is.na(meta[, 4])] = 0
meta$days = as.numeric(meta[, 3]) + as.numeric(meta[, 4])
meta = meta[, c(1:2, 5:9)]
colnames(meta) = c('ID', 'event', 'race', 'age', 'gender', 'stage', "days")
library(survival)
library(survminer)
meta$event = ifelse(meta$event == 'alive', 0, 1)
meta$age = as.numeric(meta$age)
library(stringr)
meta$stage = str_split(meta$stage, ' ', simplify = T)[, 2]
boxplot(meta$age)
meta$age_group = ifelse(meta$age > median(meta$age, na.rm = T), 'older', 'younger')
meta$time = meta$days / 30
phe = meta
phe$ID = toupper(phe$ID)
phe = phe[match(substr(colnames(exprSet), 1, 12), phe$ID),]


# 这个时候是515个病人的673个miRNA表达矩阵。
## 挑选感兴趣的基因构建coxph模型
e = t(exprSet[c('hsa-mir-31', 'hsa-mir-196b', 'hsa-mir-766', 'hsa-mir-519a-1', 'hsa-mir-375', 'hsa-mir-187', 'hsa-mir-331', 'hsa-mir-101-1'),])
e = log2(e + 1)
colnames(e) = c('miR31', 'miR196b', 'miR766', 'miR519a1', 'miR375', 'miR187', 'miR331', 'miR101')
dat = cbind(phe, e)
dat$gender = factor(dat$gender)
dat$stage = factor(dat$stage)

s = Surv(time, event) ~ miR31 +
  miR196b +
  miR766 +
  miR519a1 +
  miR375 +
  miR187 +
  miR331 +
  miR101
model <- coxph(s, data = dat)
summary(model, data = dat)
options(scipen = 1)
ggforest(model, data = dat,
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0,
         refLabel = "1", noDigits = 4)


# 对于生存数据预后模型评价很多采用C-index ，但c-index展示没有roc曲线图来的直观
new_dat = dat
fp <- predict(model, new_dat, type = "risk"); boxplot(fp)
fp <- predict(model, new_dat, type = "expected"); boxplot(fp)
plot(fp, phe$days)
fp <- predict(model, new_dat); boxplot(fp)
basehaz(model)
library(Hmisc)
options(scipen = 200)
with(new_dat, rcorr.cens(fp, Surv(time, event)))

library(cowplot)
library(pheatmap)
fp_dat = data.frame(s = 1:length(fp), v = as.numeric(sort(fp)))
sur_dat = data.frame(s = 1:length(fp),
                     t = phe[names(sort(fp)), 'time'],
                     e = phe[names(sort(fp)), 'event'])
sur_dat$e = ifelse(sur_dat$e == 0, 'alive', 'death')
exp_dat = new_dat[names(sort(fp)), 10:17]
plot.point = ggplot(fp_dat, aes(x = s, y = v)) + geom_point(); print(plot.point)
plot.sur = ggplot(sur_dat, aes(x = s, y = t)) + geom_point(aes(col = e)); print(plot.sur)

mycolors <- colorRampPalette(c("black", "green", "red"), bias = 1.2)(100)
tmp = t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
plot.h = pheatmap(tmp, col = mycolors, show_colnames = F, cluster_cols = T)
plot.h = pheatmap(tmp, col = mycolors, show_colnames = F, cluster_cols = F)
plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("A", "B", "C"),
          align = 'v', ncol = 1)
```

```R
## 8.random_forest
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
exprSet = na.omit(expr)
exprSet = na.omit(expr)
all(substring(colnames(exprSet), 1, 12) == phe$ID)

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
x = t(log2(exprSet + 1))
y = phe$event

tmp = as.vector(table(y))
num_classes = length(tmp)
min_size = tmp[order(tmp, decreasing = FALSE)[1]]
sampsizes = rep(min_size, num_classes)

rf_output = randomForest(x = x, y = y, importance = TRUE, ntree = 10001, proximity = TRUE)
rf_importances = importance(rf_output, scale = FALSE)

varImpPlot(rf_output, type = 2, n.var = 30, scale = FALSE,
           main = "Variable Importance (Gini) for top 30 predictors", cex = .7)
target_labels = as.vector(y)
MDSplot(rf_output, y, k = 2, xlab = "", ylab = "",
        pch = target_labels, palette = c("red", "blue"), main = "MDS plot")


choose_gene = rownames(tail(rf_importances[order(rf_importances[, 2]),], 50))

library(pheatmap)
choose_matrix = expr[choose_gene,]
n = t(scale(t(log2(choose_matrix + 1))))
n[n > 2] = 2
n[n < -2] = -2
annotation_col = data.frame(group_list = group_list)
rownames(annotation_col) = colnames(expr)

pheatmap(n, show_colnames = F, annotation_col = annotation_col,
         filename = 'rf_genes.heatmap.png')

library(ggfortify)
df = as.data.frame(t(choose_matrix))
df$group = group_list
png('rf_genes.pca.png', res = 120)
autoplot(prcomp(df[, 1:(ncol(df) - 1)]), data = df, colour = 'group') + theme_bw()
dev.off()

library("FactoMineR")
library("factoextra")
dat.pca <- PCA(t(choose_matrix), graph = FALSE)
fviz_pca_ind(dat.pca, repel = T,
             geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
             col.ind = group_list, # color by groups 颜色组
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses 集中成椭圆
             legend.title = "Groups"
)
```

```R
## 9.maftools
laml = read.maf(maf = '../GDC/TCGA.KIRC.mutect.somatic.maf.gz')
project = 'TCGA_KIRC'
laml@data = laml@data[!grepl('^MT-', laml@data$Hugo_Symbol),]
laml@data$t_vaf = (laml@data$t_alt_count / laml@data$t_depth)
png(paste0('plotmafSummary_', project, '.png'), res = 150, width = 1080, height = 1080)
plotmafSummary(maf = laml, rmOutlier = TRUE, showBarcodes = T,
               addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(paste0('oncoplot_top30_', project, '.png'), res = 150, width = 1080, height = 1080)
oncoplot(maf = laml, top = 30, fontSize = 12, showTumorSampleBarcodes = F)
dev.off()

oncoplot(maf = laml, top = 15, fontSize = 12,
         clinicalFeatures = c("subtype"), sortByAnnotation = T)

png(paste0('TMB_', project, '.png'), res = 150, width = 1080, height = 1080)
laml.mutload = tcgaCompare(maf = laml, cohortName = project)
dev.off()

png(paste0('Vaf_', project, '.png'), res = 150, width = 1080, height = 1080)
plotVaf(maf = laml, top = 20)
dev.off()


dir.create(paste0('vaf_clust_', project))
lapply(unique(laml@data$Tumor_Sample_Barcode), function(x) {
  png(paste0('vaf_clust_', project, '/', x, '_vaf_clust.png'), res = 120, width = 1080, height = 1080)
  het = inferHeterogeneity(maf = laml, tsb = x, vafCol = 't_vaf')
  print(het$clusterMeans)
  plotClusters(clusters = het)
  dev.off()
})

laml@data$t_vaf = (laml@data$t_alt_count / laml@data$t_depth)
mut = laml@data[, c("Hugo_Symbol", "Chromosome", "Start_Position", "Tumor_Sample_Barcode", "t_vaf")]
mut$pos = paste(mut$Chromosome, mut$Start_Position, sep = ':')
```

```R
## 10.boxplot
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
dat = data.frame(gene = log2(exprSet['hsa-mir-10b',] + 1),
                 stage = phe$stage)
boxplot(dat$gene ~ dat$stage)
if (require('ggpubr')) {
  library(ggpubr)
  p <- ggboxplot(dat, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  p + stat_compare_means()
}
if (require('ggstatsplot')) {
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = stage, y = gene)
}
res.aov <- aov(gene ~ stage, data = dat)
VHL_mut = substr(as.character(
  as.data.frame(mut[mut$Hugo_Symbol == 'VHL', 'Tumor_Sample_Barcode'])[, 1]),
                 1, 12)
library(dplyr)
mut %>%
  filter(Hugo_Symbol == 'VHL') %>%
  as.data.frame() %>%
  pull(Tumor_Sample_Barcode) %>%
  as.character() %>%
  substr(1, 12)
dat = data.frame(gene = log2(exprSet['hsa-mir-10b',]),
                 mut = substr(colnames(exprSet), 1, 12) %in% VHL_mut)
if (require('ggpubr')) {
  library(ggpubr)
  p <- ggboxplot(dat, x = "mut", y = "gene",
                 color = "mut", palette = "jco",
                 add = "jitter")
  p + stat_compare_means(method = "t.test")
}
if (require('ggstatsplot')) {
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = mut, y = gene)
}
if (require('ggplot2')) {
  library(ggplot2)
  ggplot(dat, aes(x = mut, y = gene)) +
    geom_boxplot() +
    geom_jitter() +
    geom_violin() +
    theme_bw()
}

```

```R
## 11.correlion
group_list = ifelse(as.numeric(substr(colnames(expr), 14, 15)) < 10, 'tumor', 'normal')
dat = data.frame(gene1 = log2(exprSet['hsa-mir-10b',] + 1),
                 gene2 = log2(exprSet['hsa-mir-143',] + 1),
                 stage = phe$stage)
library(ggpubr)
dat$stage = as.factor(dat$stage)
sp <- ggscatter(dat, x = "gene1", y = "gene2",
                add = "reg.line",  # Add regressin line 
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
sp + stat_cor(method = "pearson", label.x = 15, label.y = 20)
sp <- ggscatter(dat, x = "gene1", y = "gene2",
                color = "stage", palette = "jco",
                add = "reg.line", conf.int = TRUE)
sp + stat_cor(aes(color = stage), label.x = 15)
```

```R
## 12 split-corr

```

```R
## 13 roc
```

```R
## 14.choose_lncrna
```