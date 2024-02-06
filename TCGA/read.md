## TCGA数据

```properties
## 背景:
https://portal.gdc.cancer.gov/

```

#### 基本流程

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
```
```R
## 3.


```