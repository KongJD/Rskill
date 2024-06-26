## 转录组
```properties
背景资料 
# https://mp.weixin.qq.com/s/a3y46NNNO-wardO3XWwh0w
```

#### 基本流程
```properties
流程：
1. 参考基因组、基因注释文件 下载 
2. 质控，需要fastqc及multiqc等
3. 比对：hisat2
4. 计数和归一化
5. 差异分析


流程参考：https://f1000research.com/articles/4-1070/v1
        
```

```shell
step1:
cat SRR_Acc_List.txt |while read id ;do (prefetch  ${id} &);done
ls /public/project/RNA/airway/sra/*  |while read id;do ( nohup fastq-dump --gzip --split-3 -O ./ ${id} & );done
```

```shell
step2:
ls *gz |xargs fastqc -t 10
multiqc ./ 
```

```shell
step3: 去除低质量的和有接头的
trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir  $fq1 $fq2 &
```

```shell
step4:
lls *gz|cut -d"_" -f 1 |sort -u |while read id;do
ls -lh ${id}_1_val_1.fq.gz   ${id}_2_val_2.fq.gz 
hisat2 -p 10 -x /public/reference/index/hisat/hg38/genome -1 ${id}_1_val_1.fq.gz   -2 ${id}_2_val_2.fq.gz  -S ${id}.hisat.sam
done

ls *.sam|while read id ;do (samtools sort -O bam -@ 5  -o $(basename ${id} ".sam").bam   ${id});done
rm *.sam 
ls *.bam |xargs -i samtools index {}
# ls *.bam |xargs -i samtools flagstat -@ 10  {}  > 
ls *.bam |while read id ;do ( nohup samtools flagstat -@ 1 $id >  $(basename ${id} ".bam").flagstat  & );done

```

```shell
step5:
gtf="/public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz"   
featureCounts -T 5 -p -t exon -g gene_id  -a $gtf -o  all.id.txt  1>counts.id.log 2>&1 &
|htseq-count -f bam -r pos ../clean/sample.bam  /public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz > SRR1039508.count.txt 
```

#### 差异分析

```R
### ---------------
###
### Firstly run DEseq2 
###
### ---------------
suppressMessages(library(DESeq2))
colData <- data.frame(row.names = colnames(exprSet), group_list = group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~group_list)
dds <- DESeq(dds)

res <- results(dds, contrast = c("group_list", "treat_2", "control"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG_treat_2 = as.data.frame(resOrdered)
DEG_treat_2 = na.omit(DEG_treat_2)

write.csv(DEG_treat_2, "DEG_treat_2_deseq2.results.csv")


# 热图
library(pheatmap)
choose_gene = head(rownames(need_DEG), 50) ## 50 maybe better
choose_matrix = exprSet[choose_gene,]
choose_matrix = t(scale(t(choose_matrix)))
pheatmap(choose_matrix, filename = paste0(n, '_need_DEG_top50_heatmap.png'))


logFC_cutoff <- with(need_DEG, mean(abs(log2FoldChange)) + 2 * sd(abs(log2FoldChange)))
# logFC_cutoff=1

need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                   ifelse(need_DEG$log2FoldChange > logFC_cutoff, 'UP', 'DOWN'), 'NOT')
)
this_tile <- paste0('Cutoff for logFC is ', round(logFC_cutoff, 3),
                    '\nThe number of up gene is ', nrow(need_DEG[need_DEG$change == 'UP',]),
                    '\nThe number of down gene is ', nrow(need_DEG[need_DEG$change == 'DOWN',])
)

#火山图
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


if (F) {
  png("DESeq2_qc_dispersions.png", 1000, 1000, pointsize = 20)
  plotDispEsts(dds, main = "Dispersion plot")
  dev.off()

  rld <- rlogTransformation(dds)
  exprMatrix_rlog = assay(rld)

  x = apply(exprMatrix_rlog, 1, mean)
  y = apply(exprMatrix_rlog, 1, mad)
  plot(x, y)

  png("DESeq2_RAWvsNORM.png", height = 800, width = 800)
  par(cex = 0.7)
  n.sample = ncol(exprSet)
  if (n.sample > 40) par(cex = 0.5)
  cols <- rainbow(n.sample * 1.2)
  par(mfrow = c(2, 2))
  boxplot(exprSet, col = cols, main = "expression value", las = 2)
  boxplot(exprMatrix_rlog, col = cols, main = "expression value", las = 2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()

}

```

```R
### ---------------
###
### Then run edgeR 
###
### ---------------
library(edgeR)
g = factor(group_list)
g = relevel(g, g1)
d <- DGEList(counts = exprSet, group = g)
keep <- rowSums(cpm(d) > 1) >= 2
table(keep)
d <- d[keep, , keep.lib.sizes = FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples

dge = d
design <- model.matrix(~0 + factor(group_list))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group_list))

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)

lrt <- glmLRT(fit, contrast = c(-1, 1))
nrDEG = topTags(lrt, n = nrow(dge))
nrDEG = as.data.frame(nrDEG)
head(nrDEG)
# DEG_edgeR = nrDEG
# nrDEG = DEG_edgeR[, c(1, 5)]
# colnames(nrDEG) = c('log2FoldChange', 'pvalue')
# draw_h_v(exprSet, nrDEG, paste0(pro, '_edgeR'))
```

```R
### ---------------
###
### Then run limma 
###
### --------------- 
suppressMessages(library(limma))
design <- model.matrix(~0 + factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(exprSet)
design

dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

v <- voom(dge, design, plot = TRUE, normalize = "quantile")
fit <- lmFit(v, design)

group_list
con = paste0(g2, '-', g1)
cat(con)
cont.matrix = makeContrasts(contrasts = c(con), levels = design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

tempOutput = topTable(fit2, coef = con, n = Inf)
DEG_limma_voom = na.omit(tempOutput)
head(DEG_limma_voom)
```

