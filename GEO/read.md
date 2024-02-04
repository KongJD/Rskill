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
tmp <- by(dat_final, ids$symbol, function (x) {
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

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                               levels = design)

deg <- function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  return(nrDEG)
}

res <- deg(dat_final_res,design,contrast.matrix)

## volcano 
library(ggpubr)
plot(res$logFC,-log10(res$P.Value))
res$v= -log10(res$P.Value)
ggscatter(res, x = "logFC", y = "v",size=0.5)
res$g=ifelse(res$P.Value>0.01,'stable', 
            ifelse( res$logFC >2,'up', 
                    ifelse( res$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
res$name=rownames(res)
ggscatter(res, x = "logFC", y = "v",size=0.5,color = 'g')
ggscatter(res, x = "logFC", y = "v", color = "g",size = 0.5,
          label = "name", repel = T,
          #label.select = rownames(df)[df$g != 'stable'] ,
          label.select = head(rownames(res)), #挑选一些基因在图中显示出来
          palette = c("#00AFBB", "#E7B800", "#FC4E07") )

ggscatter(res, x = "AveExpr", y = "logFC",size = 0.2)
res$p_c = ifelse(res$P.Value<0.001,'p<0.001',
                ifelse(res$P.Value<0.01,'0.001<p<0.01','p>0.01'))
ggscatter(res,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
          palette = c("green", "red", "black") )


### pheatmap
library(pheatmap)

x=res$logFC 
names(x)=rownames(res) 
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))
pheatmap(dat[cg,],show_colnames =F,show_rownames = F)
n=t(scale(t(dat[cg,]))) #通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
pheatmap(n,show_colnames =F,show_rownames = F)
n[n>2]=2
n[n< -2]= -2
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
pheatmap(n,show_colnames =F,
         show_rownames = F,
         cluster_cols = F, 
         annotation_col=ac,) #列名注释信息为ac即分组信息
```

```R
## 4.Go和Kegg


```