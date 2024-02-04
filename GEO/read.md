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

```
