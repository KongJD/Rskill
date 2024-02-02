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
        https://f1000research.com/articles/5-1438/v1 
        
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
```

#### 差异分析

```R

 suppressMessages(library(DESeq2)) 
 colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) 
 dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)
 
  res <- results(dds, contrast=c("group_list","treat_2","control"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_treat_2=as.data.frame(resOrdered)
  write.csv(DEG_treat_2,"DEG_treat_2_deseq2.results.csv")

```
