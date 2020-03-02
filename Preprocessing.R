###################
#   TCGA-BRCA
#   Pancheng Wu
#   2020-01-15
#   PUMC
###################

library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(readr)
library(plyr)
library(DESeq)
library(DESeq2)
library(edgeR)
a=fread('TCGA-BRCA.htseq_counts.tsv.gz',sep = '\t',header = T)
class(a)
a=as.data.frame(a)
a[1:4,1:4] 
rownames(a)=a[,1]
a=a[,-1]
genes=rownames(a)
a[1:4,1:4]
a=2^a-1
a[1:4,1:4]
a <- ceiling(a)
a[1:4,1:4]
b <- a ##作为备份

## 去除极低表达量基因,暂时可不用
geneLists <- rownames(a)
keepGene = rowSums(edgeR::cpm(a)>0) >=2
table(keepGene)
dim(a)
dim(a[keepGene,])
a=a[keepGene,]
rownames(a)=geneLists[keepGene]

## 这个基因需要在超过1/3的样本中大于0
index <- rowSums(b > 0) > (ncol(b) / 3) ## 去除25446个，保留35042个
table(index)
b <- b[index,]
## 基因在每个样本中平均表达量为1就被过滤，暂时不用
index1 <- rowSums(b) < ncol(b)
table(index1)

## 准备标准化
exprSet <- cbind("gene_id" = rownames(b), b)
rownames(exprSet) <- NULL
exprSet$gene_id <- as.character(exprSet$gene_id)
## 查看尾部，发现有问题
tail(exprSet$gene_id,10)
## 删除3个,保留35039个 gene symbol
exprSet <- exprSet[1:(nrow(exprSet)-3),]
## 去除 点号
expr_df_nopoint <- exprSet %>% 
  tidyr::separate(gene_id,into = c("gene_id"),sep="\\.") 
save(expr_df_nopoint, file = "expr_df_nopoint_nonlog.Rda") ## 为去log化的数据

## 提取 mRNA 数据
load("C:/Rlearn/tongyong/gtf_df.Rda")
exprSet_m <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>% #筛选gene,和编码指标
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(expr_df_nopoint,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_biotype,sep = " | ")
write.csv(exprSet_m, file = "exprSet_m.csv")
## 分离
gene <- exprSet_m$gene_id
gene1 <-  strsplit(gene, "|", fixed = T)
gene2 <- as.data.frame(gene1)
gene3 <- gene2[1,]
gene3 <- t(gene3)
gene3 <- as.data.frame(gene3)
gene3 <- gene3$`1`
gene3 <- as.character(gene3)

exprSet_m <- exprSet_m %>%
  group_by(gene_id) %>%
  summarise_all(max)
save(exprSet_m, file = "exprSet_m_LUAD_nolog.Rda")
rm(list = ls())
## Desq2 标准化
## 准备metadata 文件

library(DESeq2)
TCGA_id <- colnames(exprSet_m)[-1]
table(substring(TCGA_id,14,15))
TCGA_id <- TCGA_id[substring(TCGA_id, 14,15) != "06"]
gene_id <- exprSet_m$gene_id
exprSet_m <- exprSet_m[,colnames(exprSet_m) %in% TCGA_id]
exprSet_m <- cbind("gene_id" = gene_id, exprSet_m)
exprSet_m$gene_id <- as.character(exprSet_m$gene_id)
metadata <- data.frame(TCGA_id = colnames(exprSet_m)[-1])
sample <- ifelse(substring(metadata$TCGA_id,14,15)=="01","Cancer","Normal")
## 这里的factor是为了dds的需要
metadata$sample <- as.factor(sample)
exprSet_m[,c(2:1211)] <- ceiling(exprSet_m[,c(2:1211)])
## 
mycounts <- exprSet_m
dds <-DESeqDataSetFromMatrix(countData=mycounts, 
                             colData=metadata, 
                             design=~sample,
                             tidy=TRUE)

dds <- DESeq(dds)

save(dds,file="exprSet_m_dds_sample.Rdata")
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, "sample")

exprSet_m_vst <- as.data.frame(assay(vsd))
library(stringr)
e <- exprSet_m_vst
gene <- rownames(exprSet_m_vst)
gene1 <-  strsplit(gene, "|", fixed = T)
gene2 <- as.data.frame(gene1)
gene3 <- gene2[1,]
gene3 <- t(gene3)
gene3 <- as.data.frame(gene3)
gene3 <- gene3$`1`
gene3 <- as.character(gene3)
gene3 <- str_squish(gene3)
expr_m_vst <- exprSet_m_vst
rownames(expr_m_vst) <- gene3
exprSet_m_vst <- cbind("gene" = gene3, exprSet_m_vst)
rownames(exprSet_m_vst) <- NULL
exprSet_m_vst$gene <- as.character(exprSet_m_vst$gene)
exprSet_m_vst <- exprSet_m_vst[!duplicated(exprSet_m_vst$gene),] ## 17917到17915个gene symbol 
save(exprSet_m_vst,file = "exprSet_m_vst_17915.Rda")

##
load("exprSet_m_vst_17915.Rda")


