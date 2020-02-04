###################
# Downloading GEo data sets
# 2020-02-04
# Pancheng Wu
# PUMC
#
##################
##

## 
all <- read.csv("breast_cancer_survival.csv")
all_GPL570 <- subset(all, all$gpl == "GPL570")
all_GPL96 <- subset(all, all$gpl == "GPL96")


## GEO 数据集下载

## GSE21653 85个LUAD样本，61个LUSC样本 // DFS
GSE21653 {
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE21653", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ##
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  #save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
} 

## GSE20685
GSE20685{
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE20685", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ##
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  #save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
}

## GSE20711 
GSE20711{
  library(GEOquery)
  library(stringr)
  library(limma)
  library(dplyr)
  library('Biobase')
  library(ggstatsplot)
  library(tidyr)
  library(ggpubr)
  
  gset = getGEO("GSE20711", destdir = ".", getGPL = F)
  class(gset)
  str(gset)
  gset
  gset <- gset[[1]]
  ## 表达矩阵
  exprSet = exprs(gset)
  a <- exprSet
  exprSet <- normalizeBetweenArrays(exprSet)  ## 标准化一番
  exprSet <- as.data.frame(exprSet) ## 转为数据框
  gene_id <- rownames(exprSet)
  exprSet <- cbind("probe_id" = gene_id, exprSet)
  exprSet$probe_id <- as.character(exprSet$probe_id)
  
  ##
  load("C:/Rlearn/GEO/GPL-probe2symbol/probe2symbol-GPL570.Rda")
  exprSet <- exprSet%>%
    inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
    select(-probe_id) %>% #去掉多余信息
    select(symbol, everything()) %>% #重新排列，
    mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #求出平均数(这边的.真的是画龙点睛)
    arrange(desc(rowMean)) %>% #把表达量的平均值按从大到小排序
    distinct(symbol,.keep_all = T) %>% # symbol留下第一个
    select(-rowMean) %>% #反向选择去除rowMean这一列
    tibble::column_to_rownames(colnames(.)[1])
  
  #save(exprSet, file = "exprSet_.Rda")  ## 带gene symbol 表达矩阵
  
  pdata = pData(gset)
  table(pdata$characteristics_ch1.3)
  pdata1 <- subset(pdata, pdata$characteristics_ch1.3 ==  "histology: ADC")
  pdata1 <- pdata1[,c(2,41:50)]
  pdata2 <- subset(pdata, pdata$characteristics_ch1.3 == "histology: SQC")
  pdata2 <- pdata2[,c(2,41:50)]
  str(pdata1)
  str(pdata2)
  
  ## 选取LUAD数据,1为LUAD,2为LUSC(以后可能会用到)
  GSM <- pdata1$geo_accession
  a <- exprSet
  gene <- rownames(exprSet)
  GSM <- pdata1$geo_accession
  expr1 <- exprSet[,colnames(exprSet) %in% GSM]
  exprSet <- a
  GSM <- pdata2$geo_accession
  expr2 <- exprSet[,colnames(exprSet) %in% GSM]
  save(expr1,expr2,pdata1,pdata2, file = "GSE30219_expr_clinic.Rda")
}
