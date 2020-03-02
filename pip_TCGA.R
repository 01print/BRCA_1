###################
# Pancheng Wu
# created in 2020-02-11/
# 03-02/
# TCGA-BRCA analysis 
# PUMC
# XINKAILU HUOTNG 80#, JIANGUOMEN ROAD, DONGCHENG DISTRICT, Beijing, China
#
#
##################
#
load("exprSet_m_vst_17915.Rda")
load("innate_imm.Rda")
imm <- innate_imm$`Gene Symbol`
a <- exprSet_m_vst
ge <- a$gene
rownames(a) <- ge
a <- a[,-1]
a <- as.data.frame(t(a))
a <- a[,colnames(a) %in% imm]  ## 909

## 
TCGA_id <- rownames(a)
table(substring(TCGA_id,14,15)) ## 01/1097 // 11/113
TCGA_id <- TCGA_id[substring(TCGA_id,14,15)== "01"]
a <- a[rownames(a) %in% TCGA_id,] ## 1097

## clinical information 
phen <- read.csv("Pan_cancer_phen.csv")
BRCA <- subset(phen, phen$cancer.type.abbreviation == "BRCA")
BRCA <- BRCA[,c(1:7,9,12,24,26,27)]
BRCA <- BRCA[!duplicated(BRCA$sample),] ## 去除重复值，为0
save(BRCA, file = "BRCA_clinical.Rda")
## import clinical data

a <- expr_vst_imm
TCGA_id <- rownames(a)
id <- substring(TCGA_id,1,15)
a <- cbind("TCGA_id" = id, a)
a$TCGA_id <- as.character(a$TCGA_id)
rownames(a) <- NULL
names(LUAD_phen)[2] <- "TCGA_id"
expr_imm_clinic <- merge(a, LUAD_phen,by = "TCGA_id")  ## 511个样本
expr_imm_clinic <- expr_imm_clinic[!duplicated(expr_imm_clinic$TCGA_id),]  ## 508样本/ 时刻谨记样本名要去除重复值
## 查看并去除存在正常组织样本
TCGA_id1 <- expr_imm_clinic$TCGA_id
table(substring(TCGA_id1,14,15))
TCGA_id <- TCGA_id1[substring(TCGA_id1,14,15) != "11"]   ## 6个正常组织去掉，剩余502个病灶
rownames(expr_imm_clinic) <- TCGA_id1
expr_imm_clinic <- expr_imm_clinic[,-1]
expr_imm_clinic <- expr_imm_clinic[rownames(expr_imm_clinic) %in% TCGA_id,] ## 对这502个样本进行训练和验证

## 对数据进一步清洗，去除生存时间为NA的,去除生存时间小于30天的
a <- expr_imm_clinic
TCGA_id <- rownames(a)
a <- cbind("TCGA_id" = TCGA_id, a)
a$TCGA_id <- as.character(a$TCGA_id)
rownames(a) <- NULL
## 去除生存时间为NA的行.去除9个，剩余493个
a <- a[complete.cases(a$OS.time),] 
## 去除生存时间小于30天的样本。从493到479个
a <- subset(a, a$OS.time > 30) 
rownames(a) <- 1:479
expr_imm_clinic_final <- a
save(expr_imm_clinic_final, file = "expr_imm_clinic_final_new_diff.Rda")



library(readr)
a <- read.csv("Pan_cancer_phen.csv")
colnames(a)
### 提取LUAD的
LUAD <- subset(a,a$cancer.type.abbreviation == "LUAD")
LUAD <- LUAD[,c(1:5,7,14,26,27)]

a <- read.csv("TCGA_479.csv")
set.seed(2)
ind <- sample(2, nrow(a),
              replace = TRUE, prob = c(0.7,0.3)) 
trainset <- a[ind == 1,] ## trainset 335个
validationset <- a[ind == 2,] ## 144个

a <- trainset ## 备份
a$TCGA_id <- as.character(a$TCGA_id)
library(survival)
library(plyr)
Basurv <- Surv(time = a$OS.time, event = a$os)
a$Basurv <- with(a, Basurv)
## use function
Ucox <- function(x){
  FML <- as.formula(paste0('Basurv~',x))
  GCox <- coxph(FML,data = a)
  GSum <- summary(GCox)
  CI <- paste0(round(GSum$conf.int[,3:4],3),collapse = '-')
  HR <- round(GSum$coefficients[,2],3)
  Pvalue <- round(GSum$coefficients[,5],3)
  Unicox <- data.frame('gene_id' = x,
                       'HR' = HR,
                       '95%CI' = CI,
                       'p' = Pvalue)
  return(Unicox)
}
varNames <- colnames(a)[2:911]
Univar <- lapply(varNames, Ucox)
Univar <- ldply(Univar,data.frame)
Univar$gene_id <- as.character(Univar$gene_id)
Univar$X95.CI <- as.character(Univar$X95.CI)
# save(Univar, file = "Univar_pall.Rda") ## 保存p值结果
table(Univar$p < 0.05) ## 212个
table(Univar$p < 0.01) ## 82个
table(Univar$p < 0.005) ## 46个
table(Univar$p < 0.001) ## 16个
Univar_0.005 <- subset(Univar,Univar$p < 0.005)
gene <- Univar_0.005$gene_id
TCGA_P_0.05_x1 <- subset(Univar_0.05, Univar_0.05$HR < 1)
TCGA_P_0.05_x1 <- TCGA_P_0.05_x1$gene_id ## 129个

## trainset analysis 
trainset{
  a <- trainset
  OS <- a[,c(1,912:926)]
  a <- a[,colnames(a) %in% gene]
  a <- cbind(OS, a)
  a$X_PATIENT <- as.character(a$X_PATIENT)
  a$cancer.type.abbreviation <- as.character(a$cancer.type.abbreviation)
  
  ## lasso 回归分析
  ## 制作x1,y1文件，进行lasso回归分析
  TCGA_id <- a$TCGA_id
  x1 <- a[,c(17:62)]
  x1 <- as.matrix(x1)
  y1 <- a[,c(9,8)]
  names(y1) <- c("time","status")
  y1 <- as.matrix(y1)
  ## 进行分析
  library(glmnet)
  fit <- glmnet(x1, y1, family = "cox")
  plot(fit)
  cvfit = cv.glmnet(x1, y1, family = "cox")
  plot(cvfit)
  
  coef.min = coef(cvfit, s = "lambda.min")
  active.min = which(coef.min != 0)
  index.min = coef.min[active.min]
  index.min
  coef.min
  row.names(coef.min)[active.min]  
  
  ## 系数和所筛选到的26个基因
  bt <- index.min
  gene_14 <- row.names(coef.min)[active.min]
  save(bt,gene_14, file = "gene_coeffi_14.Rda")
  ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
  a <- trainset
  a <- a[,colnames(a) %in% gene_14]
  a <- cbind(OS,a)
  
  ## 多因素回归，将筛选的17基因进一步缩小
  Multi{
    library(survival)
    library(survminer)
    b <- a
    Basurv <- Surv(time = a$OS.time, event = a$os)
    fml <-as.formula(paste0('Basurv~',paste0(colnames(a)[17:44],collapse = '+')))
    MultiCox <- coxph(fml, data = a,ties = "breslow")
    MultiSum <- summary(MultiCox)
    MultiName <- as.character(colnames(a)[17:44])
    MHR <- round(MultiSum$coefficients[,2],3)
    MPV <- round(MultiSum$coefficients[,5],3)
    MCIL <- round(MultiSum$conf.int[,3],3)
    MCIU <- round(MultiSum$conf.int[,4],3)
    MCI <- paste0(MCIL,'-',MCIU)
    Mulcox <- data.frame('Chracteristics' = MultiName,
                         'HR' = MHR,
                         '95%CI' = MCI,
                         'P' = MPV)
    
    table(Mulcox$P < 0.05) ## 5个基因
    MUlcox_0.05 <- subset(Mulcox, Mulcox$P < 0.05)
    gene_5 <- MUlcox_0.05$Chracteristics
    gene_5 <- as.character(gene_5)
    coeff_log <- log(MUlcox_0.05$HR)
    coeff_lasso <- c(-0.347794879,0.050660468,0.195827737,-0.303941288,0.364098876)
    save(gene_5, coeff_log,file = "gene_coeff_212_1_1000.Rda")
    save(bt,gene_17, file = "gene_coeff_212_L1_1000.Rda")
  }
  ## 用多因素COX回归结果来构建模型
  
  ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
  a$TCGA_id <- as.character(a$TCGA_id)
  a$X_PATIENT <- as.character(a$X_PATIENT)
  a$cancer.type.abbreviation <- a$cancer.type.abbreviation
  
  a$riskscore <- NA
  for (i in 1:335){
    a[i,31] <- bt[1]*a[i,17]+bt[2]*a[i,18]+bt[3]*a[i,19]+bt[4]*a[i,20]+bt[5]*a[i,21]+bt[6]*a[i,22]+bt[7]*a[i,23]+bt[8]*a[i,24]+
      bt[9]*a[i,25]+bt[10]*a[i,26]+bt[11]*a[i,27]+bt[12]*a[i,28]+bt[13]*a[i,29]+bt[14]*a[i,30]
  }
  save(a, file = "train_05_0.005.Rda") ## 单因素cox--lasso--多因素cox
  #+bt[15]*a[i,32]+bt[16]*a[i,33]+bt[17]*a[i,34]+
  bt[18]*a[i,35]+bt[19]*a[i,36]+bt[20]*a[i,36]+bt[21]*a[i,37]+bt[22]*a[i,38]+bt[23]*a[i,39]+
    bt[24]*a[i,40]+bt[25]*a[i,41]+bt[26]*a[i,42]+bt[27]*a[i,43]+bt[28]*a[i,44]
  ##生存分析
  library(survival)
  library(survminer)
  a$OS.time1 <- a$OS.time/30
  group <- ifelse(a$riskscore>median(a$riskscore),'High risk','Low risk')
  table(group)
  sfit <- survfit(Surv(OS.time1, OS)~group, data=a)  ## 构建生存对象，进行数据处理
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE,
             legend.title = "Trainset",
             xlab = "Time(months)",
             risk.table = TRUE,
             tables.height = 0.2,
             tables.theme = theme_cleantable(),
             # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
             # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
             palette = c("#E7B800", "#2E9FDF"),
             ggtheme = theme_bw()
  )
  ##
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  a$OS.time1 <- a$OS.time/30
  a$status <- a$os
  ## Divide the paitents into high and low risk group by the surv.cut fuction.
  ## 按照函数进行分组
  new_function{
    surv.cut <- surv_cutpoint(
      a,
      time = "OS.time1",
      event = "status",
      variables = "riskscore"
    )
    summary(surv.cut)
    plot(surv.cut, "riskscore", palette = "npg")
    
    
    surv.cat <- surv_categorize(surv.cut)
    
    surv.fit <- survfit(Surv(OS.time1, status) ~ riskscore,
                        data = surv.cat)
    ggsurvplot(
      surv.fit,                     # survfit object with calculated statistics.
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      # show confidence intervals for 
      # point estimaes of survival curves.
      #xlim = c(0,3000),        # present narrower X axis, but not affect
      # survival estimates.
      break.time.by = 24,    # break X axis in time intervals by 500. 
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = FALSE # show bars instead of names in text annotations
      # in legend of risk table
    )
    ggsurvplot(surv.fit, conf.int=F, pval=TRUE,
               legend.title = "Training set",
               xlab = "Time(months)",
               legend.labs = c("High risk", "Low risk"),
               risk.table = TRUE,
               break.time.by = 24,
               xlim = c(0,168), ## 2020-02-09学到，自定义X轴
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw()
    )
  }
  ## 绘制多个时间点ROC曲线，评估模型效能
  library("timeROC")
  a$OS.time2 <- a$OS.time/365
  ROC_a <- timeROC(T = a$OS.time2,delta=a$status,
                   marker=a$riskscore,cause=1,
                   weighting="marginal",
                   times=c(1,3,5),ROC=TRUE)
  ROC_a
  plot(ROC_a,time=1,title=FALSE,lwd=2)
  plot(ROC_a,time=3,col="blue",add=TRUE,title=FALSE,lwd=2)
  plot(ROC_a,time=5,col="black",add=TRUE,title=FALSE,lwd=2)
  legend("bottomright",
         c(paste0("AUC at 1 years: ",round(ROC_a$AUC[1],2)),
           paste0("AUC at 3 years: ",round(ROC_a$AUC[2],2)),
           paste0("AUC at 5 years: ",round(ROC_a$AUC[3],2))),
         col=c("red","blue","black"),lwd=2,bty = "n")
  
  ##风险因子分布图绘制，包括3部分：riskscore散点图、生存状态散点图、和热图
  ##library(ComplexHeatmap) 以后找时间再学习，现在没时间学这个包了，花了好几天，没达到想要的结果
  ## 绘制散点图
  library('ggplot2')
  load("train_17_riskscore.Rda")
  ID <- 1:355
  b <- a
  ## 按照riskscore大小进行排列,从小往大进行排列,方便搜索
  b <- b[order(b$riskscore,decreasing = FALSE),]
  b$OS <- as.factor(b$OS)
  TCGA_id <- b$TCGA_id
  risk <- b$riskscore
  risk <- cbind("ID" = ID, risk)
  colnames(risk)[2] <- "Risk score"
  risk <- as.data.frame(risk)
  colnames(risk)[1] <- "Patients (increasing risk score)"
  risk$Group <- ifelse(risk$`Risk score`>median(risk$`Risk score`),'High risk','Low risk')
  ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
  p1 <- ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
  
  ## 生存状态散点图
  OS <- b[,c(2,4)]
  OS$Status <- ifelse(OS$OS == 0, 'Alive','Dead')
  OS <- cbind("Patients (increasing risk score)" = ID,OS)
  OS$time2 <- OS$OS.time/365
  colnames(OS)[5] <- "Survival time (years)" 
  ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
    geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
  p2 <- ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
    geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
  
  ## 热图绘制
  save(heat,OS,risk, file = "train_17_3_picture.Rda")
  
  ## 筛选基因热图绘制，用pheatmap作图
  library(pheatmap)
  library(cowplot)
  Heat <- b[,c(5:21)]
  rownames(Heat) <- ID
  Heat <- t(Heat)
  Risk <- risk[,3]
  Risk <- as.data.frame(Risk)
  Risk$Risk <- ifelse(Risk$Risk == "high risk", "High", "Low") 
  rownames(Risk) <- c(1:355) 
  ###自定义颜色
  ann_colors = list(
    Risk = c(Low = "#00BFC4", High = "#F8766D"))
  pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
           cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
  p3 <- pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
                 cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
  save(risk, OS, Heat,Risk, file = "Train_17_3figure.Rda")
  ## 整合, 不能很好的将3幅图合在一起。目前，只能这样了。
  plot_grid(p1, p2,
            labels = c("A", "B"),
            align = 'v',ncol = 1)
}

## Internal validation 
validationset{
  ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
  a <- validationset
  OS <- a[,c(1,912:926)]
  a <- a[,colnames(a) %in% gene_14]
  a <- cbind(OS,a)
  
  ## 多因素回归，将筛选的17基因进一步缩小
  Multi{
    library(survival)
    library(survminer)
    b <- a
    Basurv <- Surv(time = a$OS.time, event = a$os)
    fml <-as.formula(paste0('Basurv~',paste0(colnames(a)[17:44],collapse = '+')))
    MultiCox <- coxph(fml, data = a,ties = "breslow")
    MultiSum <- summary(MultiCox)
    MultiName <- as.character(colnames(a)[17:44])
    MHR <- round(MultiSum$coefficients[,2],3)
    MPV <- round(MultiSum$coefficients[,5],3)
    MCIL <- round(MultiSum$conf.int[,3],3)
    MCIU <- round(MultiSum$conf.int[,4],3)
    MCI <- paste0(MCIL,'-',MCIU)
    Mulcox <- data.frame('Chracteristics' = MultiName,
                         'HR' = MHR,
                         '95%CI' = MCI,
                         'P' = MPV)
    
    table(Mulcox$P < 0.05) ## 5个基因
    MUlcox_0.05 <- subset(Mulcox, Mulcox$P < 0.05)
    gene_5 <- MUlcox_0.05$Chracteristics
    gene_5 <- as.character(gene_5)
    coeff_log <- log(MUlcox_0.05$HR)
    coeff_lasso <- c(-0.347794879,0.050660468,0.195827737,-0.303941288,0.364098876)
    save(gene_5, coeff_log,file = "gene_coeff_212_1_1000.Rda")
    save(bt,gene_17, file = "gene_coeff_212_L1_1000.Rda")
  }
  ## 用多因素COX回归结果来构建模型
  
  ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
  a$TCGA_id <- as.character(a$TCGA_id)
  a$X_PATIENT <- as.character(a$X_PATIENT)
  a$cancer.type.abbreviation <- a$cancer.type.abbreviation
  
  a$riskscore <- NA
  for (i in 1:144){
    a[i,31] <- bt[1]*a[i,17]+bt[2]*a[i,18]+bt[3]*a[i,19]+bt[4]*a[i,20]+bt[5]*a[i,21]+bt[6]*a[i,22]+bt[7]*a[i,23]+bt[8]*a[i,24]+
      bt[9]*a[i,25]+bt[10]*a[i,26]+bt[11]*a[i,27]+bt[12]*a[i,28]+bt[13]*a[i,29]+bt[14]*a[i,30]
  }
  
  #+bt[15]*a[i,32]+bt[16]*a[i,33]+bt[17]*a[i,34]+
  bt[18]*a[i,35]+bt[19]*a[i,36]+bt[20]*a[i,36]+bt[21]*a[i,37]+bt[22]*a[i,38]+bt[23]*a[i,39]+
    bt[24]*a[i,40]+bt[25]*a[i,41]+bt[26]*a[i,42]+bt[27]*a[i,43]+bt[28]*a[i,44]
  ##生存分析
  library(survival)
  library(survminer)
  a$OS.time1 <- a$OS.time/30
  group <- ifelse(a$riskscore>median(a$riskscore),'High risk','Low risk')
  table(group)
  sfit <- survfit(Surv(OS.time1, OS)~group, data=a)  ## 构建生存对象，进行数据处理
  sfit
  summary(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE,
             legend.title = "Trainset",
             xlab = "Time(months)",
             risk.table = TRUE,
             tables.height = 0.2,
             tables.theme = theme_cleantable(),
             # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
             # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
             palette = c("#E7B800", "#2E9FDF"),
             ggtheme = theme_bw()
  )
  ##
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
  a$OS.time1 <- a$OS.time/30
  a$status <- a$os
  ## Divide the paitents into high and low risk group by the surv.cut fuction.
  ## 按照函数进行分组
  new_function{
    surv.cut <- surv_cutpoint(
      a,
      time = "OS.time1",
      event = "status",
      variables = "riskscore"
    )
    summary(surv.cut)
    plot(surv.cut, "riskscore", palette = "npg")
    
    
    surv.cat <- surv_categorize(surv.cut)
    
    surv.fit <- survfit(Surv(OS.time1, status) ~ riskscore,
                        data = surv.cat)
    ggsurvplot(
      surv.fit,                     # survfit object with calculated statistics.
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      # show confidence intervals for 
      # point estimaes of survival curves.
      #xlim = c(0,3000),        # present narrower X axis, but not affect
      # survival estimates.
      break.time.by = 24,    # break X axis in time intervals by 500. 
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = FALSE # show bars instead of names in text annotations
      # in legend of risk table
    )
    ggsurvplot(surv.fit, conf.int=F, pval=TRUE,
               legend.title = "Training set",
               xlab = "Time(months)",
               legend.labs = c("High risk", "Low risk"),
               risk.table = TRUE,
               break.time.by = 24,
               xlim = c(0,96), ## 2020-02-09学到，自定义X轴
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw()
    )
  }
  ## 绘制多个时间点ROC曲线，评估模型效能
  library("timeROC")
  a$OS.time2 <- a$OS.time/365
  ROC_a <- timeROC(T = a$OS.time2,delta=a$status,
                   marker=a$riskscore,cause=1,
                   weighting="marginal",
                   times=c(1,3,5),ROC=TRUE)
  ROC_a
  plot(ROC_a,time=1,title=FALSE,lwd=2)
  plot(ROC_a,time=3,col="blue",add=TRUE,title=FALSE,lwd=2)
  plot(ROC_a,time=5,col="black",add=TRUE,title=FALSE,lwd=2)
  legend("bottomright",
         c(paste0("AUC at 1 years: ",round(ROC_a$AUC[1],2)),
           paste0("AUC at 3 years: ",round(ROC_a$AUC[2],2)),
           paste0("AUC at 5 years: ",round(ROC_a$AUC[3],2))),
         col=c("red","blue","black"),lwd=2,bty = "n")
  
  ##风险因子分布图绘制，包括3部分：riskscore散点图、生存状态散点图、和热图
  ##library(ComplexHeatmap) 以后找时间再学习，现在没时间学这个包了，花了好几天，没达到想要的结果
  ## 绘制散点图
  library('ggplot2')
  load("train_17_riskscore.Rda")
  ID <- 1:355
  b <- a
  ## 按照riskscore大小进行排列,从小往大进行排列,方便搜索
  b <- b[order(b$riskscore,decreasing = FALSE),]
  b$OS <- as.factor(b$OS)
  TCGA_id <- b$TCGA_id
  risk <- b$riskscore
  risk <- cbind("ID" = ID, risk)
  colnames(risk)[2] <- "Risk score"
  risk <- as.data.frame(risk)
  colnames(risk)[1] <- "Patients (increasing risk score)"
  risk$Group <- ifelse(risk$`Risk score`>median(risk$`Risk score`),'High risk','Low risk')
  ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
  p1 <- ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
  
  ## 生存状态散点图
  OS <- b[,c(2,4)]
  OS$Status <- ifelse(OS$OS == 0, 'Alive','Dead')
  OS <- cbind("Patients (increasing risk score)" = ID,OS)
  OS$time2 <- OS$OS.time/365
  colnames(OS)[5] <- "Survival time (years)" 
  ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
    geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
  p2 <- ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
    geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
  
  ## 热图绘制
  save(heat,OS,risk, file = "train_17_3_picture.Rda")
  
  ## 筛选基因热图绘制，用pheatmap作图
  library(pheatmap)
  library(cowplot)
  Heat <- b[,c(5:21)]
  rownames(Heat) <- ID
  Heat <- t(Heat)
  Risk <- risk[,3]
  Risk <- as.data.frame(Risk)
  Risk$Risk <- ifelse(Risk$Risk == "high risk", "High", "Low") 
  rownames(Risk) <- c(1:355) 
  ###自定义颜色
  ann_colors = list(
    Risk = c(Low = "#00BFC4", High = "#F8766D"))
  pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
           cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
  p3 <- pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
                 cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
  save(risk, OS, Heat,Risk, file = "Train_17_3figure.Rda")
  ## 整合, 不能很好的将3幅图合在一起。目前，只能这样了。
  plot_grid(p1, p2,
            labels = c("A", "B"),
            align = 'v',ncol = 1)
}



## 0.5,0.5
TCGA_half{
  preprocess{
    a <- read.csv("TCGA_479.csv")
    set.seed(2)
    tmp <- sample(1:479)
    id1 <- tmp[1:240]
    id2 <- tmp[241:479]
    trainset <- a[rownames(a) %in% id1,]
    validationset <- a[rownames(a) %in% id2,]
    a <- trainset ## 备份
    a$TCGA_id <- as.character(a$TCGA_id)
    library(survival)
    library(plyr)
    Basurv <- Surv(time = a$OS.time, event = a$os)
    a$Basurv <- with(a, Basurv)
    ## use function
    Ucox <- function(x){
      FML <- as.formula(paste0('Basurv~',x))
      GCox <- coxph(FML,data = a)
      GSum <- summary(GCox)
      CI <- paste0(round(GSum$conf.int[,3:4],3),collapse = '-')
      HR <- round(GSum$coefficients[,2],3)
      Pvalue <- round(GSum$coefficients[,5],3)
      Unicox <- data.frame('gene_id' = x,
                           'HR' = HR,
                           '95%CI' = CI,
                           'p' = Pvalue)
      return(Unicox)
    }
    varNames <- colnames(a)[2:911]
    Univar <- lapply(varNames, Ucox)
    Univar <- ldply(Univar,data.frame)
    Univar$gene_id <- as.character(Univar$gene_id)
    Univar$X95.CI <- as.character(Univar$X95.CI)
    # save(Univar, file = "Univar_pall.Rda") ## 保存p值结果
    table(Univar$p < 0.05) ## 121个
    table(Univar$p < 0.01) ## 25个
    table(Univar$p < 0.005) ## 18个
    table(Univar$p < 0.001) ## 0个
    Univar_0.01 <- subset(Univar,Univar$p < 0.01)
    Univar_0.005 <- subset(Univar, Univar$p < 0.005)
    Univar_0.05 <- subset(Univar, Univar$p < 0.05)
    gene <- Univar_0.01$gene_id
    gene <- Univar_0.005$gene_id
    gene <- Univar_0.05$gene_id
    
  }
  
  ## trainset analysis 
  trainset{
    a <- trainset
    OS <- a[,c(1,912:926)]
    a <- a[,colnames(a) %in% gene]
    a <- cbind(OS, a)
    a$X_PATIENT <- as.character(a$X_PATIENT)
    a$cancer.type.abbreviation <- as.character(a$cancer.type.abbreviation)
    
    ## lasso 回归分析
    ## 制作x1,y1文件，进行lasso回归分析
    TCGA_id <- a$TCGA_id
    x1 <- a[,c(17:34)]
    x1 <- as.matrix(x1)
    y1 <- a[,c(9,8)]
    names(y1) <- c("time","status")
    y1 <- as.matrix(y1)
    ## 进行分析
    library(glmnet)
    fit <- glmnet(x1, y1, family = "cox")
    plot(fit)
    cvfit = cv.glmnet(x1, y1, family = "cox", maxit = 1000)
    plot(cvfit)
    
    coef.min = coef(cvfit, s = "lambda.min")
    active.min = which(coef.min != 0)
    index.min = coef.min[active.min]
    index.min
    coef.min
    row.names(coef.min)[active.min]  
    
    ## 系数和所筛选到的26个基因
    bt <- index.min
    gene_12 <- row.names(coef.min)[active.min]
    save(bt,gene_12, file = "gene_coeffi_12.Rda")
    ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
    a <- trainset
    a <- a[,colnames(a) %in% gene_12]
    a <- cbind(OS,a)
    
    ## 多因素回归，将筛选的17基因进一步缩小
    Multi{
      library(survival)
      library(survminer)
      b <- a
      Basurv <- Surv(time = a$OS.time, event = a$os)
      fml <-as.formula(paste0('Basurv~',paste0(colnames(a)[17:44],collapse = '+')))
      MultiCox <- coxph(fml, data = a,ties = "breslow")
      MultiSum <- summary(MultiCox)
      MultiName <- as.character(colnames(a)[17:44])
      MHR <- round(MultiSum$coefficients[,2],3)
      MPV <- round(MultiSum$coefficients[,5],3)
      MCIL <- round(MultiSum$conf.int[,3],3)
      MCIU <- round(MultiSum$conf.int[,4],3)
      MCI <- paste0(MCIL,'-',MCIU)
      Mulcox <- data.frame('Chracteristics' = MultiName,
                           'HR' = MHR,
                           '95%CI' = MCI,
                           'P' = MPV)
      
      table(Mulcox$P < 0.05) ## 5个基因
      MUlcox_0.05 <- subset(Mulcox, Mulcox$P < 0.05)
      gene_5 <- MUlcox_0.05$Chracteristics
      gene_5 <- as.character(gene_5)
      coeff_log <- log(MUlcox_0.05$HR)
      coeff_lasso <- c(-0.347794879,0.050660468,0.195827737,-0.303941288,0.364098876)
      save(gene_5, coeff_log,file = "gene_coeff_212_1_1000.Rda")
      save(bt,gene_17, file = "gene_coeff_212_L1_1000.Rda")
    }
    ## 用多因素COX回归结果来构建模型
    
    ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
    a$TCGA_id <- as.character(a$TCGA_id)
    a$X_PATIENT <- as.character(a$X_PATIENT)
    a$cancer.type.abbreviation <- a$cancer.type.abbreviation
    
    a$riskscore <- NA
    for (i in 1:240){
      a[i,29] <- bt[1]*a[i,17]+bt[2]*a[i,18]+bt[3]*a[i,19]+bt[4]*a[i,20]+bt[5]*a[i,21]+bt[6]*a[i,22]+bt[7]*a[i,23]+bt[8]*a[i,24]+
        bt[9]*a[i,25]+bt[10]*a[i,26]+bt[11]*a[i,27]+bt[12]*a[i,28]
    }
    ##生存分析
    library(survival)
    library(survminer)
    a$OS.time1 <- a$OS.time/30
    group <- ifelse(a$riskscore>median(a$riskscore),'High risk','Low risk')
    table(group)
    sfit <- survfit(Surv(OS.time1, OS)~group, data=a)  ## 构建生存对象，进行数据处理
    sfit
    summary(sfit)
    ggsurvplot(sfit, conf.int=F, pval=TRUE,
               legend.title = "Trainset",
               xlab = "Time(months)",
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw()
    )
    ##
    library(dplyr)
    library(survival)
    library(survminer)
    library(ggplot2)
    a$OS.time1 <- a$OS.time/30
    a$status <- a$os
    ## Divide the paitents into high and low risk group by the surv.cut fuction.
    ## 按照函数进行分组
    new_function{
      surv.cut <- surv_cutpoint(
        a,
        time = "OS.time1",
        event = "status",
        variables = "riskscore"
      )
      summary(surv.cut)
      plot(surv.cut, "riskscore", palette = "npg")
      
      
      surv.cat <- surv_categorize(surv.cut)
      
      surv.fit <- survfit(Surv(OS.time1, status) ~ riskscore,
                          data = surv.cat)
      ggsurvplot(
        surv.fit,                     # survfit object with calculated statistics.
        risk.table = TRUE,       # show risk table.
        pval = TRUE,             # show p-value of log-rank test.
        # show confidence intervals for 
        # point estimaes of survival curves.
        #xlim = c(0,3000),        # present narrower X axis, but not affect
        # survival estimates.
        break.time.by = 24,    # break X axis in time intervals by 500. 
        risk.table.y.text.col = T, # colour risk table text annotations.
        risk.table.y.text = FALSE # show bars instead of names in text annotations
        # in legend of risk table
      )
      ggsurvplot(surv.fit, conf.int=F, pval=TRUE,
                 legend.title = "Training set",
                 xlab = "Time(months)",
                 legend.labs = c("High risk", "Low risk"),
                 risk.table = TRUE,
                 break.time.by = 24,
                 xlim = c(0,96), ## 2020-02-09学到，自定义X轴
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                 # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                 palette = c("#E7B800", "#2E9FDF"),
                 ggtheme = theme_bw()
      )
    }
    ## 绘制多个时间点ROC曲线，评估模型效能
    library("timeROC")
    a$OS.time2 <- a$OS.time/365
    ROC_a <- timeROC(T = a$OS.time2,delta=a$status,
                     marker=a$riskscore,cause=1,
                     weighting="marginal",
                     times=c(1,3,5),ROC=TRUE)
    ROC_a
    plot(ROC_a,time=1,title=FALSE,lwd=2)
    plot(ROC_a,time=3,col="blue",add=TRUE,title=FALSE,lwd=2)
    plot(ROC_a,time=5,col="black",add=TRUE,title=FALSE,lwd=2)
    legend("bottomright",
           c(paste0("AUC at 1 years: ",round(ROC_a$AUC[1],2)),
             paste0("AUC at 3 years: ",round(ROC_a$AUC[2],2)),
             paste0("AUC at 5 years: ",round(ROC_a$AUC[3],2))),
           col=c("red","blue","black"),lwd=2,bty = "n")
    
    ##风险因子分布图绘制，包括3部分：riskscore散点图、生存状态散点图、和热图
    ##library(ComplexHeatmap) 以后找时间再学习，现在没时间学这个包了，花了好几天，没达到想要的结果
    ## 绘制散点图
    library('ggplot2')
    load("train_17_riskscore.Rda")
    ID <- 1:355
    b <- a
    ## 按照riskscore大小进行排列,从小往大进行排列,方便搜索
    b <- b[order(b$riskscore,decreasing = FALSE),]
    b$OS <- as.factor(b$OS)
    TCGA_id <- b$TCGA_id
    risk <- b$riskscore
    risk <- cbind("ID" = ID, risk)
    colnames(risk)[2] <- "Risk score"
    risk <- as.data.frame(risk)
    colnames(risk)[1] <- "Patients (increasing risk score)"
    risk$Group <- ifelse(risk$`Risk score`>median(risk$`Risk score`),'High risk','Low risk')
    ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
    p1 <- ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
    
    ## 生存状态散点图
    OS <- b[,c(2,4)]
    OS$Status <- ifelse(OS$OS == 0, 'Alive','Dead')
    OS <- cbind("Patients (increasing risk score)" = ID,OS)
    OS$time2 <- OS$OS.time/365
    colnames(OS)[5] <- "Survival time (years)" 
    ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
      geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
    p2 <- ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
      geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
    
    ## 热图绘制
    save(heat,OS,risk, file = "train_17_3_picture.Rda")
    
    ## 筛选基因热图绘制，用pheatmap作图
    library(pheatmap)
    library(cowplot)
    Heat <- b[,c(5:21)]
    rownames(Heat) <- ID
    Heat <- t(Heat)
    Risk <- risk[,3]
    Risk <- as.data.frame(Risk)
    Risk$Risk <- ifelse(Risk$Risk == "high risk", "High", "Low") 
    rownames(Risk) <- c(1:355) 
    ###自定义颜色
    ann_colors = list(
      Risk = c(Low = "#00BFC4", High = "#F8766D"))
    pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
             cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
    p3 <- pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
                   cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
    save(risk, OS, Heat,Risk, file = "Train_17_3figure.Rda")
    ## 整合, 不能很好的将3幅图合在一起。目前，只能这样了。
    plot_grid(p1, p2,
              labels = c("A", "B"),
              align = 'v',ncol = 1)
  }
  
  ## Internal validation 
  validationset{
    ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
    a <- validationset
    OS <- a[,c(1,912:926)]
    a <- a[,colnames(a) %in% gene_12]
    a <- cbind(OS,a)
    
    ## 多因素回归，将筛选的17基因进一步缩小
    Multi{
      library(survival)
      library(survminer)
      b <- a
      Basurv <- Surv(time = a$OS.time, event = a$os)
      fml <-as.formula(paste0('Basurv~',paste0(colnames(a)[17:44],collapse = '+')))
      MultiCox <- coxph(fml, data = a,ties = "breslow")
      MultiSum <- summary(MultiCox)
      MultiName <- as.character(colnames(a)[17:44])
      MHR <- round(MultiSum$coefficients[,2],3)
      MPV <- round(MultiSum$coefficients[,5],3)
      MCIL <- round(MultiSum$conf.int[,3],3)
      MCIU <- round(MultiSum$conf.int[,4],3)
      MCI <- paste0(MCIL,'-',MCIU)
      Mulcox <- data.frame('Chracteristics' = MultiName,
                           'HR' = MHR,
                           '95%CI' = MCI,
                           'P' = MPV)
      
      table(Mulcox$P < 0.05) ## 5个基因
      MUlcox_0.05 <- subset(Mulcox, Mulcox$P < 0.05)
      gene_5 <- MUlcox_0.05$Chracteristics
      gene_5 <- as.character(gene_5)
      coeff_log <- log(MUlcox_0.05$HR)
      coeff_lasso <- c(-0.347794879,0.050660468,0.195827737,-0.303941288,0.364098876)
      save(gene_5, coeff_log,file = "gene_coeff_212_1_1000.Rda")
      save(bt,gene_17, file = "gene_coeff_212_L1_1000.Rda")
    }
    ## 用多因素COX回归结果来构建模型
    
    ## 将筛选到的17个基因，重新匹配训练集，计算风险评分
    a$TCGA_id <- as.character(a$TCGA_id)
    a$X_PATIENT <- as.character(a$X_PATIENT)
    a$cancer.type.abbreviation <- a$cancer.type.abbreviation
    
    a$riskscore <- NA
    for (i in 1:239){
      a[i,29] <- bt[1]*a[i,17]+bt[2]*a[i,18]+bt[3]*a[i,19]+bt[4]*a[i,20]+bt[5]*a[i,21]+bt[6]*a[i,22]+bt[7]*a[i,23]+bt[8]*a[i,24]+
        bt[9]*a[i,25]+bt[10]*a[i,26]+bt[11]*a[i,27]+bt[12]*a[i,28]
    }
    
    ##生存分析
    library(survival)
    library(survminer)
    a$OS.time1 <- a$OS.time/30
    group <- ifelse(a$riskscore>median(a$riskscore),'High risk','Low risk')
    table(group)
    sfit <- survfit(Surv(OS.time1, OS)~group, data=a)  ## 构建生存对象，进行数据处理
    sfit
    summary(sfit)
    ggsurvplot(sfit, conf.int=F, pval=TRUE,
               legend.title = "Trainset",
               xlab = "Time(months)",
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw()
    )
    ##
    library(dplyr)
    library(survival)
    library(survminer)
    library(ggplot2)
    a$OS.time1 <- a$OS.time/30
    a$status <- a$os
    ## Divide the paitents into high and low risk group by the surv.cut fuction.
    ## 按照函数进行分组
    new_function{
      surv.cut <- surv_cutpoint(
        a,
        time = "OS.time1",
        event = "status",
        variables = "riskscore"
      )
      summary(surv.cut)
      plot(surv.cut, "riskscore", palette = "npg")
      
      
      surv.cat <- surv_categorize(surv.cut)
      
      surv.fit <- survfit(Surv(OS.time1, status) ~ riskscore,
                          data = surv.cat)
      ggsurvplot(
        surv.fit,                     # survfit object with calculated statistics.
        risk.table = TRUE,       # show risk table.
        pval = TRUE,             # show p-value of log-rank test.
        # show confidence intervals for 
        # point estimaes of survival curves.
        #xlim = c(0,3000),        # present narrower X axis, but not affect
        # survival estimates.
        break.time.by = 24,    # break X axis in time intervals by 500. 
        risk.table.y.text.col = T, # colour risk table text annotations.
        risk.table.y.text = FALSE # show bars instead of names in text annotations
        # in legend of risk table
      )
      ggsurvplot(surv.fit, conf.int=F, pval=TRUE,
                 legend.title = "validation set",
                 xlab = "Time(months)",
                 legend.labs = c("High risk", "Low risk"),
                 risk.table = TRUE,
                 break.time.by = 24,
                 xlim = c(0,96), ## 2020-02-09学到，自定义X轴
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                 # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                 palette = c("#E7B800", "#2E9FDF"),
                 ggtheme = theme_bw()
      )
    }
    ## 绘制多个时间点ROC曲线，评估模型效能
    library("timeROC")
    a$OS.time2 <- a$OS.time/365
    ROC_a <- timeROC(T = a$OS.time2,delta=a$status,
                     marker=a$riskscore,cause=1,
                     weighting="marginal",
                     times=c(1,3,5),ROC=TRUE)
    ROC_a
    plot(ROC_a,time=1,title=FALSE,lwd=2)
    plot(ROC_a,time=3,col="blue",add=TRUE,title=FALSE,lwd=2)
    plot(ROC_a,time=5,col="black",add=TRUE,title=FALSE,lwd=2)
    legend("bottomright",
           c(paste0("AUC at 1 years: ",round(ROC_a$AUC[1],2)),
             paste0("AUC at 3 years: ",round(ROC_a$AUC[2],2)),
             paste0("AUC at 5 years: ",round(ROC_a$AUC[3],2))),
           col=c("red","blue","black"),lwd=2,bty = "n")
    
    ##风险因子分布图绘制，包括3部分：riskscore散点图、生存状态散点图、和热图
    ##library(ComplexHeatmap) 以后找时间再学习，现在没时间学这个包了，花了好几天，没达到想要的结果
    ## 绘制散点图
    library('ggplot2')
    load("train_17_riskscore.Rda")
    ID <- 1:355
    b <- a
    ## 按照riskscore大小进行排列,从小往大进行排列,方便搜索
    b <- b[order(b$riskscore,decreasing = FALSE),]
    b$OS <- as.factor(b$OS)
    TCGA_id <- b$TCGA_id
    risk <- b$riskscore
    risk <- cbind("ID" = ID, risk)
    colnames(risk)[2] <- "Risk score"
    risk <- as.data.frame(risk)
    colnames(risk)[1] <- "Patients (increasing risk score)"
    risk$Group <- ifelse(risk$`Risk score`>median(risk$`Risk score`),'High risk','Low risk')
    ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
    p1 <- ggplot(risk,aes(x = `Patients (increasing risk score)`, y = `Risk score`, color = Group)) + geom_point()
    
    ## 生存状态散点图
    OS <- b[,c(2,4)]
    OS$Status <- ifelse(OS$OS == 0, 'Alive','Dead')
    OS <- cbind("Patients (increasing risk score)" = ID,OS)
    OS$time2 <- OS$OS.time/365
    colnames(OS)[5] <- "Survival time (years)" 
    ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
      geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
    p2 <- ggplot(OS, aes(x = `Patients (increasing risk score)`, y = `Survival time (years)`, color = Status))+
      geom_point(size = 1.5)+scale_color_manual(values =c('#00BFC4','#F8766D'))
    
    ## 热图绘制
    save(heat,OS,risk, file = "train_17_3_picture.Rda")
    
    ## 筛选基因热图绘制，用pheatmap作图
    library(pheatmap)
    library(cowplot)
    Heat <- b[,c(5:21)]
    rownames(Heat) <- ID
    Heat <- t(Heat)
    Risk <- risk[,3]
    Risk <- as.data.frame(Risk)
    Risk$Risk <- ifelse(Risk$Risk == "high risk", "High", "Low") 
    rownames(Risk) <- c(1:355) 
    ###自定义颜色
    ann_colors = list(
      Risk = c(Low = "#00BFC4", High = "#F8766D"))
    pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
             cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
    p3 <- pheatmap(Heat,annotation_col = Risk,annotation_legend = FALSE,annotation_colors = ann_colors,
                   cluster_cols = FALSE,show_colnames = F,cluster_rows = FALSE)
    save(risk, OS, Heat,Risk, file = "Train_17_3figure.Rda")
    ## 整合, 不能很好的将3幅图合在一起。目前，只能这样了。
    plot_grid(p1, p2,
              labels = c("A", "B"),
              align = 'v',ncol = 1)
  }
  
  ## total 
  total{
    a <- read.csv("TCGA_479.csv")
    OS <- a[,c(1,912:926)]
    a <- a[,colnames(a) %in% gene_12]
    a <- cbind(OS,a)
    a$riskscore <- NA
    for (i in 1:479){
      a[i,29] <- bt[1]*a[i,17]+bt[2]*a[i,18]+bt[3]*a[i,19]+bt[4]*a[i,20]+bt[5]*a[i,21]+bt[6]*a[i,22]+bt[7]*a[i,23]+bt[8]*a[i,24]+
        bt[9]*a[i,25]+bt[10]*a[i,26]+bt[11]*a[i,27]+bt[12]*a[i,28]
    }
    
    ##生存分析
    library(survival)
    library(survminer)
    a$OS.time1 <- a$OS.time/30
    group <- ifelse(a$riskscore>median(a$riskscore),'High risk','Low risk')
    table(group)
    sfit <- survfit(Surv(OS.time1, OS)~group, data=a)  ## 构建生存对象，进行数据处理
    sfit
    summary(sfit)
    ggsurvplot(sfit, conf.int=F, pval=TRUE,
               legend.title = "Trainset",
               xlab = "Time(months)",
               risk.table = TRUE,
               tables.height = 0.2,
               tables.theme = theme_cleantable(),
               # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
               # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw()
    )
    ##
    library(dplyr)
    library(survival)
    library(survminer)
    library(ggplot2)
    a$OS.time1 <- a$OS.time/30
    a$status <- a$os
    ## Divide the paitents into high and low risk group by the surv.cut fuction.
    ## 按照函数进行分组
    new_function{
      surv.cut <- surv_cutpoint(
        a,
        time = "OS.time1",
        event = "status",
        variables = "riskscore"
      )
      summary(surv.cut)
      plot(surv.cut, "riskscore", palette = "npg")
      
      
      surv.cat <- surv_categorize(surv.cut)
      
      surv.fit <- survfit(Surv(OS.time1, status) ~ riskscore,
                          data = surv.cat)
      ggsurvplot(
        surv.fit,                     # survfit object with calculated statistics.
        risk.table = TRUE,       # show risk table.
        pval = TRUE,             # show p-value of log-rank test.
        # show confidence intervals for 
        # point estimaes of survival curves.
        #xlim = c(0,3000),        # present narrower X axis, but not affect
        # survival estimates.
        break.time.by = 24,    # break X axis in time intervals by 500. 
        risk.table.y.text.col = T, # colour risk table text annotations.
        risk.table.y.text = FALSE # show bars instead of names in text annotations
        # in legend of risk table
      )
      ggsurvplot(surv.fit, conf.int=F, pval=TRUE,
                 legend.title = "Total set",
                 xlab = "Time(months)",
                 legend.labs = c("High risk", "Low risk"),
                 risk.table = TRUE,
                 break.time.by = 24,
                 xlim = c(0,128), ## 2020-02-09学到，自定义X轴
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(),
                 # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
                 # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
                 palette = c("#E7B800", "#2E9FDF"),
                 ggtheme = theme_bw()
      )
    }
    ## 绘制多个时间点ROC曲线，评估模型效能
    library("timeROC")
    a$OS.time2 <- a$OS.time/365
    ROC_a <- timeROC(T = a$OS.time2,delta=a$status,
                     marker=a$riskscore,cause=1,
                     weighting="marginal",
                     times=c(1,3,5),ROC=TRUE)
    ROC_a
    plot(ROC_a,time=1,title=FALSE,lwd=2)
    plot(ROC_a,time=3,col="blue",add=TRUE,title=FALSE,lwd=2)
    plot(ROC_a,time=5,col="black",add=TRUE,title=FALSE,lwd=2)
    legend("bottomright",
           c(paste0("AUC at 1 years: ",round(ROC_a$AUC[1],2)),
             paste0("AUC at 3 years: ",round(ROC_a$AUC[2],2)),
             paste0("AUC at 5 years: ",round(ROC_a$AUC[3],2))),
           col=c("red","blue","black"),lwd=2,bty = "n")
  }
  save(tmp,id1,id2, file = "sample_1_1.Rda")
}

## 提取临床信息文件
tip{
  library(readr)
  a <- read.csv("Pan_cancer_phen.csv")
  colnames(a)
  ### 提取LUAD的
  LUAD <- subset(a,a$cancer.type.abbreviation == "LUAD")
  LUAD <- LUAD[,c(1:5,7,14,26,27)]
  ##
  load("train_17_riskscore.Rda")
  dat <- merge(a,LUAD,by = "X_PATIENT")
  dat1 <- dat[!duplicated(dat$X_PATIENT),]
  
  library(mygene)
  ?query()
  entrez <- queryMany(gene_17, scopes="symbol", fields="entrezgene", species="human")
  gene1 <- entrez$entrezgene
  gene1
  
  ego <- enrichGO(gene1,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL")
  head(ego)
  
  ekegg <- enrichKEGG(gene1,
                      organism = "hsa",
                      keyType = "Kegg")
}
