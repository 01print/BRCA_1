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