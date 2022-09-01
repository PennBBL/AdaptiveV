# Coming up with crowdsourced norms in order to look at ICCs for AdaptiveV CNB vs CAT-CNB analysis
# based on CNB-CAT_time_tables_24May2021.R (but changed the time tables to be accuracy tables)
#
# 06.10.22 --- Akira Di Sandro


library(psych)
library(PerFit)
library(mirt)
library(difR)
library(tidyverse)
library(splitstackshape)
library(dplyr)
library(stats)

# (1) Import Data #################################################################
# Here we read in a bunch of data from files scored by Greg Arbet, plus full 
# files uploaded by Allie & Kosha.

ADT30_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/ADT30_Scored.csv")
ADT30_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/ADT30_Scored.csv")
ADT30_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/ADT30_Scored.csv")
CPF_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 1/CPF_Scored.csv")
CPF_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/CPF_Scored.csv")
CPF_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/CPF_Scored.csv")
CPW_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 1/CPW_Scored.csv")
CPW_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/CPW_Scored.csv")
CPW_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/CPW_Scored.csv")
ER40_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/ER40_Scored.csv")
ER40_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/ER40_Scored.csv")
ER40_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/ER40_Scored.csv")
MEDF_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/MEDF_Scored.csv")
MEDF_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/MEDF_Scored.csv")
MEDF_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/MEDF_Scored.csv")
PLOT24_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 1/PLOT24_Scored.csv")
PLOT24_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/PLOT24_Scored.csv")
PLOT24_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/PLOT24_Scored.csv")
PMAT24_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 1/PMAT24_Scored.csv")
PMAT24_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/PMAT24_Scored.csv")
PMAT24_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/PMAT24_Scored.csv")
PVRT25_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 1/PVRT25_Scored.csv")
PVRT25_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/PVRT25_Scored.csv")
PVRT25_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/PVRT25_Scored.csv")
SVOLT_1 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 1/SVOLT_Scored.csv")
SVOLT_2 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 2/SVOLT_Scored.csv")
SVOLT_3 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 3/SVOLT_Scored.csv")
SVOLT_4 <- read.csv("../CNB-CAT/Data/ScoredFiles_6_21 - Batteries 1-4/Battery 4/SVOLT_Scored.csv")

# pull from battery 1
x <- read.csv("../CNB-CAT/Data/cnb-cat_scored_files/cs_battery_1_05272020.csv")

acc_menu <- x[,grepl("_resp",colnames(x)) | grepl("_corr",colnames(x)) | grepl("_CORR",colnames(x))] 
RT_menu <- x[,grepl("_ttr",colnames(x)) | grepl("_TTR",colnames(x))]

gng_acc <- acc_menu[,1:150]
edisc_acc <- acc_menu[,331:364] - 1 

gng_rt <- RT_menu[,1:150]
edisc_rt <- RT_menu[,384:417] 

GNG_1 <- data.frame(gng_acc,gng_rt)
EDISC_1 <- data.frame(edisc_acc,edisc_rt)

DIGSYM_1 <- x %>% dplyr::select(matches("DIGSYM.DSCOR|DIGSYM.DSMEMCR|DIGSYM.DSMCRRT"))
AIM_1 <- x %>% dplyr::select(matches("AIM.AIMTOT"))

# pull from battery 2
x <- read.csv("../CNB-CAT/Data/cnb-cat_scored_files/cs_battery_2_05272020.csv")

acc_menu <- x[,grepl("_resp",colnames(x)) | grepl("_corr",colnames(x)) | grepl("_CORR",colnames(x))| grepl(".q_",colnames(x))] 
RT_menu <- x[,grepl("_ttr",colnames(x)) | grepl("_TTR",colnames(x))  | grepl("trr",colnames(x))]

lnb_acc <- acc_menu[,41:130]
ddisc_acc <- acc_menu[,131:164] - 1
rdisc_acc <- acc_menu[,165:198] - 1

lnb_rt <- RT_menu[,101:190]
ddisc_rt <- RT_menu[,191:224]
rdisc_rt <- RT_menu[,333:366]

LNB_1 <- data.frame(lnb_acc,lnb_rt)
DDISC_1 <- data.frame(ddisc_acc,ddisc_rt)
RDISC_1 <- data.frame(rdisc_acc,rdisc_rt)

# pull from battery 3
x <- read.csv("../CNB-CAT/Data/cnb-cat_scored_files/cs_battery_3_05272020.csv")

acc_menu <- x[,grepl("_resp",colnames(x)) | grepl("_corr",colnames(x)) | grepl("_CORR",colnames(x)) | grepl(".q_",colnames(x)) & !grepl("ttr",colnames(x))] 
RT_menu <- x[,grepl("_ttr",colnames(x)) | grepl("_TTR",colnames(x))  | grepl("trr",colnames(x))]

ddisc_acc <- acc_menu[,1:34] - 1
edisc_acc <- acc_menu[,135:168] - 1
rdisc_acc <- acc_menu[,169:202] - 1

ddisc_rt <- RT_menu[,169:202]
edisc_rt <- RT_menu[,303:336]
rdisc_rt <- RT_menu[,337:370]

DDISC_2 <- data.frame(ddisc_acc,ddisc_rt)
EDISC_2 <- data.frame(edisc_acc,edisc_rt)
RDISC_2 <- data.frame(rdisc_acc,rdisc_rt)

DIGSYM_3 <- x %>% dplyr::select(matches("DIGSYM.DSCOR|DIGSYM.DSMEMCR|DIGSYM.DSMCRRT"))
AIM_3 <- x %>% dplyr::select(matches("AIM.B.AIMTOT"))

# pull from battery 4,5,6
x <- read.csv("../CNB-CAT/Data/cnb-cat_scored_files/cs_battery_4_05272020.csv")
x5 <- read.csv("../CNB-CAT/Data/cnb-cat_scored_files/cs_battery_5_11182020.csv")
x6 <- read.csv("../CNB-CAT/Data/cnb-cat_scored_files/cs_battery_6_11182020.csv")
acc5 <- x5[,grepl("CORR",colnames(x5))]
acc6 <- x6[,grepl("CORR",colnames(x6))]
spe5 <- x5[,grepl("TTR",colnames(x5))]
spe6 <- x6[,grepl("TTR",colnames(x6))]
cpf4 <- data.frame(acc5[,1:40],spe5[,450:489])
cpf5 <- data.frame(acc6[,1:40],spe6[,1:40])
cpw2 <- data.frame(acc5[,41:88],spe5[,201:248])
cpw3 <- data.frame(acc6[,41:88],spe6[,241:288])
volt2 <- data.frame(acc5[,89:108],spe5[,490:509])
volt4 <- data.frame(acc6[,89:108],spe6[,490:509])

acc_menu <- x[,grepl("_resp",colnames(x)) | grepl("_corr",colnames(x)) | grepl("_CORR",colnames(x)) | grepl(".q_",colnames(x)) & !grepl("ttr",colnames(x))] 
RT_menu <- x[,grepl("_ttr",colnames(x)) | grepl("_TTR",colnames(x))  | grepl("trr",colnames(x))]

ctom_resp <- acc_menu[,1:42]
ctom_key <- c(3,1,3,1,2,1,3,2,2,1,3,2,1,3,2,2,3,1,3,1,1,3,2,1,1,1,3,3,2,3,2,1,2,3,2,3,1,3,1,1,3,2)

ctom_acc <- matrix(NA,nrow(ctom_resp),ncol(ctom_resp))
for (i in 1:nrow(ctom_resp)) {ctom_acc[i,] <- as.numeric(ctom_resp[i,] == ctom_key)}
ddisc_acc <- acc_menu[,43:76] - 1
edisc_acc <- acc_menu[,177:210] - 1
rdisc_acc <- acc_menu[,211:244] - 1
cpt_acc <- acc_menu[,245:424]

ctom_rt <- acc_menu[,1:42]   # WARNING WARNING WARNING WARNING WARNING
ddisc_rt <- RT_menu[,166:199]
edisc_rt <- RT_menu[,300:333]
rdisc_rt <- RT_menu[,334:367]
cpt_rt <- RT_menu[,436:615]

CTOM_1 <- data.frame(ctom_acc,ctom_rt)
DDISC_3 <- data.frame(ddisc_acc,ddisc_rt)
EDISC_3 <- data.frame(edisc_acc,edisc_rt)
RDISC_3 <- data.frame(rdisc_acc,rdisc_rt)
CPT_1 <- data.frame(cpt_acc,cpt_rt)

CPF_4 <- cpf4                    # lower-case are from other script 3November
CPF_4[CPF_4 == TRUE] <- 1
CPF_4[CPF_4 == FALSE] <- 0
CPF_5 <- cpf5
CPF_5[CPF_5 == TRUE] <- 1
CPF_5[CPF_5 == FALSE] <- 0
 
CPW_2 <- CPW_2[,c(grep("_CORR",colnames(CPW_2)),grep("_TTR",colnames(CPW_2)))]
colnames(CPW_2) <- colnames(cpw2)
CPW_2 <- rbind(CPW_2,cpw2)

CPW_3 <- CPW_3[,c(grep("_CORR",colnames(CPW_3)),grep("_TTR",colnames(CPW_3)))]
colnames(CPW_3) <- colnames(cpw3)
CPW_3 <- rbind(CPW_3,cpw3) 
 
SVOLT_2 <- SVOLT_2[,c(grep("_CORR",colnames(SVOLT_2)),grep("_TTR",colnames(SVOLT_2)))]
colnames(SVOLT_2) <- colnames(volt2)
SVOLT_2 <- rbind(SVOLT_2,volt2) 

SVOLT_4 <- SVOLT_4[,c(grep("_CORR",colnames(SVOLT_4)),grep("_TTR",colnames(SVOLT_4)))]
colnames(SVOLT_4) <- colnames(volt4)
SVOLT_4 <- rbind(SVOLT_4,volt4)


# 7/22/22, using ptp-fpf for cpt_acc to match what we have from adaptiveV data.
cpt_all <- x %>% dplyr::select(matches("CPT")) %>% mutate(cpt_acc = (SPCPTNL.SCPT_TP/60)-(SPCPTNL.SCPT_FP/120))
cpt_all$cpt_acc <- ifelse(cpt_all$cpt_acc < 0, 0, cpt_all$cpt_acc)

# many lines below  -  marked {{Earlier}}  -  we go through the performance 
# validity estimator, but the lines immediately below make the correction after 
# the second round of crowdsourcing, affecting only some of the tests in 
# the {{Earlier}} script.  This new one includes the new versions of the 
# CPF ("4" and "5"), as well as adding cases to the CPW and SVOLT.  Now start 
# the new lines (not the original):
# * (a) CPF4,5 ----

dat <-CPF_4[,c(grep("_CORR",colnames(CPF_4)),grep("_TTR",colnames(CPF_4)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPF_4<- data.frame(CPF_4,sc)
# get rid of bottom 5%
x <-CPF_4
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPF_4 <- x[,colnames(x) != "PFscores"] # only 62.253% of respondents left

num_na <- function(out, acc, pfit1, pfit2) {
  nout <- length(out[is.na(out)])
  nacc <- length(acc[is.na(acc)])
  npf1 <- length(pfit1[is.na(pfit1)])
  npf2 <- length(pfit2[is.na(pfit2)])
  na_list <- c(nout, nacc, npf1, npf2)
  print(na_list)
}

dat <-CPF_5[,c(grep("_CORR",colnames(CPF_5)),grep("_TTR",colnames(CPF_5)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPF_5<- data.frame(CPF_5,sc)
x <-CPF_5
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPF_5 <- x[,colnames(x) != "PFscores"] # only 78.571% of respondents left

# * (b) CPW2,3 ----
dat <-CPW_2[,c(grep("_CORR",colnames(CPW_2)),grep("_TTR",colnames(CPW_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPW_2<- data.frame(CPW_2,sc)
x <-CPW_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPW_2 <- x[,colnames(x) != "PFscores"] # only 76.080%, the NAs seem to be coming from both pfit scores

dat <-CPW_3[,c(grep("_CORR",colnames(CPW_3)),grep("_TTR",colnames(CPW_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items] # 297 NaNs
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 297 NAs
pfit1 <- r.pbis(dat2)$PFscores # 373 NAs
pfit2 <- E.KB(dat2)$PFscores # 373 NAs
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPW_3<- data.frame(CPW_3,sc)
x <-CPW_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPW_3 <- x[,colnames(x) != "PFscores"]

# * (c) SVOLT2,4 ----
dat <-SVOLT_2[,c(grep("_CORR",colnames(SVOLT_2)),grep("_TTR",colnames(SVOLT_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) #523 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) #528 NAs
pfit1 <- r.pbis(dat2)$PFscores # 586 NAs
pfit2 <- E.KB(dat2)$PFscores # 586 NAs
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
SVOLT_2<- data.frame(SVOLT_2,sc)
x <-SVOLT_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
SVOLT_2 <- x[,colnames(x) != "PFscores"]

dat <-SVOLT_4[,c(grep("_CORR",colnames(SVOLT_4)),grep("_TTR",colnames(SVOLT_4)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) #351 NaN
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 352 NA
pfit1 <- r.pbis(dat2)$PFscores # 394 NA
pfit2 <- E.KB(dat2)$PFscores #394 NA
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
SVOLT_4 <- data.frame(SVOLT_4,sc)
x <-SVOLT_4
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
SVOLT_4 <- x[,colnames(x) != "PFscores"]


# (2) Performance Validity #####################################################
# {{Earlier}} Now we run our performance validity estimator on each test and add those 
# performance validity scores to each matrix 

# * (a) ADT30 ----
dat <- ADT30_1[,c(grep("_CORR",colnames(ADT30_1)),grep("_TTR",colnames(ADT30_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
ADT30_1 <- data.frame(ADT30_1,sc)
# get rid of the bottom 5%
x <- ADT30_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
ADT30_1 <- x[, colnames(x) != "PFscores"]

dat <-ADT30_2[,c(grep("_CORR",colnames(ADT30_2)),grep("_TTR",colnames(ADT30_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
ADT30_2<- data.frame(ADT30_2,sc)
# get rid of bottom 5%
x <-ADT30_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
ADT30_2 <- x[,colnames(x) != "PFscores"]

dat <-ADT30_3[,c(grep("_CORR",colnames(ADT30_3)),grep("_TTR",colnames(ADT30_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
ADT30_3<- data.frame(ADT30_3,sc)
# get rid of bottom 5%
x <-ADT30_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
ADT30_3 <- x[,colnames(x) != "PFscores"]

# * (b) CPF ----
dat <-CPF_1[,c(grep("_CORR",colnames(CPF_1)),grep("_TTR",colnames(CPF_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPF_1<- data.frame(CPF_1,sc)
x <-CPF_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPF_1 <- x[,colnames(x) != "PFscores"]

dat <-CPF_2[,c(grep("_CORR",colnames(CPF_2)),grep("_TTR",colnames(CPF_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPF_2<- data.frame(CPF_2,sc)
x <-CPF_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPF_2 <- x[,colnames(x) != "PFscores"]

dat <-CPF_3[,c(grep("_CORR",colnames(CPF_3)),grep("_TTR",colnames(CPF_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 14 NA
pfit2 <- E.KB(dat2)$PFscores # 14 NA
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPF_3<- data.frame(CPF_3,sc)
x <-CPF_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPF_3 <- x[,colnames(x) != "PFscores"]

# * (c) CPW ----
dat <-CPW_1[,c(grep("_CORR",colnames(CPW_1)),grep("_TTR",colnames(CPW_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 47 NA
pfit2 <- E.KB(dat2)$PFscores # 47 NA
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
CPW_1<- data.frame(CPW_1,sc)
x <-CPW_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
CPW_1 <- x[,colnames(x) != "PFscores"]

# dat <-CPW_2[,c(grep("_CORR",colnames(CPW_2)),grep("_TTR",colnames(CPW_2)))]
# items <- ncol(dat)/2
# dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
# res <- matrix(NA,dim(dat)[1],items)
# for (j in 1:items) {
# mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
# res[,j] <- scale(residuals(mod,na.action=na.exclude))}
# res2 <- res
# res2[abs(res2) < 2] <- 0
# res2[abs(res2) > 2] <- 1
# outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
# dat2 <- dat[,1:items]
# acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
# pfit1 <- r.pbis(dat2)$PFscores
# pfit2 <- E.KB(dat2)$PFscores
# sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
# CPW_2<- data.frame(CPW_2,sc)
# x <-CPW_2
# qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
# x <- x[which(x$PFscores > qu),]
# CPW_2 <- x[,colnames(x) != "PFscores"]

# dat <-CPW_3[,c(grep("_CORR",colnames(CPW_3)),grep("_TTR",colnames(CPW_3)))]
# items <- ncol(dat)/2
# dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
# res <- matrix(NA,dim(dat)[1],items)
# for (j in 1:items) {
# mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
# res[,j] <- scale(residuals(mod,na.action=na.exclude))}
# res2 <- res
# res2[abs(res2) < 2] <- 0
# res2[abs(res2) > 2] <- 1
# outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
# dat2 <- dat[,1:items]
# acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
# pfit1 <- r.pbis(dat2)$PFscores
# pfit2 <- E.KB(dat2)$PFscores
# sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
# CPW_3<- data.frame(CPW_3,sc)
# x <-CPW_3
# qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
# x <- x[which(x$PFscores > qu),]
# CPW_3 <- x[,colnames(x) != "PFscores"]

# * (d) ER40 ----
dat <-ER40_1[,c(grep("_CORR",colnames(ER40_1)),grep("_TTR",colnames(ER40_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
ER40_1<- data.frame(ER40_1,sc)
x <-ER40_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
ER40_1 <- x[,colnames(x) != "PFscores"]

dat <-ER40_2[,c(grep("_CORR",colnames(ER40_2)),grep("_TTR",colnames(ER40_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
ER40_2<- data.frame(ER40_2,sc)
x <- ER40_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
ER40_2 <- x[,colnames(x) != "PFscores"]

dat <-ER40_3[,c(grep("_CORR",colnames(ER40_3)),grep("_TTR",colnames(ER40_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
ER40_3<- data.frame(ER40_3,sc)
x <- ER40_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
ER40_3 <- x[,colnames(x) != "PFscores"]

# * (e) MEDF ----
dat <-MEDF_1[,c(grep("_CORR",colnames(MEDF_1)),grep("_TTR",colnames(MEDF_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
MEDF_1<- data.frame(MEDF_1,sc)
x <- MEDF_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
MEDF_1 <- x[,colnames(x) != "PFscores"]

dat <-MEDF_2[,c(grep("_CORR",colnames(MEDF_2)),grep("_TTR",colnames(MEDF_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
MEDF_2<- data.frame(MEDF_2,sc)
x <- MEDF_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
MEDF_2 <- x[,colnames(x) != "PFscores"]

dat <-MEDF_3[,c(grep("_CORR",colnames(MEDF_3)),grep("_TTR",colnames(MEDF_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores
pfit2 <- E.KB(dat2)$PFscores
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
MEDF_3<- data.frame(MEDF_3,sc)
x <- MEDF_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
MEDF_3 <- x[,colnames(x) != "PFscores"]

# * (f) PLOT24 ----
dat <-PLOT24_1[,c(grep("_CORR",colnames(PLOT24_1)),grep("_TTR",colnames(PLOT24_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 22 NAs
pfit2 <- E.KB(dat2)$PFscores # 22 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PLOT24_1<- data.frame(PLOT24_1,sc)
x <- PLOT24_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PLOT24_1 <- x[,colnames(x) != "PFscores"]

dat <-PLOT24_2[,c(grep("_CORR",colnames(PLOT24_2)),grep("_TTR",colnames(PLOT24_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 10 NAs
pfit2 <- E.KB(dat2)$PFscores # 10 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PLOT24_2<- data.frame(PLOT24_2,sc)
x <- PLOT24_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PLOT24_2 <- x[,colnames(x) != "PFscores"]

dat <-PLOT24_3[,c(grep("_CORR",colnames(PLOT24_3)),grep("_TTR",colnames(PLOT24_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 16 NAs
pfit2 <- E.KB(dat2)$PFscores # 16 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PLOT24_3<- data.frame(PLOT24_3,sc)
x <- PLOT24_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PLOT24_3 <- x[,colnames(x) != "PFscores"]

# * (g) PMAT24 ----
dat <-PMAT24_1[,c(grep("_CORR",colnames(PMAT24_1)),grep("_TTR",colnames(PMAT24_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 10 NAs
pfit2 <- E.KB(dat2)$PFscores # 10 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PMAT24_1<- data.frame(PMAT24_1,sc)
x <- PMAT24_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PMAT24_1 <- x[,colnames(x) != "PFscores"]

dat <-PMAT24_2[,c(grep("_CORR",colnames(PMAT24_2)),grep("_TTR",colnames(PMAT24_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 7 NAs
pfit2 <- E.KB(dat2)$PFscores # 7 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PMAT24_2<- data.frame(PMAT24_2,sc)
x <- PMAT24_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PMAT24_2 <- x[,colnames(x) != "PFscores"]

dat <-PMAT24_3[,c(grep("_CORR",colnames(PMAT24_3)),grep("_TTR",colnames(PMAT24_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 8 NAs
pfit2 <- E.KB(dat2)$PFscores # 8 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PMAT24_3<- data.frame(PMAT24_3,sc)
x <- PMAT24_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PMAT24_3 <- x[,colnames(x) != "PFscores"]

# * (h) PVRT25 ----
dat <-PVRT25_1[,c(grep("_CORR",colnames(PVRT25_1)),grep("_TTR",colnames(PVRT25_1)))]          # grab all the columns that include "_CORR" (correct or not) and "_TTR" (time to respond)
items <- ncol(dat)/2                                                                          # items is the number of columns in data divided by two (probably because there are as many ttr as corr & we mostly just care about the corr)
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])                                   # replacing the second half of all the rows with ln() of the original values in those rows
res <- matrix(NA,dim(dat)[1],items)                                                           # create a matrix of all NA entries with dimension nrow(dat), items (half of ncol(dat))
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)                              # linear model of ttr depending on corr
res[,j] <- scale(residuals(mod,na.action=na.exclude))}                                        # in jth row, store scaled residuals into res
res2 <- res                                                                                   # copy of res into res2
res2[abs(res2) < 2] <- 0                                                                      # everything within -2 and 2 gets changed into a 0
res2[abs(res2) > 2] <- 1                                                                      # everything <-2 and >2 gets changed into a 1, -2 and 2 get left alone
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)                                           # list made up of one minus mean of each row of res2
dat2 <- dat[,1:items]                                                                         # only store the corr (first half) aspects of dat to dat2
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # take the row means of the rows of dat2 where the col means are greater than min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))] --- the minimum of the three largest means from the columns of dat2
pfit1 <- r.pbis(dat2)$PFscores                                                                # donlon and fischer's personal biserial statistic, assessing personal fit
pfit2 <- E.KB(dat2)$PFscores                                                                  # kane and brennan's person-fit statistics, specifically the dependability statistic
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)                     # sc is a column vector 
PVRT25_1<- data.frame(PVRT25_1,sc)                                                            # replace PVRT25_1 with a df of itself and sc
x <- PVRT25_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PVRT25_1 <- x[,colnames(x) != "PFscores"]

dat <-PVRT25_2[,c(grep("_CORR",colnames(PVRT25_2)),grep("_TTR",colnames(PVRT25_2)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 14 NAs
pfit2 <- E.KB(dat2)$PFscores # 14 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PVRT25_2<- data.frame(PVRT25_2,sc)
x <- PVRT25_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PVRT25_2 <- x[,colnames(x) != "PFscores"]

dat <-PVRT25_3[,c(grep("_CORR",colnames(PVRT25_3)),grep("_TTR",colnames(PVRT25_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 10 NAs
pfit2 <- E.KB(dat2)$PFscores # 10 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
PVRT25_3<- data.frame(PVRT25_3,sc)
x <- PVRT25_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
PVRT25_3 <- x[,colnames(x) != "PFscores"]

# * (i) SVOLT ----
dat <-SVOLT_1[,c(grep("_CORR",colnames(SVOLT_1)),grep("_TTR",colnames(SVOLT_1)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 8 NAs
pfit2 <- E.KB(dat2)$PFscores # 8 NAs
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
SVOLT_1<- data.frame(SVOLT_1,sc)
x <- SVOLT_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
SVOLT_1 <- x[,colnames(x) != "PFscores"]

# dat <-SVOLT_2[,c(grep("_CORR",colnames(SVOLT_2)),grep("_TTR",colnames(SVOLT_2)))]
# items <- ncol(dat)/2
# dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
# res <- matrix(NA,dim(dat)[1],items)
# for (j in 1:items) {
# mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
# res[,j] <- scale(residuals(mod,na.action=na.exclude))}
# res2 <- res
# res2[abs(res2) < 2] <- 0
# res2[abs(res2) > 2] <- 1
# outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
# dat2 <- dat[,1:items]
# acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
# pfit1 <- r.pbis(dat2)$PFscores
# pfit2 <- E.KB(dat2)$PFscores
# sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
# SVOLT_2<- data.frame(SVOLT_2,sc)
# x <- SVOLT_2
# qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
# x <- x[which(x$PFscores > qu),]
# SVOLT_2 <- x[,colnames(x) != "PFscores"]

dat <-SVOLT_3[,c(grep("_CORR",colnames(SVOLT_3)),grep("_TTR",colnames(SVOLT_3)))]
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
pfit1 <- r.pbis(dat2)$PFscores # 34 NAs
pfit2 <- E.KB(dat2)$PFscores # 34 NAs
sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
SVOLT_3<- data.frame(SVOLT_3,sc)
x <- SVOLT_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
SVOLT_3 <- x[,colnames(x) != "PFscores"]

# dat <-SVOLT_4[,c(grep("_CORR",colnames(SVOLT_4)),grep("_TTR",colnames(SVOLT_4)))]
# items <- ncol(dat)/2
# dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
# res <- matrix(NA,dim(dat)[1],items)
# for (j in 1:items) {
# mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
# res[,j] <- scale(residuals(mod,na.action=na.exclude))}
# res2 <- res
# res2[abs(res2) < 2] <- 0
# res2[abs(res2) > 2] <- 1
# outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE)
# dat2 <- dat[,1:items]
# acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
# pfit1 <- r.pbis(dat2)$PFscores
# pfit2 <- E.KB(dat2)$PFscores
# sc <- (0.42*outlier_score_2cut) + (0.02*acc3e) + (0.05*pfit1) + (0.50*pfit2)
# SVOLT_4<- data.frame(SVOLT_4,sc)
# x <- SVOLT_4
# qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
# x <- x[which(x$PFscores > qu),]
# SVOLT_4 <- x[,colnames(x) != "PFscores"]

# * (j) D,E,RDISC ----
dat <-DDISC_1
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 388 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 388 NAs
pfit1 <- r.pbis(dat2)$PFscores # 465 NAs
pfit2 <- E.KB(dat2)$PFscores # 465 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
DDISC_1<- data.frame(DDISC_1,sc)
x <- DDISC_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
DDISC_1 <- x[,colnames(x) != "PFscores"]

dat <-EDISC_1
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 414 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 414 NAs
pfit1 <- r.pbis(dat2)$PFscores # 953 NAs
pfit2 <- E.KB(dat2)$PFscores # 953 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
EDISC_1<- data.frame(EDISC_1,sc)
x <- EDISC_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
EDISC_1 <- x[,colnames(x) != "PFscores"]

dat <-RDISC_1
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 430 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 430 NAs
pfit1 <- r.pbis(dat2)$PFscores # 507 NAs
pfit2 <- E.KB(dat2)$PFscores # 507 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
RDISC_1<- data.frame(RDISC_1,sc)
x <- RDISC_1
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
RDISC_1 <- x[,colnames(x) != "PFscores"]

dat <-DDISC_2
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 276 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])  # 276 NaNs
pfit1 <- r.pbis(dat2)$PFscores # 361 NAs
pfit2 <- E.KB(dat2)$PFscores # 361 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
DDISC_2<- data.frame(DDISC_2,sc)
x <- DDISC_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
DDISC_2 <- x[,colnames(x) != "PFscores"] # 65.831% rather than 95%

dat <-EDISC_2
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 294 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 294 NAs
pfit1 <- r.pbis(dat2)$PFscores # 590 NAs
pfit2 <- E.KB(dat2)$PFscores # 590 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
EDISC_2<- data.frame(EDISC_2,sc)
x <- EDISC_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
EDISC_2 <- x[,colnames(x) != "PFscores"]

dat <-RDISC_2
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 296 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 296 NAs
pfit1 <- r.pbis(dat2)$PFscores # 415 NAs
pfit2 <- E.KB(dat2)$PFscores # 415 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
RDISC_2<- data.frame(RDISC_2,sc)
x <- RDISC_2
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
RDISC_2 <- x[,colnames(x) != "PFscores"]

dat <-DDISC_3
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 190 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 190 NAs
pfit1 <- r.pbis(dat2)$PFscores # 291 NAs
pfit2 <- E.KB(dat2)$PFscores # 291 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
DDISC_3<- data.frame(DDISC_3,sc)
x <- DDISC_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
DDISC_3 <- x[,colnames(x) != "PFscores"]

dat <-EDISC_3
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 204 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 204 NAs
pfit1 <- r.pbis(dat2)$PFscores # 590 NAs
pfit2 <- E.KB(dat2)$PFscores # 590 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
EDISC_3<- data.frame(EDISC_3,sc)
x <- EDISC_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
EDISC_3 <- x[,colnames(x) != "PFscores"]

dat <-RDISC_3
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 209 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 209 NAs
pfit1 <- r.pbis(dat2)$PFscores # 353 NAs
pfit2 <- E.KB(dat2)$PFscores # 353 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
RDISC_3<- data.frame(RDISC_3,sc)
x <- RDISC_3
qu <- quantile(x$PFscores,0.05,na.rm=TRUE)
x <- x[which(x$PFscores > qu),]
RDISC_3 <- x[,colnames(x) != "PFscores"]

# * (k) GNG,LNB,CTOM,CPT ----
dat <-GNG_1                      # use ALL data
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 337 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 308 NAs
pfit1 <- r.pbis(dat2)$PFscores # 360 NAs
pfit2 <- E.KB(dat2)$PFscores # 360 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
GNG_1 <- data.frame(GNG_1,sc)
GNG_1 <-GNG_1[,grepl("_CORR",colnames(GNG_1))]
#write.csv(GNG_1,"Data/full_files/GNG_ttr_CNB-CAT.csv",na="",row.names=FALSE)

# dat <-LNB_1                      # use ALL data
# items <- ncol(dat)/2
# dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
# res <- matrix(NA,dim(dat)[1],items)
# for (j in 1:items) {
# mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
# res[,j] <- scale(residuals(mod,na.action=na.exclude))}
# res2 <- res
# res2[abs(res2) < 2] <- 0
# res2[abs(res2) > 2] <- 1
# outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 393 NaNs
# dat2 <- dat[,1:items]
# acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) #383 NAs
# pfit1 <- r.pbis(dat2)$PFscores # 540 NAs
# pfit2 <- E.KB(dat2)$PFscores # 540 NAs
# sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
# LNB_1 <- data.frame(LNB_1,sc)
# LNB_1 <- LNB_1[,grepl("_CORR",colnames(LNB_1))]
#write.csv(LNB_1,"Data/full_files/LNB_ttr_CNB-CAT.csv",na="",row.names=FALSE)


dat <-CPT_1                      # use ALL data
items <- ncol(dat)/2
dat[,(items+1):(2*items)] <- log(dat[,(items+1):(2*items)])
res <- matrix(NA,dim(dat)[1],items)
for (j in 1:items) {
mod <- lm(dat[,(j+items)]~dat[,j],data=dat,na.action=na.exclude)
res[,j] <- scale(residuals(mod,na.action=na.exclude))}
res2 <- res
res2[abs(res2) < 2] <- 0
res2[abs(res2) > 2] <- 1
outlier_score_2cut <- 1 - rowMeans(res2,na.rm=TRUE) # 304 NaNs
dat2 <- dat[,1:items]
acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))]) # 297 NAs
pfit1 <- r.pbis(dat2)$PFscores # 307 NAs
pfit2 <- E.KB(dat2)$PFscores # 307 NAs
sc <- (0.34*outlier_score_2cut) + (0*acc3e) + (0.22*pfit1) + (0.44*pfit2)
CPT_1 <- data.frame(CPT_1,sc)
CPT_1 <- CPT_1[,grepl("_CORR",colnames(CPT_1))]
#write.csv(CPT_1,"Data/full_files/CPT_ttr_CNB-CAT.csv",na="",row.names=FALSE)

# dat <-CTOM_1                  # no SPEED scores
# items <- ncol(dat)/2
# dat2 <- dat[,1:items]
# acc3e <- rowMeans(dat2[,colMeans(dat2,na.rm=TRUE) >= min(tail(sort(colMeans(dat2,na.rm=TRUE)),3))])
# pfit1 <- r.pbis(dat2)$PFscores
# pfit2 <- E.KB(dat2)$PFscores
# sc <- (0.333333*acc3e) + (0.333333*pfit1) + (0.333333*pfit2)
# CTOM_1<- data.frame(CTOM_1,sc)


# 6.16.22 I think for ICCs I just need mean and sd of tot scores in general
{ # acc/score
  ADT30_1$ADT30_1_CORR <- rowSums(ADT30_1 %>% dplyr::select(matches("CORR")))
  ADT30_2$ADT30_2_CORR <- rowSums(ADT30_2 %>% dplyr::select(matches("CORR")))
  ADT30_3$ADT30_3_CORR <- rowSums(ADT30_3 %>% dplyr::select(matches("CORR")))
  
  ADT_all <- c(ADT30_1$ADT30_1_CORR,ADT30_2$ADT30_2_CORR,ADT30_3$ADT30_3_CORR)
  score_sum <- data.frame(mean = mean(ADT_all,na.rm = T),sd = sd(ADT_all,na.rm = T))
  
  # rt/speed
  # for RT, i don't wanna combine by taking medians, i want to link all items (use linking scripts + google doc for these)
  # first 12 of each form are the same items
  ADT_1tojoin <- data.frame(ADT30_1 %>% dplyr::select(matches("TTR")),matrix(NA, nrow = nrow(ADT30_1), ncol = 30))
  ADT_2tojoin <- data.frame(ADT30_2 %>% dplyr::select(matches("TTR")) %>% dplyr::select(ADT30_QID001_TTR:ADT30_QID015_TTR),
                            matrix(NA, nrow = nrow(ADT30_2), ncol = 15),
                            ADT30_2 %>% dplyr::select(matches("TTR")) %>% dplyr::select(ADT30_QID016_TTR:ADT30_QID030_TTR),
                            matrix(NA, nrow = nrow(ADT30_2), ncol = 15))
  ADT_3tojoin <- data.frame(ADT30_3 %>% dplyr::select(matches("TTR")) %>% dplyr::select(ADT30_QID001_TTR:ADT30_QID015_TTR),
                            matrix(NA, nrow = nrow(ADT30_3), ncol = 30),
                            ADT30_3 %>% dplyr::select(matches("TTR")) %>% dplyr::select(ADT30_QID016_TTR:ADT30_QID030_TTR))
  
  names <- c("ADT60_20","ADT60_56","ADT60_60","ADT60_25","ADT60_13","ADT60_4","ADT60_28","ADT60_59","ADT60_12","ADT60_22","ADT60_16","ADT60_11","ADT60_33","ADT60_17","ADT60_54","ADT60_39","ADT60_55","ADT60_37","ADT60_50","ADT60_18","ADT60_14","ADT60_26","ADT60_58","ADT60_48","ADT60_41","ADT60_57","ADT60_27","ADT60_44","ADT60_47","ADT60_2","ADT60_34","ADT60_51","ADT60_30","ADT60_36","ADT60_43","ADT60_19","ADT60_29","ADT60_21","ADT60_3","ADT60_1","ADT60_5","ADT60_35","ADT60_8","ADT60_42","ADT60_40","ADT60_15","ADT60_7","ADT60_32","ADT60_38","ADT60_53","ADT60_10","ADT60_6","ADT60_9","ADT60_52","ADT60_24","ADT60_31","ADT60_45","ADT60_46","ADT60_49","ADT60_23")
  names(ADT_1tojoin) <- names   # ADT30_1 items are named ADT30_QID001 - ADT30_QID030
  names(ADT_2tojoin) <- names   # ADT30_1 items are named ADT30_QID001 - ADT30_QID012, ADT30_QID031 - ADT30_QID048
  names(ADT_3tojoin) <- names   # ADT30_1 items are named ADT30_QID001 - ADT30_QID012, ADT30_QID049 - ADT30_QID066
  ADT_all <- rbind(ADT_1tojoin,ADT_2tojoin,ADT_3tojoin)
  
  # add a QA, limiting RT to something in the range 100 ms - 120000 ms
  ADT_all[ADT_all > 120000 | ADT_all < 100] <- NA
  
  ADT_rt <- data.frame(matrix(NA,nrow = ncol(ADT_all),ncol = 2))
  row.names(ADT_rt) <- names(ADT_all)
  names(ADT_rt) <- c("mean","sd")
  
  for (i in 1:ncol(ADT_all)) {
    # SD and mean for each item
    ADT_rt[i,] <- c(mean(ADT_all[,i],na.rm=T),sd(ADT_all[,i],na.rm=T))
  }
}

{ # acc/score
  AIM_all <- c(AIM_1$AIM.AIMTOT,AIM_3$AIM.B.AIMTOT)
  score_sum[2,] <- data.frame(mean = mean(AIM_all,na.rm = T),sd = sd(AIM_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  CPF_1$CPF_1_CORR <- rowSums(CPF_1)
  CPF_4$CPF_4_CORR <- rowSums(CPF_4)
  CPF_5$CPF_5_CORR <- rowSums(CPF_5)
  
  CPF_all <- c(CPF_1$CPF_1_CORR,CPF_4$CPF_4_CORR,CPF_5$CPF_5_CORR)
  score_sum[3,] <- data.frame(mean = mean(CPF_all,na.rm = T),sd = sd(CPF_all,na.rm = T))
  
  # rt/speed
  CPF_1tojoin <- data.frame(CPF_1 %>% dplyr::select(matches("TTR")),matrix(NA, nrow = nrow(CPF_1), ncol = 55))
  CPF_4tojoin <- data.frame(matrix(NA,nrow = nrow(CPF_4),ncol = 95))
  CPF_4tojoin[,c(1,2,22,36,5,14,7,19,13,21,11,12,41:68)] <- CPF_4 %>% dplyr::select(matches("TTR"))
  CPF_5tojoin <- data.frame(matrix(NA,nrow = nrow(CPF_5),ncol = 95))
  CPF_5tojoin[,c(1:2,22,36,5,14,7,19,13,21,11:12,69:72,45,73:95)] <- CPF_5 %>% dplyr::select(matches("TTR"))
  
  names <- c(paste0("CPF_v1_Item_",1:40),"cpf-b_10",paste0("cpf-b_",14:40),"cpf-b_3","cpf-b_4","cpf-b_1","cpf-b_2","cpf-b_11","cpf-b_5","cpf-b_6","cpf-b_8","cpf-b_9","cpf-b_12","cpf-b_13","cpfd-a_11","cpfd-a_1","cpfd-b_37","cpfd-a_36","cpfd-b_12","cpfd-a_6","cpfd-b_11","cpfd-a_8","cpfd-a_31","cpfd-b_21","cpfd-a_38","cpfd-a_5","cpfd-a_40","cpfd-a_35","cpfd-b-28","cpfd-a_23")
  names(CPF_1tojoin) <- names   
  names(CPF_4tojoin) <- names   
  names(CPF_5tojoin) <- names   
  CPF_all <- rbind(CPF_1tojoin,CPF_4tojoin,CPF_5tojoin)
  
  # add a QA, limiting RT to something in the range 100 ms - 120000 ms
  CPF_all[CPF_all > 120000 | CPF_all < 100] <- NA
  
  CPF_rt <- data.frame(matrix(NA,nrow = ncol(CPF_all),ncol = 2))
  row.names(CPF_rt) <- names(CPF_all)
  names(CPF_rt) <- c("mean","sd")
  
  for (i in 1:ncol(CPF_all)) {
    # SD and mean for each item
    CPF_rt[i,] <- c(mean(CPF_all[,i],na.rm=T),sd(CPF_all[,i],na.rm=T))
  }
}

{ # acc/score
  CPT_all <- cpt_all$cpt_acc
  score_sum[4,] <- data.frame(mean = mean(CPT_all,na.rm = T),sd = sd(CPT_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  CPW_1$CPW_1_CORR <- rowSums(CPW_1)
  CPW_2$CPW_2_CORR <- rowSums(CPW_2)
  CPW_3$CPW_3_CORR <- rowSums(CPW_3)
  
  CPW_all <- c(CPW_1$CPW_1_CORR,CPW_2$CPW_2_CORR,CPW_3$CPW_3_CORR)
  score_sum[5,] <- data.frame(mean = mean(CPW_all,na.rm = T),sd = sd(CPW_all,na.rm = T))
  
  # rt/speed
  CPW_1tojoin <- data.frame(CPW_1 %>% dplyr::select(matches("TTR")),matrix(NA, nrow = nrow(CPW_1), ncol = 72))
  CPW_2tojoin <- data.frame(matrix(NA,nrow = nrow(CPW_2),ncol = 120))
  CPW_2tojoin[,c(49:52,6,14,53,24,54:57,3,58:61,38,62:66,36,67:70,27,71:72,44,73:76,1,45,77:79,33,80:84,21)] <- CPW_2 %>% dplyr::select(matches("TTR"))
  CPW_3tojoin <- data.frame(matrix(NA,nrow = nrow(CPW_3),ncol = 120))
  CPW_3tojoin[,c(85:88,38,45,89,14,90:95,27,96:99,1,100:104,21,105:106,24,107:113,36,114:115,33,116,6,117,3,118:120,44)] <- CPW_3 %>% dplyr::select(matches("TTR"))
  
  names <- c("Egg","Law","Nurse","Madness","Truth","Hobby","Tower","Pride","Flag","Hint","Frog","Mayor","Cheese","Misconception","Vanity","Ring","Socialist","Soul","Hour","Umbrella","Virtue","Figment","Memory","Economy","Distraction","Laundry","Brick","Unbeliever","Namesake","Smile","Style","Scale","Nation","Exclusion","Cost","Event","Misery","Reaction","Gravy","Jeopardy","Night","Aptitude","Hope","Lemon","Reminder","Decree","Fantasy","Essence",
             "Impulse","Cotton","Ego","Bereavement","Lap","Thought","Proxy","Hate","Quality","Clemency","Direction","Soap","Fate","Opinion","Volume","Idiom","Nonsense","Surtax","Bloom","Mercy","Increment","Effort","Box","Equity","Saw","Illusion","Blood","Wistfulness","Welfare","Station","Miracle","Kitten","Power","Moral","Rating","Cherry",
             "Sleeve","Prestige","Bravery","Sensation","Gender","Situation","Hostility","Village","Profession","Origin","Insolence","Substitute","Obsession","Package","Bird","Greed","Episode","Confidence","Fork","Hindrance","Glory","Majority","Disclosure","Franchise","Outcome","Honey","Fault","Sword","Perjury","Pocket","Fact","Joke","Atmosphere","Painting","Replacement","Amount")
  names(CPW_1tojoin) <- names   
  names(CPW_2tojoin) <- names   
  names(CPW_3tojoin) <- names   
  CPW_all <- rbind(CPW_1tojoin,CPW_2tojoin,CPW_3tojoin)
  
  # add a QA, limiting RT to something in the range 100 ms - 120000 ms
  CPW_all[CPW_all > 120000 | CPW_all < 100] <- NA
  
  CPW_rt <- data.frame(matrix(NA,nrow = ncol(CPW_all),ncol = 2))
  row.names(CPW_rt) <- names(CPW_all)
  names(CPW_rt) <- c("mean","sd")
  
  for (i in 1:ncol(CPW_all)) {
    # SD and mean for each item
    CPW_rt[i,] <- c(mean(CPW_all[,i],na.rm=T),sd(CPW_all[,i],na.rm=T))
  }
}

{ # acc/score
  DDISC_1$DDISC_1_CORR <- rowSums(DDISC_1)
  DDISC_2$DDISC_2_CORR <- rowSums(DDISC_2)
  DDISC_3$DDISC_3_CORR <- rowSums(DDISC_3)
  
  DDISC_all <- c(DDISC_1$DDISC_1_CORR,DDISC_2$DDISC_2_CORR,DDISC_3$DDISC_3_CORR)
  score_sum[6,] <- data.frame(mean = mean(DDISC_all,na.rm = T),sd = sd(DDISC_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  DIGSYM_all <- c(DIGSYM_1$DIGSYM.DSCOR,DIGSYM_3$DIGSYM.DSCOR)
  score_sum[7,] <- data.frame(mean = mean(DIGSYM_all,na.rm = T),sd = sd(DIGSYM_all,na.rm = T))
  
  
  DIGSYMmem_all <- c(DIGSYM_1$DIGSYM.DSMEMCR,DIGSYM_3$DIGSYM.DSMEMCR)
  score_sum[8,] <- data.frame(mean = mean(DIGSYMmem_all,na.rm = T),sd = sd(DIGSYMmem_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  EDISC_1$EDISC_1_CORR <- rowSums(EDISC_1)
  EDISC_2$EDISC_2_CORR <- rowSums(EDISC_2)
  EDISC_3$EDISC_3_CORR <- rowSums(EDISC_3)
  
  EDISC_all <- c(EDISC_1$EDISC_1_CORR,EDISC_2$EDISC_2_CORR,EDISC_3$EDISC_3_CORR)
  score_sum[9,] <- data.frame(mean = mean(EDISC_all,na.rm = T),sd = sd(EDISC_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  ER40_1$ER40_1_CORR <- rowSums(ER40_1)
  ER40_2$ER40_2_CORR <- rowSums(ER40_2)
  ER40_3$ER40_3_CORR <- rowSums(ER40_3)
  
  ER40_all <- c(ER40_1$ER40_1_CORR,ER40_2$ER40_2_CORR,ER40_3$ER40_3_CORR)
  score_sum[10,] <- data.frame(mean = mean(ER40_all,na.rm = T),sd = sd(ER40_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  GNG_1$GNG_1_CORR <- rowSums(GNG_1)
  
  GNG_all <- GNG_1$GNG_1_CORR
  score_sum[11,] <- data.frame(mean = mean(GNG_all,na.rm = T),sd = sd(GNG_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  MEDF_1$MEDF_1_CORR <- rowSums(MEDF_1)
  MEDF_2$MEDF_2_CORR <- rowSums(MEDF_2)
  MEDF_3$MEDF_3_CORR <- rowSums(MEDF_3)
  
  MEDF_all <- c(MEDF_1$MEDF_1_CORR,MEDF_2$MEDF_2_CORR,MEDF_3$MEDF_3_CORR)
  score_sum[12,] <- data.frame(mean = mean(MEDF_all,na.rm = T),sd = sd(MEDF_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  PLOT24_1$PLOT24_1_CORR <- rowSums(PLOT24_1)
  PLOT24_2$PLOT24_2_CORR <- rowSums(PLOT24_2)
  PLOT24_3$PLOT24_3_CORR <- rowSums(PLOT24_3)
  
  PLOT24_all <- c(PLOT24_1$PLOT24_1_CORR,PLOT24_2$PLOT24_2_CORR,PLOT24_3$PLOT24_3_CORR)
  score_sum[13,] <- data.frame(mean = mean(PLOT24_all,na.rm = T),sd = sd(PLOT24_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  PMAT24_1$PMAT24_1_CORR <- rowSums(PMAT24_1)
  PMAT24_2$PMAT24_2_CORR <- rowSums(PMAT24_2)
  PMAT24_3$PMAT24_3_CORR <- rowSums(PMAT24_3)
  
  PMAT24_all <- c(PMAT24_1$PMAT24_1_CORR,PMAT24_2$PMAT24_2_CORR,PMAT24_3$PMAT24_3_CORR)
  score_sum[14,] <- data.frame(mean = mean(PMAT24_all,na.rm = T),sd = sd(PMAT24_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  PVRT25_1$PVRT25_1_CORR <- rowSums(PVRT25_1)
  PVRT25_2$PVRT25_2_CORR <- rowSums(PVRT25_2)
  PVRT25_3$PVRT25_3_CORR <- rowSums(PVRT25_3)
  
  PVRT25_all <- c(PVRT25_1$PVRT25_1_CORR,PVRT25_2$PVRT25_2_CORR,PVRT25_3$PVRT25_3_CORR)
  score_sum[15,] <- data.frame(mean = mean(PVRT25_all,na.rm = T),sd = sd(PVRT25_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  RDISC_1$RDISC_1_CORR <- rowSums(RDISC_1)
  RDISC_2$RDISC_2_CORR <- rowSums(RDISC_2)
  RDISC_3$RDISC_3_CORR <- rowSums(RDISC_3)
  
  RDISC_all <- c(RDISC_1$RDISC_1_CORR,RDISC_2$RDISC_2_CORR,RDISC_3$RDISC_3_CORR)
  score_sum[16,] <- data.frame(mean = mean(RDISC_all,na.rm = T),sd = sd(RDISC_all,na.rm = T))
  
  # rt/speed
  
}

{ # acc/score
  SVOLT_1$SVOLT_1_CORR <- rowSums(SVOLT_1)
  SVOLT_2$SVOLT_2_CORR <- rowSums(SVOLT_2)
  SVOLT_3$SVOLT_3_CORR <- rowSums(SVOLT_3)
  SVOLT_4$SVOLT_4_CORR <- rowSums(SVOLT_4)
  
  SVOLT_all <- c(SVOLT_1$SVOLT_1_CORR,SVOLT_2$SVOLT_2_CORR,SVOLT_3$SVOLT_3_CORR)
  score_sum[17,] <- data.frame(mean = mean(SVOLT_all,na.rm = T),sd = sd(SVOLT_all,na.rm = T))
  
  # rt/speed
  
}



rownames(score_sum) <- c("ADT","AIM","CPF","CPT","CPW","DDISC","DIGSYM","DIGSYM.MEM","EDISC","ER40","GNG","MEDF","PLOT","PMAT","PVRT","RDISC","SVOLT")


write.csv(score_sum,"data/inputs/crowd_sourced_norms/for_ICCs.csv")





