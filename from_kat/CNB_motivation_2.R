#Libraries
library(psych)
library(GPArotation)
library(tidyverse)
library(corrplot)
library(reshape)
library(grid)
library(ggplot2)
library(gtools)
library(lubridate)
library(skimr)
library(janitor)
library (flextable)
library(cowplot)
library(matrixStats)

dat= read.csv("from_kat/cnb_withgroups.csv") #mostuotdate dataset

### Examine validity and remove invalid data no invalid data ----
Validity_subset= dat %>%
  dplyr::select(ends_with("_valid"))

#removing duplicate KRDISC
dat= dat[, -c(534:615)]

### DISC--> create mean score for accuracy & median for RT ----
dat_KDDISC= dat[ , 111:178] %>%dplyr::select(starts_with("KDDISC.q_"))
dat_KRDISC= dat[ , 180:261] %>%dplyr::select(starts_with("KRDISC.q_"))
dat_EDISC= dat[ , 463:530] %>%dplyr::select(ends_with("_resp"))
dat$KDDISC_mean_acc= rowMeans(dat_KDDISC, na.rm = TRUE)
dat$KRDISC_mean_acc= rowMeans(dat_KRDISC, na.rm=TRUE)
dat$EDISC_mean_acc= rowMeans(dat_EDISC, na.rm=TRUE) #excludes the first 100


dat_KDDISC_RT= dat[ , 111:178] %>%dplyr::select(starts_with("KDDISC.trr_"))
dat_KRDISC_RT= dat[ , 180:261] %>%dplyr::select(starts_with("KRDISC.trr_"))
dat_EDISC_RT= dat[ , 463:530] %>%dplyr::select(ends_with("_ttr"))
dat$KDDISC_med_RT= rowMedians(as.matrix(dat_KDDISC_RT, na.rm=TRUE)) #, na.rm=TRUE???? check 
dat$KRDISC_med_RT= rowMedians(as.matrix(dat_KRDISC_RT, na.rm=TRUE))
dat$EDISC_med_RT= rowMedians(as.matrix(dat_EDISC_RT, na.rm=TRUE))



### look at histogram for each variable----
dat_scores_acc= dat[, c(38,50,56,59,63,68,73,77,81,85,93,94,99,103,534:536)]
dat_scores_rt= dat[, c(39,51,54,60,64,69,74,78,82,86,95:96,100,104,537:539)]



pdf("CNBacc_histograms_101022.pdf",height=3,width=3)
colNames= names(dat_scores_acc)
for(i in colNames){
  plt <- ggplot(dat_scores_acc, aes_string(x=i)) +
    geom_histogram(bins=30)+
    stat_count(width = 0.5)
  print(plt)
  Sys.sleep(2)
}
dev.off()

pdf("CNB_RT_histograms_101022.pdf",height=3,width=3)
colNames= names(dat_scores_rt)
for(i in colNames){
  plt <- ggplot(dat_scores_rt, aes_string(x=i)) +
    geom_histogram(bins=30)+
    stat_count(width = 0.5)
  print(plt)
  Sys.sleep(2)
}
dev.off()

plt <- print(ggplot(dat, aes(x=dscor)) +
               geom_histogram(bins=30)+
               stat_count(width = 0.5))
dat= dat [, -c(111:178, 180:261, 263:533)]

#Tyler's comments for revising histograms for accuracy

dat$cpt_ptp <- ifelse(dat$cpt_ptp <0.5, NA, dat$cpt_ptp)
dat$cpt_ptp <- winsor(dat$cpt_ptp,trim=0.025)

dat$er40_cr <- ifelse(dat$er40_cr <20, NA, dat$er40_cr)
dat$er40_cr <- winsor(dat$er40_cr,trim=0.025)

dat$pvrt_cr <- ifelse(dat$pvrt_cr ==0, NA, dat$pvrt_cr) #code the cpt_ptp at 0.50 or less as missing and then winsorize the rest at 0.025
dat$pvrt_cr <- winsor(dat$pvrt_cr,trim=0.025)

dat$pmat_pc <- ifelse(dat$pmat_pc <20, NA, dat$pvrt_cr) #code the cpt_ptp at 0.50 or less as missing and then winsorize the rest at 0.025
dat$pmat_pc <- winsor(dat$pmat_pc,trim=0.025)

dat$volt_cr <- ifelse(dat$volt_cr <=10, NA, dat$volt_cr) #anyone with 10 or fewer should be coded missing 
dat$volt_cr <- winsor(dat$volt_cr,trim=0.025)

dat$plot_pc <- winsor(dat$plot_pc,trim=0.025)

dat$cpf_cr <- ifelse(dat$cpf_cr <=20, NA, dat$cpf_cr) #anyone with 20 or fewer should be coded missing
dat$cpf_cr <- winsor(dat$cpf_cr,trim=0.025)

dat$medf_pc <- ifelse(dat$medf_pc <33, NA, dat$medf_pc)#medf_pc - anyone < 33 should be coded missing
dat$medf_pc <- winsor(dat$medf_pc,trim=0.025)

dat$adt_pc <- ifelse(dat$adt_pc <33, NA, dat$adt_pc)#adt_pc - anyone < 33 should be coded missing
dat$adt_pc <- winsor(dat$adt_pc,trim=0.025)

dat$cpw_cr <- ifelse(dat$cpw_cr <=20, NA, dat$cpw_cr) #anyone with 20 or fewer should be coded missing
dat$cpw_cr <- winsor(dat$cpw_cr,trim=0.025)

dat$dscor <- ifelse(dat$dscor >100, NA, dat$dscor)# dscor - anyone >150 should be coded missing + the other top responders (>90) should be winsorized to 90
          temp <- dat[,93]  
          temp[temp > 90] <- 90
          dat[,93] <- temp
dat$dscor <- winsor(dat$dscor,trim=0.025)

dat$gng_cr <- ifelse(dat$gng_cr <100, NA, dat$gng_cr)# gng_cr - anyone <100 should be coded missing + the remaining < 120 should be winsorized to 120
          temp1 <- dat[,99]  
          temp1[temp1 < 120] <- 120
          dat[,99] <- temp1
dat$gng_cr <- winsor(dat$gng_cr,trim=0.025)
          
dat$aim_tot <- ifelse(dat$aim_tot <=30, NA, dat$aim_tot)#aim_tot - anyone 30 or below 
dat$aim_tot <- winsor(dat$aim_tot,trim=0.025)

dat$KDDISC_mean_acc <- winsor(dat$KDDISC_mean_acc,trim=0.025)
dat$KRDISC_mean_acc <- winsor(dat$KRDISC_mean_acc,trim=0.025)
dat$EDISC_mean_acc <- winsor(dat$EDISC_mean_acc,trim=0.025)




#Tyler's comments for revising histograms for Speed
#One is to remove responses < 1000ms (code missing) and winsorize the rest at 2.5%.  
#I'd apply this rule to: PMAT, PLOT, and all DISC tasks
dat$pmat_rtcr <- ifelse(dat$pmat_rtcr <1000, NA, dat$pmat_rtcr)
dat$pmat_rtcr <- winsor(dat$pmat_rtcr,trim=0.025)

dat$plot_rtcr <- ifelse(dat$plot_rtcr <1000, NA, dat$plot_rtcr)
dat$plot_rtcr <- winsor(dat$plot_rtcr,trim=0.025)

dat$KDDISC_med_RT <- ifelse(dat$KDDISC_med_RT <1000, NA, dat$KDDISC_med_RT)
dat$KDDISC_med_RT <- winsor(dat$KDDISC_med_RT,trim=0.025)

dat$KRDISC_med_RT <- ifelse(dat$KRDISC_med_RT <1000, NA, dat$KRDISC_med_RT)
dat$KRDISC_med_RT <- winsor(dat$KRDISC_med_RT,trim=0.025)

dat$EDISC_med_RT <- ifelse(dat$EDISC_med_RT <1000, NA, dat$EDISC_med_RT)
dat$EDISC_med_RT <- winsor(dat$EDISC_med_RT,trim=0.025)

dat$pvrt_rtcr <- ifelse(dat$pvrt_rtcr <1000, NA, dat$pvrt_rtcr)
dat$pvrt_rtcr <- winsor(dat$pvrt_rtcr,trim=0.025)

#The other rule is to remove responses < 500ms (code missing) and winsorize the rest at 2.5% 
#I'd apply this to: ER40, VOLT, CPF, MEDF, ADT, CPW, digit-symbol, digit-symbol memory, and AIM
dat$er40_rtcr <- ifelse(dat$er40_rtcr <500, NA, dat$er40_rtcr)
dat$er40_rtcr <- winsor(dat$er40_rtcr,trim=0.025)

dat$volt_rtcr <- ifelse(dat$volt_rtcr <500, NA, dat$volt_rtcr)
dat$volt_rtcr <- winsor(dat$volt_rtcr,trim=0.025)

dat$cpf_rtcr <- ifelse(dat$cpf_rtcr <500, NA, dat$cpf_rtcr)
dat$cpf_rtcr <- winsor(dat$cpf_rtcr,trim=0.025)

dat$medf_rtcr <- ifelse(dat$medf_rtcr <500, NA, dat$medf_rtcr)
dat$medf_rtcr <- winsor(dat$medf_rtcr,trim=0.025)

dat$adt_rtcr <- ifelse(dat$adt_rtcr <500, NA, dat$adt_rtcr)
dat$adt_rtcr <- winsor(dat$adt_rtcr,trim=0.025)

dat$cpw_rtcr <- ifelse(dat$cpw_rtcr <500, NA, dat$cpw_rtcr)
dat$cpw_rtcr <- winsor(dat$cpw_rtcr,trim=0.025)

dat$dscorrt <- ifelse(dat$dscorrt <500, NA, dat$dscorrt)
dat$dscorrt <- winsor(dat$dscorrt,trim=0.025)

dat$dsmcrrt <- ifelse(dat$dsmcrrt <500, NA, dat$dsmcrrt)
dat$dsmcrrt <- winsor(dat$dsmcrrt,trim=0.025)

dat$aim_totrt <- ifelse(dat$aim_totrt <500, NA, dat$aim_totrt)
dat$aim_totrt <- winsor(dat$aim_totrt,trim=0.025)


#Check histos again
dat_scores_acc= dat[, c(38,50,56,59,63,68,73,77,81,85,93:94,99,103,113:115)]
dat_scores_rt= dat[, c(39,51,54,60,64,69,74,78,82,86,95,96,100,104,116:118)]

pdf("CNBacc_histograms_092822.pdf",height=3,width=3)
colNames= names(dat_scores_acc)
for(i in colNames){
  plt <- ggplot(dat_scores_acc, aes_string(x=i)) +
    geom_histogram(bins=30)+
    stat_count(width = 0.5)
  print(plt)
  Sys.sleep(2)
}
dev.off()

pdf("CNB_RT_histograms_092822.pdf",height=3,width=3)
colNames= names(dat_scores_rt)
for(i in colNames){
  plt <- ggplot(dat_scores_rt, aes_string(x=i)) +
    geom_histogram(bins=30)+
    stat_count(width = 0.5)
  print(plt)
  Sys.sleep(2)
}
dev.off()


##### regress DScor out of the DSmem ----
mod1 <- lm(dsmemcr~dscor,data=dat_scores_acc,na.action=na.exclude)
DS_acc_reg <- residuals(mod1,na.action=na.exclude) + mean(dat_scores_acc$dsmemcr,na.rm=TRUE)
dat_scores_acc$DS_mem_acc_reg=DS_acc_reg
dat_scores_acc=dat_scores_acc[, -c(12)]

mod2 <- lm(dsmcrrt~dscorrt,data=dat_scores_rt,na.action=na.exclude)
DS_rt_reg <- residuals(mod2,na.action=na.exclude) + mean(dat_scores_rt$dsmcrrt,na.rm=TRUE)
dat_scores_rt$DS_mem_rt_reg=DS_rt_reg
dat_scores_rt= dat_scores_rt[, -c(12)]


# check quantity of variables
describe(dat_scores_acc)
describe(dat_scores_rt)


##### create linear and nonlinear age variables (age, age2, age3)
age <- dat$age_enroll
age_squared <- scale(age)^2
age_cubed <- scale(age)^3


##### regress age out (age, age2, age3)
colNames= names(dat_scores_acc)
newdat_acc_age_reg <- matrix(NA,nrow(dat_scores_acc),ncol(dat_scores_acc))
for (i in 1:length(colNames)) {
  mod<- lm(dat_scores_acc[,i]~ age + age_squared + age_cubed, na.action=na.exclude)
  newdat_acc_age_reg[,i] <- scale(residuals(mod, na.action=na.exclude))
}
colnames(newdat_acc_age_reg) <- colNames
newdat_acc_age_reg= data.frame(newdat_acc_age_reg)
pairs.panels(newdat_acc_age_reg)

colNames_RT= names(dat_scores_rt)
newdat_rt_age_reg <- matrix(NA,nrow(dat_scores_rt),ncol(dat_scores_rt))
for (i in 1:length(colNames_RT)) {
  mod<- lm(dat_scores_rt[,i]~ age + age_squared + age_cubed, na.action=na.exclude)
  newdat_rt_age_reg[,i] <- scale(residuals(mod, na.action=na.exclude))
}
colnames(newdat_rt_age_reg) <- colNames_RT
newdat_rt_age_reg= data.frame(newdat_rt_age_reg)
pairs.panels(dat_scores_rt)


#### EFA on Efficiency----
accuracy_z <- scale(newdat_acc_age_reg)
RT_z <- scale(newdat_rt_age_reg)
efficiency <- scale(accuracy_z - RT_z)
efficiency = data.frame(efficiency)

o= cor(efficiency, use = "pairwise")
corrplot(o ,order = "hclust", tl.srt=45, method= 'number',
         tl.pos = 'd', number.cex=.55,cl.cex=1.1, tl.col="black", tl.cex =.56,
         pointsize=25)
VSS.scree(o)
fa.parallel(efficiency)# the number of factors =  4
nfactors(efficiency) #Empirical BIC= 2, MAP=1

fa.diagram(fa(efficiency,2,rotate="promax"))
fa.diagram(fa(efficiency,3,rotate="promax"))
fa.diagram(fa(efficiency,4,rotate="promax"))
fa.diagram(fa(efficiency,5,rotate="promax"))
fa.diagram(fa(efficiency,6,rotate="promax"),digits=2,cut=0.2)
m=fa.sort(fa(efficiency,6,rotate="promax"))

omega.diagram(omega(efficiency,5,rotate="promax", plot=FALSE),digits=2,
              cut=0.25,gcut=0.0001)  

write.csv(efficiency, "efficiency_CNB.csv")


#### FA and correlation Accuracy ----
m= cor(newdat_acc_age_reg, use = "pairwise")
corrplot(m ,order = "hclust", tl.srt=45,  addrect = 2, method= 'number',
         tl.pos = 'd', number.cex=.55,cl.cex=1.1, tl.col="black", tl.cex =.56,
         pointsize=25)
VSS.scree(m) #4
fa.parallel(newdat_acc_age_reg) # the number of factors =  3, components=2
nfactors(newdat_acc_age_reg) # Empirical BIC suggest= 2; MAP= 2

# fa.sort(fa(newdat_acc_age_reg,5,rotate="promax"))

fa.diagram(fa(newdat_acc_age_reg,2,rotate="promax"),cut=0.2)
fa.sort(fa(newdat_acc_age_reg,5,rotate="promax"))
fa.diagram(fa(newdat_acc_age_reg,3,rotate="promax"),cut=0.2)
fa.diagram(fa(newdat_acc_age_reg,4,rotate="promax"),cut=0.2)
fa.diagram(fa(newdat_acc_age_reg,5,rotate="promax"),cut=0.2)
fa.diagram(fa(newdat_acc_age_reg,6,rotate="promax"),cut=0.2)

omega.diagram(omega(newdat_acc_age_reg,5,rotate="promax", plot=FALSE),digits=2,
              cut=0.25,gcut=0.0001)  



##### RT ----
n= cor(newdat_rt_age_reg, use = "pairwise")
corrplot(n,order = "hclust", tl.srt=45,  addrect = 1, method= 'number',
         tl.pos = 'd', number.cex=.55,cl.cex=1.1, tl.col="black", tl.cex =.56,
         pointsize=25)
VSS.scree(n)
fa.parallel(newdat_rt_age_reg) # the number of factors =  4
nfactors(newdat_rt_age_reg) # Empirical BIC suggest= 4; MAP= 1


fa.diagram(fa(newdat_rt_age_reg,2,rotate="promax"),digits=2,cut=0.2)
fa.diagram(fa(newdat_rt_age_reg,3,rotate="promax"),digits=2,cut=0.2)
fa.diagram(fa(newdat_rt_age_reg,4,rotate="promax"),digits=2,cut=0.2)
fa.diagram(fa(newdat_rt_age_reg,5,rotate="promax"),digits=2,cut=0.2)
fa.diagram(fa(newdat_rt_age_reg,6,rotate="promax"),digits=2,cut=0.2)


omega.diagram(omega(newdat_rt_age_reg,5,rotate="promax", plot=FALSE),digits=2,
              cut=0.25,gcut=0.0001)





#### CFA in Lavaan ----
library(sem)
library(lavaan)
library(missForest)
library(semPlot)

set.seed(2022)
#Imputation of matrix using MissForst
efficiency[,1:17] <- missForest(efficiency[,1:17])$ximp

CNB.Fullmodel <- 
  ' MR1  =~ gng_cr + cpt_ptp  
              MR2  =~ pvrt_cr + pmat_pc  + plot_pc  + aim_tot + dscor
              MR3  =~ KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc 
              MR4  =~ cpw_cr + volt_cr + DS_mem_acc_reg 
              MR5  =~ adt_pc + medf_pc + er40_cr + cpf_cr'


CNB.Fullmodel2 <- 
  ' MR1  =~ pvrt_cr + pmat_pc  + plot_pc
              MR2  =~ KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc
              MR3  =~ aim_tot + gng_cr + cpt_ptp  + dscor
              MR4  =~ cpf_cr + cpw_cr + volt_cr + DS_mem_acc_reg 
              MR5  =~ adt_pc + medf_pc + er40_cr'

CNB.Fullmodel3 <- 
  ' MR1  =~ pvrt_cr + pmat_pc  + plot_pc + aim_tot + dscor
              MR2  =~ KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc
              MR3  =~ gng_cr + cpt_ptp
              MR4  =~ cpw_cr + volt_cr + DS_mem_acc_reg 
              MR5  =~ adt_pc + medf_pc + er40_cr + cpf_cr
              G  =~ pvrt_cr + pmat_pc  + plot_pc+
                      KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc+
                      aim_tot + gng_cr + cpt_ptp  + dscor+
                      cpf_cr + cpw_cr + volt_cr + DS_mem_acc_reg+
                      adt_pc + medf_pc + er40_cr
                      
                # fix covariance between f1 and f2 to zero
                MR1~~0*MR2
                MR1~~0*MR3
                MR1~~0*MR4
                MR1~~0*MR5
                MR1~~ 0*G
                MR2~~0*MR3
                MR2~~0*MR4
                MR2~~0*MR5
                MR2~~ 0*G
                MR3~~0*MR4
                MR3~~0*MR5
                MR3~~ 0*G
                MR4~~0*MR5
                MR4~~ 0*G
                MR5~~ 0*G'


CNB.Fullmodel4 <- 
  ' MR1  =~ pvrt_cr + pmat_pc  + plot_pc + aim_tot + dscor
              MR2  =~ KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc
              MR3  =~ v2*gng_cr + v2*cpt_ptp
              MR4  =~ cpw_cr + volt_cr + DS_mem_acc_reg 
              MR5  =~ adt_pc + medf_pc + er40_cr + cpf_cr
              G  =~ pvrt_cr + pmat_pc  + plot_pc+
                      KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc+
                      aim_tot + gng_cr + cpt_ptp  + dscor+
                      cpf_cr + cpw_cr + volt_cr + DS_mem_acc_reg+
                      adt_pc + medf_pc + er40_cr
                      
                # fix covariance between f1 and f2 to zero
                MR1~~0*MR2
                MR1~~0*MR3
                MR1~~0*MR4
                MR1~~0*MR5
                MR1~~ 0*G
                MR2~~0*MR3
                MR2~~0*MR4
                MR2~~0*MR5
                MR2~~ 0*G
                MR3~~0*MR4
                MR3~~0*MR5
                MR3~~ 0*G
                MR4~~0*MR5
                MR4~~ 0*G
                MR5~~ 0*G'

CNB.Fullmodel5 <- 
  ' MR1  =~ pvrt_cr + pmat_pc  + plot_pc + aim_tot + dscor
              MR2  =~ KDDISC_mean_acc + KRDISC_mean_acc + EDISC_mean_acc
              MR3  =~ v2*gng_cr + v2*cpt_ptp
              MR4  =~ cpw_cr + volt_cr + DS_mem_acc_reg 
              MR5  =~ adt_pc + medf_pc + er40_cr + cpf_cr
                # fix covariance between f1 and f2 to zero
                MR1~~0*MR2
                MR1~~0*MR3
                MR1~~0*MR4
                MR1~~0*MR5
                MR2~~0*MR3
                MR2~~0*MR4
                MR2~~0*MR5
                MR3~~0*MR4
                MR3~~0*MR5
                MR4~~0*MR5
'


# fit the model
fit <- lavaan::cfa(CNB.Fullmodel, data = efficiency)
fit <- lavaan::cfa(CNB.Fullmodel2, data = efficiency)
fit3 <- lavaan::cfa(CNB.Fullmodel3, data = efficiency)
fit4 <- lavaan::cfa(CNB.Fullmodel4, data = efficiency)
fit5 <- lavaan::cfa(CNB.Fullmodel5, data = efficiency)

# display summary output
# CFI should be > .90.
# RMSEA should be < .08 or < .05
# SRMR should be < .08

summary(fit2, fit.measures = TRUE)
standardizedSolution(fit) #provides CI of estimates

Original_efficiency= fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"))

print(fitMeasures(fit5, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                  output = "text"), add.h0 = TRUE)


# bootstrap to get bootstrap test statistics
bootstrap_efficiency <- bootstrapLavaan(fit, R=1000, type="ordinary",
                                        FUN=fitMeasures, fit.measures=c("chisq", "cfi", "tli", "rmsea", "srmr"))

# compute a bootstrap based p-value
pvalue.boot <- length(which(bootstrap_efficiency > Original_efficiency))/length(bootstrap_efficiency)


hist(bootstrap_efficiency[,2])
hist(bootstrap_efficiency[,3])
hist(bootstrap_efficiency[,4])
hist(bootstrap_efficiency[,5])


quantile(bootstrap_efficiency[,2],0.05)
quantile(bootstrap_efficiency[,3],0.05)
quantile(bootstrap_efficiency[,4],0.95)
quantile(bootstrap_efficiency[,5],0.95)



semPaths(fit4, edge.label.cex=.85,
         whatLabels= "std", style= "lisrel",curve = 1.8,
         layout = "tree2", curvePivot = TRUE, rotation=2, 
         nCharNodes=6,sizeMan=7,sizeInt=5, label.prop = .95,
         what= "col", residuals= FALSE, intercepts= FALSE,curveAdjacent = TRUE,
         fixedStyle = c("darkblue",3), optimizeLatRes=TRUE)


