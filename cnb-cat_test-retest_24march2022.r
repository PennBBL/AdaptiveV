
library(psych)
library(qgraph)

x <- read.csv("data/outputs/full_CAT_cnb_23Mar22.csv")

###############################################################
# x <- x[which(x$group == "Healthy Controls"),]               # limit to controls?
###############################################################

temp <- x$gng_cr
temp[temp<100] <- NA
x$gng_cr <- temp

temp <- x$GNG60.GNG60_CR
temp[temp<50] <- NA
x$GNG60.GNG60_CR <- temp

cpt_acc <- x$cpt_ptp - x$cpt_pfp

acc <- x[,c(28,34,37,41,46,51,55,59,63,67,68,73,77)]

er40_cat <- ((x$er40.1.00.cat_neutral + x$er40.1.00.cat_emotive + x$er40.1.00.cat_emotive + x$er40.1.00.cat_emotive)/4)*sqrt(4)
medf_cat <- ((x$medf.1.00.cat_same + x$medf.1.00.cat_different + x$medf.1.00.cat_different + x$medf.1.00.cat_different)/4)*sqrt(4)
adt_cat <- ((x$adt.1.00.cat_same + x$adt.1.00.cat_different + x$adt.1.00.cat_different + x$adt.1.00.cat_different)/4)*sqrt(4)

new_er40_cat <- ((x$er40.1.00.cat_neutral)*(2/5) + (x$er40.1.00.cat_emotive)*(3/5))*sqrt(5)
new_medf_cat <- ((x$medf.1.00.cat_same)*(5/17) + (x$medf.1.00.cat_different)*(12/17))*sqrt(17)
new_adt_cat <- ((x$adt.1.00.cat_same)*(5/17) + (x$adt.1.00.cat_different)*(12/17))*sqrt(17)

cpf_cat <- ((x$cpf.1.00.v1.cat_target + x$cpf.1.00.v1.cat_foil)/2)*sqrt(2)
cpw_cat <- ((x$cpw.1.00.v1.cat_target + x$cpw.1.00.v1.cat_foil)/2)*sqrt(2)
volt_cat <- ((x$volt.1.00.v1.cat_targets + x$volt.1.00.v1.cat_foils)/2)*sqrt(2)
cpt_cat <- (x$CPT108.CATCPTT_TP/36) - (x$CPT108.CATCPTL_FP/72)

cat_acc <- x[,c(96,103,108:111,136,142,151,161)]

x99 <- data.frame(x[,c(6,7,13,84:93)],acc,cpt_acc,cat_acc,er40_cat,medf_cat,adt_cat,cpf_cat,cpw_cat,volt_cat,cpt_cat)

sc <- matrix(NA,nrow(x99),31)

for (i in 1:31) {
mod <- lm(x99[,(i+13)]~proto_3,data=x99,na.action=na.exclude)
sc[,i] <- scale(residuals(mod,na.action=na.exclude))
}

colnames(sc) <- paste0(colnames(x99[,14:44]),"_Oreg")

x <- data.frame(x99,sc)

# write.csv(x,"CNB-CAT_test-retest_with_order-regressed.csv",na="",row.names=FALSE)

pdf("data/outputs/scatters/oldCNB-CAT_test-retest_scatter_matrices_CONTROLS050522.pdf",height=9,width=12)
pairs.panels(x[,c(46,47,52,59,61,60)],lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x[,c(45,50,51,69,70,71)],lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x[,c(48,49,53,55,74,72,73,68)],lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x[,c(54,56,57,58,67,66,65,75)],lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()













# old

xcor <- cor_auto(x99[,4:36])
# write.csv(xcor,"CNB-CAT_test-retest_validity_cor.csv")



fa.sort(fa(x99[,6:19]))

fa.sort(fa(x99[,c(20:22,26:36)],4,fm="ml"))


