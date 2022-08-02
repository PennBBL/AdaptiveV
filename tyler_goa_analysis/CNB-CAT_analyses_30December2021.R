
library(psych)
library(mirt)

##############################################################################
# latest as of 12January2022
##############################################################################

x <- read.csv("adaptive_commonInterview_data_20220111.csv")

x[x=="Unk"] <- NA

mood <- data.frame(x$bblid,x$dep001_short,x$ocd007,x$gad002,x$dep002_short,x$ocd001,x$sip032,x$pan001,
x$ocd005,x$ocd016,x$pan004,x$gad001,x$pan003,x$dep006_short,x$ocd003,x$scr007,x$man007_short,x$ocd012,
x$dep004,
#x$sui001_short,
x$ocd004,x$ocd006,x$ocd011,x$ocd019,x$ocd013,x$ocd008,x$scr001,x$sip039,
x$ocd018,x$man004_short,x$sip033,x$ocd017,x$man005_short,x$ocd002,x$ocd014,x$scr006,x$ocd015,x$sep510,
x$sip038)

ext <- data.frame(x$bblid,x$odd002,
#x$cdd010,
x$cdd003,x$cdd008,x$cdd005,x$cdd001,x$cdd007,x$odd005,
#x$cdd009,
x$odd001,x$add016_short,x$cdd002,x$odd003,x$add012_short,x$cdd006,x$add011_short,x$odd006,
x$cdd004,x$add021_short,as.numeric(x$add013_short),x$add014_short,x$add022_short,x$add020_short,
as.numeric(x$add015_short),x$cdd011)

psy <- data.frame(x$bblid,x$sip012,x$sip007,x$sip010,x$sip008,x$sip011,x$sip013,x$sip005,x$sip003,
x$sip004,x$sip006,x$sip009,x$psy001,x$psy029,x$sip014,x$psy060,x$psy020,x$psy070,x$man006_short,
x$psy050,x$man003_short,x$man002_short,x$man001_short,x$psy071,x$sip027,x$sip028)

fear <- data.frame(x$bblid,x$agr006,x$agr005,x$soc004,x$agr008,x$soc003,x$agr004,x$agr001,x$soc005,
x$soc001,x$phb004,x$agr002,x$agr003,x$phb007,x$soc002,x$phb006,x$phb001,x$agr007,x$phb002,x$sep509,
x$phb003,x$phb005,x$sep508,x$phb008,x$sep500,x$sep511)

mood[mood==9] <- NA
ext[ext==9] <- NA
psy[psy==9] <- NA
fear[fear==9] <- NA

mood <- mood[rowSums(is.na(mood)) < (ncol(mood)-1),]
ext <- ext[rowSums(is.na(ext)) < (ncol(ext)-1),]
psy <- psy[rowSums(is.na(psy)) < (ncol(psy)-1),]
fear <- fear[rowSums(is.na(fear)) < (ncol(fear)-1),]

mod1 <- mirt(mood[,-1],1)
mod2 <- mirt(ext[,-1],1)
mod3 <- mirt(psy[,-1],1)
mod4 <- mirt(fear[,-1],1)

#mood <- data.frame(mood,rowMeans(mood[,-1],na.rm=TRUE))
#ext <- data.frame(ext,rowMeans(ext[,-1],na.rm=TRUE))
#psy <- data.frame(psy,rowMeans(psy[,-1],na.rm=TRUE))
#fear <- data.frame(fear,rowMeans(fear[,-1],na.rm=TRUE))

mood <- data.frame(mood,scale(fscores(mod1)))
ext <- data.frame(ext,scale(fscores(mod2)))
psy <- data.frame(psy,scale(fscores(mod3)))
fear <- data.frame(fear,scale(fscores(mod4)))

colnames(mood)[ncol(mood)] <- "Mood_Full_Test"
colnames(ext)[ncol(ext)] <- "Externalizing_Full_Test"
colnames(psy)[ncol(psy)] <- "Psychosis_Full_Test"
colnames(fear)[ncol(fear)] <- "Phobias_Full_Test"

# pull in adaptive 

x1 <- read.csv("GACAT_completed_Mood_Anxiety_v1_20220110_0946_clean.csv")
x2 <- read.csv("GACAT_completed_Psychosis_v1_20220110_0948_clean.csv")
x3 <- read.csv("GACAT_completed_Personality_v1_20220110_0950_clean.csv")
x4 <- read.csv("GACAT_completed_Phobias_v1_20220110_0949_clean.csv")
x5 <- read.csv("GACAT_completed_Externalizing_v1_20220110_0947_clean.csv")

colnames(x1) <- paste0(colnames(x1),"_MOOD")
colnames(x2) <- paste0(colnames(x2),"_PSYCHOSIS")
colnames(x3) <- paste0(colnames(x3),"_PERSONALITY")
colnames(x4) <- paste0(colnames(x4),"_PHOBIAS")
colnames(x5) <- paste0(colnames(x5),"_EXTERNALIZING")

temp <- merge(x1[,-c(2:14,16:19)],x2[,-c(2:14,16:19)],by=1,all=TRUE)
temp <- merge(temp,x3[,-c(2:14,16:19)],by=1,all=TRUE)
temp <- merge(temp,x4[,-c(2:14,16:19)],by=1,all=TRUE)
temp <- merge(temp,x5[,-c(2:14,16:19)],by=1,all=TRUE)

x <- merge(mood,temp,by=1,all=TRUE)
x <- merge(ext,x,by=1,all=TRUE)
x <- merge(psy,x,by=1,all=TRUE)
x <- merge(fear,x,by=1,all=TRUE)

x <- x[,c(27,53,77,115,119,117,120,116)]
colnames(x) <- c("Fear_Full","Psychosis_Full","Externalizing_Full","Mood_Full","Fear_Adaptive","Psychosis_Adaptive","Externalizing_Adaptive","Mood_Adaptive")

pairs.panels(x,lm=TRUE,ci=TRUE,scale=TRUE)


































##################################################################
# older, before Kosha/Mrugank created master file used above
##################################################################

x1 <- read.csv("adaptive/adhd_short.csv")
x2 <- read.csv("adaptive/agrophobia.csv")
x3 <- read.csv("adaptive/cgas.csv")
x4 <- read.csv("adaptive/conduct_disorder.csv")
x5 <- read.csv("adaptive/demographics.csv")
x6 <- read.csv("adaptive/depression_short.csv")
x7 <- read.csv("adaptive/eating_disorder.csv")
x8 <- read.csv("adaptive/general_probes.csv")
x9 <- read.csv("adaptive/mania_short.csv")
x10 <- read.csv("adaptive/ocd.csv")
x11 <- read.csv("adaptive/odd.csv")
x12 <- read.csv("adaptive/other_disorder.csv")
x13 <- read.csv("adaptive/overanxiousgad.csv")
x14 <- read.csv("adaptive/panic_disorder.csv")
x15 <- read.csv("adaptive/psychosis.csv")
x16 <- read.csv("adaptive/ptsd_short.csv")
x17 <- read.csv("adaptive/separation_anxiety.csv")
x18 <- read.csv("adaptive/sips_goa_avolition.csv")
x19 <- read.csv("adaptive/sips_goa_expemotion.csv")
x20 <- read.csv("adaptive/sips_goa_focus_attn.csv")
x21 <- read.csv("adaptive/sips_goa_ocpfunctioning.csv")
x22 <- read.csv("adaptive/sips_goa_speech_perc_emotion.csv")
x23 <- read.csv("adaptive/social_anxiety.csv")
x24 <- read.csv("adaptive/suicide_short.csv")
x25 <- read.csv("adaptive/timeline.csv")
#x26 <- read.csv("adaptive/sips_gaf.csv")

x1<- x1[rowSums(is.na(x1)) < max(rowSums(is.na(x1))),-which(names(x1) %in% c("redcapid","assessment","protocol","assessor","date"))]
x2<- x2[rowSums(is.na(x2)) < max(rowSums(is.na(x2))),-which(names(x2) %in% c("redcapid","assessment","protocol","assessor","date"))]
x3<- x3[rowSums(is.na(x3)) < max(rowSums(is.na(x3))),-which(names(x3) %in% c("redcapid","assessment","protocol","assessor","date"))]
x4<- x4[rowSums(is.na(x4)) < max(rowSums(is.na(x4))),-which(names(x4) %in% c("redcapid","assessment","protocol","assessor","date"))]
x5<- x5[rowSums(is.na(x5)) < max(rowSums(is.na(x5))),-which(names(x5) %in% c("redcapid","assessment","protocol","assessor","date"))]
x6<- x6[rowSums(is.na(x6)) < max(rowSums(is.na(x6))),-which(names(x6) %in% c("redcapid","assessment","protocol","assessor","date"))]
x7<- x7[rowSums(is.na(x7)) < max(rowSums(is.na(x7))),-which(names(x7) %in% c("redcapid","assessment","protocol","assessor","date"))]
x8<- x8[rowSums(is.na(x8)) < max(rowSums(is.na(x8))),-which(names(x8) %in% c("redcapid","assessment","protocol","assessor","date"))]
x9<- x9[rowSums(is.na(x9)) < max(rowSums(is.na(x9))),-which(names(x9) %in% c("redcapid","assessment","protocol","assessor","date"))]
x10<- x10[rowSums(is.na(x10)) < max(rowSums(is.na(x10))),-which(names(x10) %in% c("redcapid","assessment","protocol","assessor","date"))]
x11<- x11[rowSums(is.na(x11)) < max(rowSums(is.na(x11))),-which(names(x11) %in% c("redcapid","assessment","protocol","assessor","date"))]
x12<- x12[rowSums(is.na(x12)) < max(rowSums(is.na(x12))),-which(names(x12) %in% c("redcapid","assessment","protocol","assessor","date"))]
x13<- x13[rowSums(is.na(x13)) < max(rowSums(is.na(x13))),-which(names(x13) %in% c("redcapid","assessment","protocol","assessor","date"))]
x14<- x14[rowSums(is.na(x14)) < max(rowSums(is.na(x14))),-which(names(x14) %in% c("redcapid","assessment","protocol","assessor","date"))]
x15<- x15[rowSums(is.na(x15)) < max(rowSums(is.na(x15))),-which(names(x15) %in% c("redcapid","assessment","protocol","assessor","date"))]
x16<- x16[rowSums(is.na(x16)) < max(rowSums(is.na(x16))),-which(names(x16) %in% c("redcapid","assessment","protocol","assessor","date"))]
x17<- x17[rowSums(is.na(x17)) < max(rowSums(is.na(x17))),-which(names(x17) %in% c("redcapid","assessment","protocol","assessor","date"))]
x18<- x18[rowSums(is.na(x18)) < max(rowSums(is.na(x18))),-which(names(x18) %in% c("redcapid","assessment","protocol","assessor","date"))]
x19<- x19[rowSums(is.na(x19)) < max(rowSums(is.na(x19))),-which(names(x19) %in% c("redcapid","assessment","protocol","assessor","date"))]
x20<- x20[rowSums(is.na(x20)) < max(rowSums(is.na(x20))),-which(names(x20) %in% c("redcapid","assessment","protocol","assessor","date"))]
x21<- x21[rowSums(is.na(x21)) < max(rowSums(is.na(x21))),-which(names(x21) %in% c("redcapid","assessment","protocol","assessor","date"))]
x22<- x22[rowSums(is.na(x22)) < max(rowSums(is.na(x22))),-which(names(x22) %in% c("redcapid","assessment","protocol","assessor","date"))]
x23<- x23[rowSums(is.na(x23)) < max(rowSums(is.na(x23))),-which(names(x23) %in% c("redcapid","assessment","protocol","assessor","date"))]
x24<- x24[rowSums(is.na(x24)) < max(rowSums(is.na(x24))),-which(names(x24) %in% c("redcapid","assessment","protocol","assessor","date"))]
x25<- x25[rowSums(is.na(x25)) < max(rowSums(is.na(x25))),-which(names(x25) %in% c("redcapid","assessment","protocol","assessor","date"))]
#x26<- x26[rowSums(is.na(x26)) < max(rowSums(is.na(x26))),-which(names(x26) %in% c("redcapid","assessment","protocol","assessor","date"))]

x <- merge(x1,x2,by="bblid",all=TRUE)
x <- merge(x,x3,by="bblid",all=TRUE)
x <- merge(x,x4,by="bblid",all=TRUE)
x <- merge(x,x5,by="bblid",all=TRUE)
x <- merge(x,x6,by="bblid",all=TRUE)
x <- merge(x,x7,by="bblid",all=TRUE)
x <- merge(x,x8,by="bblid",all=TRUE)
x <- merge(x,x9,by="bblid",all=TRUE)
x <- merge(x,x10,by="bblid",all=TRUE)
x <- merge(x,x11,by="bblid",all=TRUE)
x <- merge(x,x12,by="bblid",all=TRUE)
x <- merge(x,x13,by="bblid",all=TRUE)
x <- merge(x,x14,by="bblid",all=TRUE)
x <- merge(x,x15,by="bblid",all=TRUE)
x <- merge(x,x16,by="bblid",all=TRUE)
x <- merge(x,x17,by="bblid",all=TRUE)
x <- merge(x,x18,by="bblid",all=TRUE)
x <- merge(x,x19,by="bblid",all=TRUE)
x <- merge(x,x20,by="bblid",all=TRUE)
x <- merge(x,x21,by="bblid",all=TRUE)
x <- merge(x,x22,by="bblid",all=TRUE)
x <- merge(x,x23,by="bblid",all=TRUE)
x <- merge(x,x24,by="bblid",all=TRUE)
x <- merge(x,x25,by="bblid",all=TRUE)

x[x=="Unk"] <- NA

mood <- data.frame(
x$bblid,
x$dep001_short,
x$ocd007,
x$gad002,
x$dep002_short,
x$ocd001,
x$sip032,
x$pan001,
x$ocd005,
x$ocd016,
x$pan004,
x$gad001,
x$pan003,
x$dep006_short,
x$ocd003,
x$scr007,
x$man007_short,
x$ocd012,
x$dep004,
x$sui001_short,
x$ocd004,
x$ocd006,
x$ocd011,
x$ocd019,
x$ocd013,
x$ocd008,
x$scr001,
x$sip039,
x$ocd018,
x$man004_short,
x$sip033,
x$ocd017,
x$man005_short,
x$ocd002,
x$ocd014,
x$scr006,
x$ocd015,
x$sep510,
x$sip038)

ext <- data.frame(
x$bblid,
x$odd002,
x$cdd010,
x$cdd003,
x$cdd008,
x$cdd005,
x$cdd001,
x$cdd007,
x$odd005,
x$cdd009,
x$odd001,
x$add016_short,
x$cdd002,
x$odd003,
x$add012_short,
x$cdd006,
x$add011_short,
x$odd006,
x$cdd004,
x$add021_short,
as.numeric(x$add013_short),
x$add014_short,
x$add022_short,
x$add020_short,
as.numeric(x$add015_short),
x$cdd011)

psy <- data.frame(
x$bblid,
#x$sip012,
#x$sip007,
#x$sip010,
#x$sip008,
#x$sip011,
#x$sip013,
#x$sip005,
#x$sip003,
#x$sip004,
#x$sip006,
#x$sip009,
x$psy001,
x$psy029,
#x$sip014,
x$psy060,
x$psy020,
x$psy070,
x$man006_short,
x$psy050,
x$man003_short,
x$man002_short,
x$man001_short,
x$psy071,
x$sip027,
x$sip028)

fear <- data.frame(
x$bblid,
x$agr006,
x$agr005,
x$soc004,
x$agr008,
x$soc003,
x$agr004,
x$agr001,
x$soc005,
x$soc001,
#x$phb004,
x$agr002,
x$agr003,
#x$phb007,
x$soc002,
#x$phb006,
#x$phb001,
x$agr007,
#x$phb002,
x$sep509,
#x$phb003,
#x$phb005,
x$sep508,
#x$phb008,
x$sep500,
x$sep511)

mood[mood==9] <- NA
ext[ext==9] <- NA
psy[psy==9] <- NA
fear[fear==9] <- NA

mood <- mood[rowSums(is.na(mood)) < (ncol(mood)-1),]
ext <- ext[rowSums(is.na(ext)) < (ncol(ext)-1),]
psy <- psy[rowSums(is.na(psy)) < (ncol(psy)-1),]
fear <- fear[rowSums(is.na(fear)) < (ncol(fear)-1),]

mod1 <- mirt(mood[,-1],1)
mod2 <- mirt(ext[,-1],1)
mod3 <- mirt(psy[,-1],1)
mod4 <- mirt(fear[,-1],1)

mood <- data.frame(mood,rowMeans(mood[,-1],na.rm=TRUE))
ext <- data.frame(ext,rowMeans(ext[,-1],na.rm=TRUE))
psy <- data.frame(psy,rowMeans(psy[,-1],na.rm=TRUE))
fear <- data.frame(fear,rowMeans(fear[,-1],na.rm=TRUE))

#temp <- merge(mood,ext,by=1,all=TRUE)
#temp <- merge(temp,psy,by=1,all=TRUE)
#temp <- merge(temp,fear,by=1,all=TRUE)



# pull in adaptive 

x1 <- read.csv("GACAT_completed_Mood_Anxiety_v1_20211028_1921_clean.csv")
x2 <- read.csv("GACAT_completed_Psychosis_v1_20211028_1923_clean.csv")
x3 <- read.csv("GACAT_completed_Personality_v1_20211028_1923_clean.csv")
x4 <- read.csv("GACAT_completed_Phobias_v1_20211028_1923_clean.csv")
x5 <- read.csv("GACAT_completed_Externalizing_v1_20211028_1923_clean.csv")

colnames(x1) <- paste0(colnames(x1),"_MOOD")
colnames(x2) <- paste0(colnames(x2),"_PSYCHOSIS")
colnames(x3) <- paste0(colnames(x3),"_PERSONALITY")
colnames(x4) <- paste0(colnames(x4),"_PHOBIAS")
colnames(x5) <- paste0(colnames(x5),"_EXTERNALIZING")

temp <- merge(x1[,-c(2,3,4,6,7,8,9)],x2[,-c(2,3,4,6,7,8,9)],by=1,all=TRUE)
temp <- merge(temp,x3[,-c(2,3,4,6,7,8,9)],by=1,all=TRUE)
temp <- merge(temp,x4[,-c(2,3,4,6,7,8,9)],by=1,all=TRUE)
temp <- merge(temp,x5[,-c(2,3,4,6,7,8,9)],by=1,all=TRUE)


x <- merge(mood,temp,by=1,all=TRUE)
x <- merge(ext,x,by=1,all=TRUE)
x <- merge(psy,x,by=1,all=TRUE)
x <- merge(fear,x,by=1,all=TRUE)

x <- x[,c(19,33,59,98,102,100,103,99)]
colnames(x) <- c("Fear_Full","Psychosis_Full","Externalizing_Full","Mood_Full","Fear_Adaptive","Psychosis_Adaptive","Externalizing_Adaptive","Mood_Adaptive")

pairs.panels(x,lm=TRUE,ci=TRUE,scale=TRUE)




#hist(rowSums(is.na(x1)) / rowSums(is.na(x1)==FALSE))