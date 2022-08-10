# Extralong test-retest as part of AdaptiveV

# Akira Di Sandro, 7.21.22


`%notin%` <- Negate(`%in%`)

# load packages ----
library(psych)
library(qgraph)
library(dplyr)
library(PerFit)
library(ggplot2)
library(ggdist)
library(sdamr)
library(kableExtra)
library(matrixStats)
library(irr)
library(lubridate)
# library(magick)



# load CSVs and create dataset ----

XL <- read.csv("data/inputs/cnb_merged_webcnp_surveys_allbblprjcts_longform.csv")  #29895 rows , 8/4/22
extralong <- XL %>% filter(test_sessions.bblid.clean>9999) %>% rename(bblid = test_sessions.bblid.clean)  # 16093 rows
# extralong <- XL %>% filter(test_sessions.bblid.clean>9999,test_sessions.siteid != "adaptive_v") %>% rename(bblid = test_sessions.bblid.clean)

# cpt_acc generation
extralong <- extralong %>% mutate(cpt_acc = cpt_ptp - cpt_pfp) %>% dplyr::select(datasetid_platform:cpt_pfp,cpt_acc,cpt_fprt:KRDISC.trr_41.1)

# disc scoring
ddisc_qs <- extralong %>% dplyr::select((names(extralong)[grep('DDISC.q',colnames(extralong))]))
extralong$ddisc_sum <- rowSums(ddisc_qs-1)
ddisc_ttrs <- extralong %>% dplyr::select((names(extralong)[grep('DDISC.trr',colnames(extralong))]))
extralong$ddisc_mcr <- rowMedians(as.matrix(ddisc_ttrs))

rdisc_qs <- extralong %>% dplyr::select((names(extralong)[grep('RDISC.q',colnames(extralong))]))
extralong$rdisc_sum <- rowSums(rdisc_qs[,1:41]-1)
rdisc_ttrs <- extralong %>% dplyr::select((names(extralong)[grep('RDISC.trr',colnames(extralong))]))
extralong$rdisc_mcr <- rowMedians(as.matrix(rdisc_ttrs[,1:41]))

edisc_qs <- extralong %>% dplyr::select((names(extralong)[grep('EDISC',colnames(extralong))])) %>% 
  dplyr::select((names(extralong)[grep('_resp',colnames(extralong))])) %>% 
  dplyr::select(EDISC.q_101_resp:EDISC.q_134_resp)
extralong$edisc_sum <- rowSums(edisc_qs-1)
edisc_ttrs <- extralong %>% dplyr::select((names(extralong)[grep('EDISC',colnames(extralong))])) %>% 
  dplyr::select((names(extralong)[grep('_ttr',colnames(extralong))])) %>% 
  dplyr::select(EDISC.q_101_ttr:EDISC.q_134_ttr)
extralong$edisc_mcr <- rowMedians(as.matrix(edisc_ttrs))

extralong <- extralong %>% dplyr::select(datasetid_platform:KDDISC.valid_code,KDDISC.test,ddisc_sum,ddisc_mcr,
                                         KRDISC.valid_code,KRDISC.test,rdisc_sum,rdisc_mcr,
                                         EDISC.valid_code,EDISC.test,edisc_sum,edisc_mcr)


# keep only those with multiple timepoints
extralong_repeat <- extralong %>% group_by(bblid) %>% dplyr::summarize(n = n()) %>% 
  left_join(extralong,.,by="bblid") %>% filter(n>1,test_sessions.datasetid %notin% c(48505,42470)) %>% 
  arrange(bblid,test_sessions_v.dotest)

extralong_repeat$test_sessions_v.dotest <- as.Date(extralong_repeat$test_sessions_v.dotest,format = "%Y-%m-%d")
extralong_repeat$days <- as.numeric(ymd(extralong_repeat$test_sessions_v.dotest))

# difference between timepoints in days
extralong_repeat$diffdays <- 0
for (i in 2:nrow(extralong_repeat)) {
  extralong_repeat$diffdays[i] <- extralong_repeat$days[i] - extralong_repeat$days[i-1]
}

temp <- extralong_repeat %>% dplyr::select("bblid","test_sessions_v.dotest")
t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")

extralong_repeat$combo <- paste(as.character(extralong_repeat$bblid),as.character(extralong_repeat$test_sessions_v.dotest),sep="_")
extralong_repeat$diffdays <- ifelse(extralong_repeat$combo %in% t.first$combo,0,extralong_repeat$diffdays)
# need to fix 92134_webcnp (there were two sessions in the same day so dif days comes out as 0 as well)
extralong_repeat[which(extralong_repeat$datasetid_platform == "92134_webcnp"),"diffdays"] <- extralong_repeat %>% filter(datasetid_platform == "92127_webcnp") %>% dplyr::select(diffdays)

# keep only non-timepoint1 datapoints to look at histograms
repeats_only <- extralong_repeat %>% filter(diffdays > 0) %>% mutate(constant = 1)  # 3265 rows as of 7/27/22


# PRA from itemwise
for_pra <- read.csv("data/inputs/athena_195_360.csv",na.strings=c(""," ","NA"))
for_pra2 <-  read.csv("data/inputs/athena_253_324.csv",na.strings=c(""," ","NA"))

PRA_iw1 <- for_pra %>% mutate(bblid = as.numeric(test_sessions_v.bblid)) %>% arrange(bblid) %>% 
  filter(bblid > 9999,!is.na(PRA_D.PRADWORDCR)) %>% 
  dplyr::select(matches("test_session|^bblid|PRA_D"))

PRA_iw1$PRA_D.PRADWORDCR <- ifelse(!is.na(PRA_iw1$PRA_D.PRADWORDCR),PRA_iw1$PRA_D.PRADWORDCR,0)

PRA_iw2 <- for_pra2 %>% mutate(bblid = as.numeric(test_sessions.bblid)) %>% arrange(bblid) %>% 
  filter(bblid > 9999,!is.na(PRA_D.PRADWORDCR)) %>% 
  dplyr::select(matches("test_session|^bblid|PRA_D"))

PRA_iw2_tomerge <- PRA_iw2 %>% dplyr::select(bblid,test_sessions.datasetid:test_sessions.famid,test_sessions_v.battery,
                                             test_sessions_v.dotest,test_sessions_v.valid_code,test_sessions_v.age,
                                             test_sessions_v.dob,test_sessions_v.education,PRA_D.test:PRA_D.PRADCR_RAW)

# all(PRA_iw$test_sessions.datasetid  %in% PRA_iw2$test_sessions.datasetid) # TRUE, so we can just stick with PRA_iw2

# athena_254_360 <- read.csv("data/inputs/athena_254_360_220713.csv",na.strings=c(""," ","NA")) # PRA cols, but empty
# athena_3360_1878 <- read.csv("data/inputs/athena_3360_1878.csv",na.strings=c(""," ","NA")) # no PRA cols
athena_3360_2096 <- read.csv("data/inputs/athena_3360_2096_220713.csv",na.strings=c(""," ","NA"))

PRA_iw3 <- athena_3360_2096 %>% mutate(bblid = as.numeric(test_sessions_v.bblid)) %>% arrange(bblid) %>% 
  filter(bblid > 9999,!is.na(PRA_D.PRADWORDCR)) %>% 
  dplyr::select(matches("test_session|^bblid|PRA_D"))

PRA_iw3_tomerge <- PRA_iw3 %>% mutate(test_sessions_v.valid_code = NA,PRA_D.test = NA,PRA_D.valid_code = NA) %>% 
  dplyr::select(bblid,test_sessions.datasetid:test_sessions.famid,test_sessions_v.battery,test_sessions_v.dotest,
                test_sessions_v.valid_code,test_sessions_v.age,test_sessions_v.dob,test_sessions_v.education,
                PRA_D.test,PRA_D.valid_code,PRA_D.AGE_MON,PRA_D.AGE_YR,PRA_D.PRADLETCR,PRA_D.PRADWORDCR,PRA_D.PRADCR_RAW)

PRA_iw <- rbind(PRA_iw2_tomerge,PRA_iw3_tomerge)

PRA_repeat <- PRA_iw %>% group_by(bblid) %>% dplyr::summarize(n = n()) %>% 
  left_join(PRA_iw,.,by="bblid") %>% filter(n>1) %>% arrange(bblid,test_sessions_v.dotest)

PRA_repeat$test_sessions_v.dotest <- as.Date(PRA_repeat$test_sessions_v.dotest,format = "%Y-%m-%d")
PRA_repeat$days <- as.numeric(ymd(PRA_repeat$test_sessions_v.dotest))

# difference between timepoints in days
PRA_repeat$diffdays <- 0
for (i in 2:nrow(PRA_repeat)) {
  PRA_repeat$diffdays[i] <- PRA_repeat$days[i] - PRA_repeat$days[i-1]
}

temp <- PRA_repeat %>% dplyr::select("bblid","test_sessions_v.dotest")
t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")

PRA_repeat$combo <- paste(as.character(PRA_repeat$bblid),as.character(PRA_repeat$test_sessions_v.dotest),sep="_")
PRA_repeat$diffdays <- ifelse(PRA_repeat$combo %in% t.first$combo,0,PRA_repeat$diffdays)

# distribution of difference in days
{
  my_plot <- ggplot(repeats_only, aes(x = constant, y = diffdays)) + 
    ggdist::stat_halfeye(
      adjust = .5, 
      width = .6,
      justification = -.2, 
      .width = 0, 
      point_colour = NA,
      alpha = 0.8,
      fill  = "aquamarine3"
    ) + 
    geom_boxplot(
      width = .12, 
      outlier.color = NA, ## `outlier.shape = NA` works as well
      alpha = 0.5,
      color = "aquamarine3"
    ) +
    geom_point(
      size = .5, 
      alpha = 0.3,
      position = position_jitternudge(
        jitter.width = .1,
        jitter.height = 0,
        nudge.x = -.2,
        nudge.y = 0,
        seed = 1
      ),
      color = "aquamarine3"
    )  + 
    coord_cartesian(xlim = c(1.2, NA)) +
    theme_minimal() + labs(title = "Distribution of Extralong Timepoint Differences in Days",
                           x = "", y = "Difference (days)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/extralong_repeat_date_dist_220721.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
  
  # same thing as above but only people with difference < 100
  dat <- repeats_only %>% filter(diffdays <= 100)
  my_plot <- ggplot(dat, aes(x = constant, y = diffdays)) + 
    ggdist::stat_halfeye(
      adjust = .5, 
      width = .6,
      justification = -.2, 
      .width = 0, 
      point_colour = NA,
      alpha = 0.8,
      fill  = "aquamarine3"
    ) + 
    geom_boxplot(
      width = .12, 
      outlier.color = NA, ## `outlier.shape = NA` works as well
      alpha = 0.5,
      color = "aquamarine3"
    ) +
    geom_point(
      size = .5, 
      alpha = 0.3,
      position = position_jitternudge(
        jitter.width = .1,
        jitter.height = 0,
        nudge.x = -.2,
        nudge.y = 0,
        seed = 1
      ),
      color = "aquamarine3"
    )  + 
    coord_cartesian(xlim = c(1.2, NA)) +
    theme_minimal() + labs(title = "Distribution of Extralong Timepoint Differences in Days",
                           x = "", y = "Difference (days)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/extralong_repeat_date_dist_100_220721.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
  
  
  
  # age restrictions (18-35 yo)
  dat <- repeats_only %>% filter(diffdays <= 100,test_sessions_v.age %in% 18:35)  # n= 176 unique bblids
  my_plot <- ggplot(dat, aes(x = constant, y = diffdays)) + 
    ggdist::stat_halfeye(
      adjust = .5, 
      width = .6,
      justification = -.2, 
      .width = 0, 
      point_colour = NA,
      alpha = 0.8,
      fill  = "aquamarine3"
    ) + 
    geom_boxplot(
      width = .12, 
      outlier.color = NA, ## `outlier.shape = NA` works as well
      alpha = 0.5,
      color = "aquamarine3"
    ) +
    geom_point(
      size = .5, 
      alpha = 0.3,
      position = position_jitternudge(
        jitter.width = .1,
        jitter.height = 0,
        nudge.x = -.2,
        nudge.y = 0,
        seed = 1
      ),
      color = "aquamarine3"
    )  + 
    coord_cartesian(xlim = c(1.2, NA)) +
    theme_minimal() + labs(title = "Distribution of Extralong Timepoint Differences in Days, restricting age to 18-35 yo", 
                           caption = paste("n =",dim(dat)[1]),
                           x = "", y = "Difference (days)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/extralong_repeat_date_dist_age_220721.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
}





# stick with diffdays <= 60 for now ----
bbl_keep <- repeats_only %>% filter(diffdays <= 60, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 134, realized this includes people who may have t1's < 18  yo
bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 115

dat60 <- extralong_repeat %>% filter(diffdays <= 60, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
  arrange(bblid,test_sessions_v.dotest) %>% dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,mpraxis_genus:EDISC.test)   # n = 277, reordered to keep diffdays in demos

demos60 <- dat60 %>% dplyr::select(datasetid_platform:combo) 

dat60_t1 <- dat60 %>% filter(diffdays == 0) %>% arrange(bblid,test_sessions_v.dotest)  # n = 115
dat60_t2 <- dat60 %>% filter(diffdays > 0) %>% arrange(bblid,test_sessions_v.dotest)   # n = 161

demos_t1 <- dat60_t1 %>% dplyr::select(datasetid_platform:combo)
demos_t2 <- dat60_t2 %>% dplyr::select(datasetid_platform:combo)


adt60 <-  dat60 %>% dplyr::select(matches("adt")) %>% cbind(demos60,.) %>% filter(!is.na(adt_pc))  # 83 rows
aim60 <-  dat60 %>% dplyr::select(matches("aim")) %>% cbind(demos60,.) %>% filter(!is.na(aim_tot))  # 64 rows
cpf60 <-  dat60 %>% dplyr::select(matches("cpf")) %>% cbind(demos60,.) %>% filter(!is.na(cpf_cr))  # 200 rows




# first look at CPW
t1_CPW <- dat60_t1 %>% dplyr::select(matches("cpw")) %>% cbind(demos_t1,.) %>% filter(!is.na(cpw_cr)) # 85 rows
t2_CPW <- dat60_t2 %>% dplyr::select(matches("cpw")) %>% cbind(demos_t2,.) %>% filter(!is.na(cpw_cr)) # 82 rows

temp <- t2_CPW %>% dplyr::select("bblid","test_sessions_v.dotest")
t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")

t2_CPW$combo <- paste(as.character(t2_CPW$bblid),as.character(t2_CPW$test_sessions_v.dotest),sep="_")
t2_CPW2 <- t2_CPW %>% filter(combo %in% t.first$combo)


combo_CPW <- left_join(t1_CPW,t2_CPW2 %>% dplyr::select(bblid,diffdays,cpw_genus:cpw_w_rtcr), by="bblid") %>% filter(!is.na(cpw_cr.y))
combo_cpwa <- combo_CPW %>% filter(cpw_genus.x == "cpw_a",cpw_genus.y == "cpw_a") %>% rename(cpw_cr_t1 = cpw_cr.x,cpw_cr_t2 = cpw_cr.y)   # only 16 rows
combo_kcpwa <- combo_CPW %>% filter(cpw_genus.x == "kcpw_a",cpw_genus.y == "kcpw_a") %>% rename(cpw_cr_t1 = cpw_cr.x,cpw_cr_t2 = cpw_cr.y)   # only 4 rows

# scatters
# pdf("data/outputs/full_full/CPW_testretest_220725.pdf",height=9,width=12)
# pairs.panels(combo_cpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(combo_kcpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
# dev.off()





# PVRT STILL NEED TO FIX CODE FOR THIS
t1_CPW <- dat60_t1 %>% dplyr::select(matches("cpw")) %>% cbind(demos_t1,.) %>% filter(!is.na(cpw_cr)) # 85 rows
t2_CPW <- dat60_t2 %>% dplyr::select(matches("cpw")) %>% cbind(demos_t2,.) %>% filter(!is.na(cpw_cr)) # 82 rows

temp <- t2_CPW %>% dplyr::select("bblid","test_sessions_v.dotest")
t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")

t2_CPW$combo <- paste(as.character(t2_CPW$bblid),as.character(t2_CPW$test_sessions_v.dotest),sep="_")
t2_CPW2 <- t2_CPW %>% filter(combo %in% t.first$combo)


combo_CPW <- left_join(t1_CPW,t2_CPW2 %>% dplyr::select(bblid,diffdays,cpw_genus:cpw_w_rtcr), by="bblid") %>% filter(!is.na(cpw_cr.y))
combo_cpwa <- combo_CPW %>% filter(cpw_genus.x == "cpw_a",cpw_genus.y == "cpw_a") %>% rename(cpw_cr_t1 = cpw_cr.x,cpw_cr_t2 = cpw_cr.y)   # only 16 rows
combo_kcpwa <- combo_CPW %>% filter(cpw_genus.x == "kcpw_a",cpw_genus.y == "kcpw_a") %>% rename(cpw_cr_t1 = cpw_cr.x,cpw_cr_t2 = cpw_cr.y)   # only 4 rows

# scatters
# pdf("data/outputs/full_full/CPW_testretest_220725.pdf",height=9,width=12)
# pairs.panels(combo_cpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(combo_kcpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
# dev.off()







# CPF
t1_CPF <- dat60_t1 %>% dplyr::select(matches("cpf")) %>% cbind(demos_t1,.) %>% filter(!is.na(cpf_cr)) # 95 rows
t2_CPF <- dat60_t2 %>% dplyr::select(matches("cpf")) %>% cbind(demos_t2,.) %>% filter(!is.na(cpf_cr)) # 103 rows

temp <- t2_CPF %>% dplyr::select("bblid","test_sessions_v.dotest")
t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")

t2_CPF$combo <- paste(as.character(t2_CPF$bblid),as.character(t2_CPF$test_sessions_v.dotest),sep="_")
t2_CPF2 <- t2_CPF %>% filter(combo %in% t.first$combo)  # 69


combo_CPF <- left_join(t1_CPF,t2_CPF2 %>% dplyr::select(bblid,diffdays,cpf_genus:cpf_w_rtcr), by="bblid") %>% filter(!is.na(cpf_cr.y))
combo_cpfa <- combo_CPF %>% filter(cpf_genus.x == "cpf_a",cpf_genus.y == "cpf_a") %>% rename(cpf_cr_t1 = cpf_cr.x,cpf_cr_t2 = cpf_cr.y)   # only 16 rows
combo_cpfb <- combo_CPF %>% filter(cpf_genus.x == "cpf_b",cpf_genus.y == "cpf_b") %>% rename(cpf_cr_t1 = cpf_cr.x,cpf_cr_t2 = cpf_cr.y)   # only 4 rows

# scatters
# pdf("data/outputs/full_full/CPF_testretest_220725.pdf",height=9,width=12)
# pairs.panels(combo_cpfa %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(combo_cpfb %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
# dev.off()




# adjust diffdays to get n >=50 for each test ----
# need to get at least n = 50 (preferably more like 100, but don't want to widen 
# diffdays range too much) per test for good comparison
diffday_order <- repeats_only %>% arrange(diffdays) %>% filter(test_sessions_v.age %in% 18:35) %>% 
  dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,mpraxis_genus:EDISC.test) 
# dim(diffday_order %>% filter(diffdays <= 100)) # 258 rows



# creating new datasets with a loop ----

# * only first two timepoints ----
# PRA separate
tests <- c("adt","aim","cpf","cpt","cpw","ddisc","digsym","edisc","er40","gng","medf","plot","pmat","pvrt","rdisc","volt")  # no pra for now
test_sums <- c("adt_pc","aim_tot_","cpf_cr","cpt_acc","cpw_cr","ddisc_sum","dscor_","edisc_sum","er40_cr","gng_cr","medf_pc","plot_pc","pmat_pc","pvrt_pc","rdisc_sum","volt_cr")

demos <- extralong_repeat %>% dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo)

# take first two timepoints for comparison
{
  for (i in 1:length(tests)) {
    test <- tests[i]
    test_sum <- test_sums[i]
    
    dat <- extralong_repeat %>% dplyr::select(matches(test)) %>% cbind(demos,.) %>% 
      filter(test_sessions_v.age %in% 18:35)
    dat <- subset(dat, !is.na(dat[,ncol(dat)]))
    
    # digsym variables are saved under ds as well as digsym
    if (test == "digsym"){
      dat <- extralong_repeat %>% dplyr::select(matches("digsym|ds")) %>% cbind(demos,.) %>% 
        filter(test_sessions_v.age %in% 18:35)
    }
    
    # the last row has empty rows that aren't NA so it keeps too many rows
    if (test %in% c("digsym","gng","plot")){
      dat <- subset(dat, !is.na(dat[,ncol(dat)-1]))
    }
    
    # getting rid of invalid codes
    dat <- subset(dat, (dat[,grepl("_valid",colnames(dat))] != "N"))
    
    n <- unique(dat$bblid)
    
    temp <- dat %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t1 <- dat %>% filter(combo %in% t.first$combo)  
    dat_tn <- dat %>% filter(combo %notin% t.first$combo)
    
    temp <- dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t2 <- dat_tn %>% filter(combo %in% t.first$combo)  
    
    names(dat_t1)[21:ncol(dat_t1)] <- paste(names(dat_t1)[21:ncol(dat_t1)],"t1",sep="_")
    names(dat_t2)[c(9,19:ncol(dat_t2))] <- paste(names(dat_t2)[c(9,19:ncol(dat_t2))],"t2",sep="_")
    
    new_dat <- left_join(dat_t1,dat_t2[c(5,9,19:ncol(dat_t2))],by="bblid")
    new_dat <- subset(new_dat, !is.na(new_dat[,ncol(new_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                           dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                           new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
    
    assign(paste0("new_",test),new_dat)
  }
  
  # special step for EDISC -- need to see if I can get n > 46, ideally n > 50
  {
    test <- tests[8]
    test_sum <- test_sums[8]
    
    dat <- extralong_repeat %>% dplyr::select(matches(test)) %>% cbind(demos,.) %>% 
      filter(test_sessions_v.age %in% 18:35)
    dat <- subset(dat, !is.na(dat[,ncol(dat)]))
    
    n <- unique(dat$bblid)
    
    temp <- dat %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t1 <- dat %>% filter(combo %in% t.first$combo)  
    dat_tn <- dat %>% filter(combo %notin% t.first$combo)
    
    temp <- dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t2 <- dat_tn %>% filter(combo %in% t.first$combo) 
    dat_tn2 <- dat_tn %>% filter(combo %notin% t.first$combo)
    
    names(dat_t1)[21:ncol(dat_t1)] <- paste(names(dat_t1)[21:ncol(dat_t1)],"t1",sep="_")
    names(dat_t2)[c(9,19:ncol(dat_t2))] <- paste(names(dat_t2)[c(9,19:ncol(dat_t2))],"t2",sep="_")
    
    new_dat <- left_join(dat_t1,dat_t2[c(5,9,19:ncol(dat_t2))],by="bblid")
    new_dat <- subset(new_dat, !is.na(new_dat[,ncol(new_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                           dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                           new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
    
    temp <- dat_tn2 %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t3 <- dat_tn2 %>% filter(combo %in% t.first$combo)  
    
    
    names(dat_t3)[c(9,19:ncol(dat_t3))] <- paste(names(dat_t3)[c(9,19:ncol(dat_t3))],"t3",sep="_")
    
    new_dat1_3 <- left_join(dat_t1,dat_t3[c(5,9,19:ncol(dat_t3))],by="bblid")
    new_dat1_3 <- subset(new_dat1_3, !is.na(new_dat1_3[,ncol(new_dat1_3)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                       dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                       new_diffdays = dotest_t3 - dotest_t1) %>% arrange(new_diffdays)
    
    new_dat2_3 <- left_join(dat_t2,dat_t3[c(5,9,19:ncol(dat_t3))],by="bblid")
    new_dat2_3 <- subset(new_dat2_3, !is.na(new_dat2_3[,ncol(new_dat2_3)])) %>% mutate(dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                       dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                       new_diffdays = dotest_t3 - dotest_t2) %>% arrange(new_diffdays)
    names(new_dat1_3) <- names(new_dat)
    names(new_dat2_3) <- names(new_dat)
    
    new_dat <- rbind(new_dat,new_dat1_3,new_dat2_3)
    
    assign(paste0("new_",test),new_dat)
  }
  
  
  # printing the pairs.panels produced above
  newtexts <- paste0("new_",tests)
  newtests <- mget(newtexts)
  
  # make a table to accompany plots 
  # specifically, showing N, range of diffdays, median and mean diffdays
  test_row <- data.frame(matrix(NA,nrow = length(tests),ncol = 5))
  names(test_row) <- c("N","Min diffdays","Max diffdays","Mean diffdays","Median diffdays")
  rownames(test_row) <- toupper(tests)
  
  # adding QA method of getting rid of people w/ > 3 SD of a difference
  for (i in 1:length(newtests)) {
    test <- tests[i]
    test_dat <- newtests[[i]]
    temp <- test_dat[,grepl(test_sums[i],colnames(test_dat))]
    
    mod <- lm(temp[,2]~temp[,1],data=temp,na.action=na.exclude)
    # ggplot(new_adt, aes(x=adt_pc_t1, y=adt_pc_t2)) + geom_point() + geom_smooth(method='lm')
    sc <- scale(residuals(mod,na.action=na.exclude))
    test_dat <- cbind(test_dat,sc)
    test_dat$drop_3sd <- ifelse((test_dat$sc > 3| test_dat$sc < (-3)),1,0)
    
    assign(paste0("new_",test),test_dat)
    
    test_row[i,] <- c(dim(test_dat)[1],min(test_dat$new_diffdays),max(test_dat$new_diffdays),round(mean(test_dat$new_diffdays),2),median(test_dat$new_diffdays))
  }
  
  newtests <- mget(newtexts)
  
  # print both regular pairs.panels as well as drop_3sd implemented pairs.panels
  # pdf("data/outputs/full_full/all_testretest_noPRA_moreQC_220802.pdf",height=9,width=12)
  # for (i in 1:length(tests)) {
  #   pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  #   pairs.panels(newtests[[i]] %>% filter(drop_3sd != 1) %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  # }
  # dev.off()
  
  # print table
  # test_row %>% 
  #   kbl(caption = "More Info for Scatters", align = rep("c", 8)) %>%
  #   kable_classic(full_width = F, html_font = "Cambria") %>%
  #   save_kable(file = "data/outputs/full_full/full_full_info_table_220802.pdf", self_contained = T)
}


# keep only one row per person, limit data to max 365 diffdays, add more rows for tests that need it.
{
  for (i in 1:length(tests)) {
    test <- tests[i]
    test_sum <- test_sums[i]
    
    dat <- extralong_repeat %>% dplyr::select(matches(test)) %>% cbind(demos,.) %>% 
      filter(test_sessions_v.age %in% 18:35)
    dat <- subset(dat, !is.na(dat[,ncol(dat)]))
    
    # digsym variables are saved under ds as well as digsym
    if (test == "digsym"){
      dat <- extralong_repeat %>% dplyr::select(matches("digsym|ds")) %>% cbind(demos,.) %>% 
        filter(test_sessions_v.age %in% 18:35)
    }
    
    # the last row has empty rows that aren't NA so it keeps too many rows
    if (test %in% c("digsym","gng","plot")){
      dat <- subset(dat, !is.na(dat[,ncol(dat)-1]))
    }
    
    # getting rid of invalid codes
    if (test %in% c("ddisc","edisc","rdisc")){
      dat <- subset(dat, (dat[,grepl("DISC.valid_code",colnames(dat))] != "N"))
    } else {
      dat <- subset(dat, (dat[,grepl("_valid",colnames(dat))] != "N"))
    }
    
    n <- unique(dat$bblid)
    
    temp <- dat %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t1 <- dat %>% filter(combo %in% t.first$combo)  
    dat_tn <- dat %>% filter(combo %notin% t.first$combo)
    
    temp <- dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t2 <- dat_tn %>% filter(combo %in% t.first$combo)  
    
    names(dat_t1)[21:ncol(dat_t1)] <- paste(names(dat_t1)[21:ncol(dat_t1)],"t1",sep="_")
    names(dat_t2)[c(9,19:ncol(dat_t2))] <- paste(names(dat_t2)[c(9,19:ncol(dat_t2))],"t2",sep="_")
    
    new_dat <- left_join(dat_t1,dat_t2[c(5,9,19:ncol(dat_t2))],by="bblid")
    new_dat <- subset(new_dat, !is.na(new_dat[,ncol(new_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                           dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                           new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
    
    assign(paste0("new_",test),new_dat)
  }
  
  newtexts <- paste0("new_",tests)
  newtests <- mget(newtexts)
  
  # limiting to diffdays <= 365
  for (i in 1:length(tests)) {
    test <- tests[i]
    test_dat <- newtests[[i]]
    
    new_testdat <- test_dat %>% filter(new_diffdays <= 365)
    
    assign(paste0("new_",test,"_365"),new_testdat)
  }
  
  # tests that need more rows: DIGSYM (7), EDISC (8), RDISC (15)
  for (i in c(7,8,15)) {
    test <- tests[i]
    test_sum <- test_sums[i]
    
    dat <- extralong_repeat %>% dplyr::select(matches(test)) %>% cbind(demos,.) %>% 
      filter(test_sessions_v.age %in% 18:35)
    dat <- subset(dat, !is.na(dat[,ncol(dat)]))
    
    # digsym variables are saved under ds as well as digsym
    if (test == "digsym"){
      dat <- extralong_repeat %>% dplyr::select(matches("digsym|ds")) %>% cbind(demos,.) %>% 
        filter(test_sessions_v.age %in% 18:35)
    }
    
    # the last row has empty rows that aren't NA so it keeps too many rows
    if (test %in% c("digsym","gng","plot")){
      dat <- subset(dat, !is.na(dat[,ncol(dat)-1]))
    }
    
    # getting rid of invalid codes
    if (test %in% c("ddisc","edisc","rdisc")){
      dat <- subset(dat, (dat[,grepl("DISC.valid_code",colnames(dat))] != "N"))
    } else {
      dat <- subset(dat, (dat[,grepl("_valid",colnames(dat))] != "N"))
    }
    
    n <- unique(dat$bblid)
    
    temp <- dat %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t1 <- dat %>% filter(combo %in% t.first$combo)  
    dat_tn <- dat %>% filter(combo %notin% t.first$combo)
    
    temp <- dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t2 <- dat_tn %>% filter(combo %in% t.first$combo)  
    dat_tn2 <- dat_tn %>% filter(combo %notin% t.first$combo)
    
    names(dat_t1)[21:ncol(dat_t1)] <- paste(names(dat_t1)[21:ncol(dat_t1)],"t1",sep="_")
    names(dat_t2)[c(9,19:ncol(dat_t2))] <- paste(names(dat_t2)[c(9,19:ncol(dat_t2))],"t2",sep="_")
    
    new_dat <- left_join(dat_t1,dat_t2[c(5,9,19:ncol(dat_t2))],by="bblid")
    new_dat <- subset(new_dat, !is.na(new_dat[,ncol(new_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                           dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                           new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
    
    temp <- dat_tn2 %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t3 <- dat_tn2 %>% filter(combo %in% t.first$combo)  
    
    
    names(dat_t3)[c(9,19:ncol(dat_t3))] <- paste(names(dat_t3)[c(9,19:ncol(dat_t3))],"t3",sep="_")
    
    new_dat1_3 <- left_join(dat_t1,dat_t3[c(5,9,19:ncol(dat_t3))],by="bblid")
    new_dat1_3 <- subset(new_dat1_3, !is.na(new_dat1_3[,ncol(new_dat1_3)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                       dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                       new_diffdays = dotest_t3 - dotest_t1) %>% arrange(new_diffdays)
    
    new_dat2_3 <- left_join(dat_t2,dat_t3[c(5,9,19:ncol(dat_t3))],by="bblid")
    new_dat2_3 <- subset(new_dat2_3, !is.na(new_dat2_3[,ncol(new_dat2_3)])) %>% mutate(dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                       dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                       new_diffdays = dotest_t3 - dotest_t2) %>% arrange(new_diffdays)
    names(new_dat1_3) <- names(new_dat)
    names(new_dat2_3) <- names(new_dat)
    
    new_dat <- rbind(new_dat,new_dat1_3,new_dat2_3) %>% arrange(new_diffdays)
    
    assign(paste0("new_",test,"_123"),new_dat)
  }
  
  newtexts123 <- paste0("new_",tests[c(7,8,15)],"_123")
  newtests123 <- mget(newtexts123)
  
  # limiting to diffdays <= 365
  for (i in 1:length(newtexts123)) {
    test <- newtexts123[i]
    test_dat <- newtests123[[i]]
    
    new_testdat <- test_dat %>% filter(new_diffdays <= 365)
    
    assign(paste0(test,"_365"),new_testdat)
  }
  
  # printing the pairs.panels produced above
  newtexts <- c(paste0("new_",tests[1:6],"_365"),paste0(newtexts123[1:2],"_365"),paste0("new_",tests[9:14],"_365"),paste0(newtexts123[3],"_365"),paste0("new_",tests[16],"_365"))
  newtests <- mget(newtexts)
  
  # make a table to accompany plots 
  # specifically, showing N, range of diffdays, median and mean diffdays
  test_row <- data.frame(matrix(NA,nrow = length(tests),ncol = 6))
  names(test_row) <- c("N","N (after QC)","Min diffdays","Max diffdays","Mean diffdays","Median diffdays")
  rownames(test_row) <- toupper(tests)
  
  # adding QA method of getting rid of people w/ > 3 SD of a difference
  for (i in 1:length(newtests)) {
    test <- tests[i]
    test_dat <- newtests[[i]]
    temp <- test_dat[,grepl(test_sums[i],colnames(test_dat))]
    
    mod <- lm(temp[,2]~temp[,1],data=temp,na.action=na.exclude)
    # ggplot(temp, aes(x=temp[,1], y=temp[,2])) + geom_point() + geom_smooth(method='lm')
    sc <- scale(residuals(mod,na.action=na.exclude))
    test_dat <- cbind(test_dat,sc)
    test_dat$drop_3sd <- ifelse((test_dat$sc > 3| test_dat$sc < (-3)),1,0)
    sum_3sd <- sum(test_dat$drop_3sd)
    
    if (i %in% c(7,8,15)){
      assign(paste0("new_",test,"_123_365"),test_dat)
    } else {
      assign(paste0("new_",test,"_365"),test_dat)
    }
    
    test_row[i,] <- c(dim(test_dat)[1],dim(test_dat)[1]-sum_3sd,min(test_dat$new_diffdays),max(test_dat$new_diffdays),round(mean(test_dat$new_diffdays),1),median(test_dat$new_diffdays))
  }
  
  newtests <- mget(newtexts)
  
  
  
  
  # all of the above, for PRA
  dat <- PRA_repeat %>% filter(test_sessions_v.age %in% 18:35)
  
  # getting rid of invalid codes
  dat <- subset(dat, (dat[,grepl("D.valid_code",colnames(dat))] != "N" | is.na(dat[,grepl("D.valid_code",colnames(dat))])))
  
  #manually fix diffdays for bblid's who have two rows with same dotest (datasetid 94953,92134)
  dat[which(dat$test_sessions.datasetid == 94953),"combo"] <- "19590_2021-08-02_1"
  dat[which(dat$test_sessions.datasetid == 92134),"combo"] <- "19704_2021-04-28_1"
  
  n <- unique(dat$bblid)
  
  temp <- dat %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  dat_t1 <- dat %>% filter(combo %in% t.first$combo)  
  dat_tn <- dat %>% filter(combo %notin% t.first$combo)
  
  temp <- dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  dat_t2 <- dat_tn %>% filter(combo %in% t.first$combo)  
  dat_tn2 <- dat_tn %>% filter(combo %notin% t.first$combo)
  
  names(dat_t1)[11:ncol(dat_t1)] <- paste(names(dat_t1)[11:ncol(dat_t1)],"t1",sep="_")
  names(dat_t2)[c(6,8,11:ncol(dat_t2))] <- paste(names(dat_t2)[c(6,8,11:ncol(dat_t2))],"t2",sep="_")
  
  new_dat <- left_join(dat_t1,dat_t2[c(1,6,8,11:ncol(dat_t2))],by="bblid")
  new_dat <- subset(new_dat, !is.na(new_dat[,ncol(new_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                         dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                         new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays) # only 28 rows, need to add more timepoints
  
  temp <- dat_tn2 %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  dat_t3 <- dat_tn2 %>% filter(combo %in% t.first$combo)  
  
  names(dat_t3)[c(6,8,11:ncol(dat_t3))] <- paste(names(dat_t3)[c(6,8,11:ncol(dat_t3))],"t3",sep="_")
  
  new_dat1_3 <- left_join(dat_t1,dat_t3[c(1,6,8,11:ncol(dat_t3))],by="bblid")
  new_dat1_3 <- subset(new_dat1_3, !is.na(new_dat1_3[,ncol(new_dat1_3)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                     dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                     new_diffdays = dotest_t3 - dotest_t1) %>% arrange(new_diffdays)
  
  new_dat2_3 <- left_join(dat_t2,dat_t3[c(1,6,8,11:ncol(dat_t3))],by="bblid")
  new_dat2_3 <- subset(new_dat2_3, !is.na(new_dat2_3[,ncol(new_dat2_3)])) %>% mutate(dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                     dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                     new_diffdays = dotest_t3 - dotest_t2) %>% arrange(new_diffdays)
  names(new_dat1_3) <- names(new_dat)
  names(new_dat2_3) <- names(new_dat)
  
  new_pra <- rbind(new_dat,new_dat1_3,new_dat2_3) %>% arrange(new_diffdays)
  
  
  # print both regular pairs.panels as well as drop_3sd implemented pairs.panels
  pdf("data/outputs/full_full/all_testretest_someTP3_365_220810.pdf",height=9,width=12)
  for (i in 1:length(tests)) {
    pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
    # pairs.panels(newtests[[i]] %>% filter(drop_3sd != 1) %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  }
  pairs.panels(new_digsym_123_365 %>% dplyr::select(matches("dsmemcr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(new_dat %>% dplyr::select(matches("PRA_D.PRADWORDCR_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(new_pra %>% dplyr::select(matches("PRA_D.PRADWORDCR_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  # print table
  # test_row %>%
  #   kbl(caption = "More Info for Scatters", align = rep("c", 8)) %>%
  #   kable_classic(full_width = F, html_font = "Cambria") %>%
  #   save_kable(file = "data/outputs/full_full/full_full_info_table_someTP3_365_220804.pdf", self_contained = T)
  
}




# * first three timepoints ----
# expand what i did for EDISC to all tests
{
  for (i in 1:length(tests)) {
    test <- tests[i]
    test_sum <- test_sums[i]
    
    dat <- extralong_repeat %>% dplyr::select(matches(test)) %>% cbind(demos,.) %>% 
      filter(test_sessions_v.age %in% 18:35)
    dat <- subset(dat, !is.na(dat[,ncol(dat)]))
    
    # digsym variables are saved under ds as well as digsym
    if (test == "digsym"){
      dat <- extralong_repeat %>% dplyr::select(matches("digsym|ds")) %>% cbind(demos,.) %>% 
        filter(test_sessions_v.age %in% 18:35)
    }
    
    # the last row has empty rows that aren't NA so it keeps too many rows
    if (test %in% c("digsym","gng","plot")){
      dat <- subset(dat, !is.na(dat[,ncol(dat)-1]))
    }
    
    # getting rid of invalid codes
    dat <- subset(dat, (dat[,grepl("_valid",colnames(dat))] != "N"))
    
    n <- unique(dat$bblid)
    
    temp <- dat %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t1 <- dat %>% filter(combo %in% t.first$combo)  
    dat_tn <- dat %>% filter(combo %notin% t.first$combo)
    
    temp <- dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t2 <- dat_tn %>% filter(combo %in% t.first$combo)  
    dat_tn2 <- dat_tn %>% filter(combo %notin% t.first$combo)
    
    names(dat_t1)[21:ncol(dat_t1)] <- paste(names(dat_t1)[21:ncol(dat_t1)],"t1",sep="_")
    names(dat_t2)[c(9,19:ncol(dat_t2))] <- paste(names(dat_t2)[c(9,19:ncol(dat_t2))],"t2",sep="_")
    
    new_dat <- left_join(dat_t1,dat_t2[c(5,9,19:ncol(dat_t2))],by="bblid")
    new_dat <- subset(new_dat, !is.na(new_dat[,ncol(new_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                           dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                           new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
    
    temp <- dat_tn2 %>% dplyr::select("bblid","test_sessions_v.dotest")
    t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
    t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
    
    dat_t3 <- dat_tn2 %>% filter(combo %in% t.first$combo)  
    
    
    names(dat_t3)[c(9,19:ncol(dat_t3))] <- paste(names(dat_t3)[c(9,19:ncol(dat_t3))],"t3",sep="_")
    
    new_dat1_3 <- left_join(dat_t1,dat_t3[c(5,9,19:ncol(dat_t3))],by="bblid")
    new_dat1_3 <- subset(new_dat1_3, !is.na(new_dat1_3[,ncol(new_dat1_3)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                       dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                       new_diffdays = dotest_t3 - dotest_t1) %>% arrange(new_diffdays)
    
    new_dat2_3 <- left_join(dat_t2,dat_t3[c(5,9,19:ncol(dat_t3))],by="bblid")
    new_dat2_3 <- subset(new_dat2_3, !is.na(new_dat2_3[,ncol(new_dat2_3)])) %>% mutate(dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                       dotest_t3 = as.numeric(ymd(test_sessions_v.dotest_t3)),
                                                                                       new_diffdays = dotest_t3 - dotest_t2) %>% arrange(new_diffdays)
    names(new_dat1_3) <- names(new_dat)
    names(new_dat2_3) <- names(new_dat)
    
    new_dat <- rbind(new_dat,new_dat1_3,new_dat2_3) %>% arrange(new_diffdays)
    
    assign(paste0("new_",test,"2"),new_dat)
  }
  
  
  # printing the pairs.panels produced above
  newtexts <- paste0("new_",tests,"2")
  newtests <- mget(newtexts)
  
  # make a table to accompany plots 
  # specifically, showing N, range of diffdays, median and mean diffdays
  test_row <- data.frame(matrix(NA,nrow = length(tests),ncol = 5))
  names(test_row) <- c("N","Min diffdays","Max diffdays","Mean diffdays","Median diffdays")
  rownames(test_row) <- toupper(tests)
  
  # adding QA method of getting rid of people w/ > 3 SD of a difference
  for (i in 1:length(newtests)) {
    test <- tests[i]
    test_dat <- newtests[[i]]
    temp <- test_dat[,grepl(test_sums[i],colnames(test_dat))]
    
    mod <- lm(temp[,2]~temp[,1],data=temp,na.action=na.exclude)
    # ggplot(temp, aes(x=temp[,1], y=temp[,2])) + geom_point() + geom_smooth(method='lm')
    sc <- scale(residuals(mod,na.action=na.exclude))
    test_dat <- cbind(test_dat,sc)
    test_dat$drop_3sd <- ifelse((test_dat$sc > 3| test_dat$sc < (-3)),1,0)
    
    assign(paste0("new_",test,"2"),test_dat)
    
    test_row[i,] <- c(dim(test_dat)[1],min(test_dat$new_diffdays),max(test_dat$new_diffdays),round(mean(test_dat$new_diffdays),2),median(test_dat$new_diffdays))
  }
  
  newtests <- mget(newtexts)
  
  # print both regular pairs.panels as well as drop_3sd implemented pairs.panels
  # pdf("data/outputs/full_full/all_testretest_noPRA_moreQC_3tp_220802.pdf",height=9,width=12)
  # for (i in 1:length(tests)) {
  #   pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  #   pairs.panels(newtests[[i]] %>% filter(drop_3sd != 1) %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  # }
  # dev.off()
  
  # print table
  # test_row %>% 
  #   kbl(caption = "More Info for Scatters", align = rep("c", 8)) %>%
  #   kable_classic(full_width = F, html_font = "Cambria") %>%
  #   save_kable(file = "data/outputs/full_full/full_full_info_table_3tp_220802.pdf", self_contained = T)
}


# collection of plots with diffdays < 365 days or some other threshold
{
  for (i in 1:length(tests)) {
    test <- tests[i]
    test_dat <- newtests[[i]]
    
    new_testdat <- test_dat[,1:(ncol(test_dat)-2)] %>% filter(new_diffdays <= 365)
    
    assign(paste0("new_",test,"3"),new_testdat)
  }
  
  
  # printing the pairs.panels produced above
  newtexts <- paste0("new_",tests,"3")
  newtests <- mget(newtexts)
  
  # make a table to accompany plots 
  # specifically, showing N, range of diffdays, median and mean diffdays
  test_row <- data.frame(matrix(NA,nrow = length(tests),ncol = 5))
  names(test_row) <- c("N","Min diffdays","Max diffdays","Mean diffdays","Median diffdays")
  rownames(test_row) <- toupper(tests)
  
  # adding QA method of getting rid of people w/ > 3 SD of a difference
  for (i in 1:length(newtests)) {
    test <- tests[i]
    test_dat <- newtests[[i]]
    temp <- test_dat[,grepl(test_sums[i],colnames(test_dat))]
    
    mod <- lm(temp[,2]~temp[,1],data=temp,na.action=na.exclude)
    # ggplot(temp, aes(x=temp[,1], y=temp[,2])) + geom_point() + geom_smooth(method='lm')
    sc <- scale(residuals(mod,na.action=na.exclude))
    test_dat <- cbind(test_dat,sc)
    test_dat$drop_3sd <- ifelse((test_dat$sc > 3| test_dat$sc < (-3)),1,0)
    
    assign(paste0("new_",test,"3"),test_dat)
    
    test_row[i,] <- c(dim(test_dat)[1],min(test_dat$new_diffdays),max(test_dat$new_diffdays),round(mean(test_dat$new_diffdays),2),median(test_dat$new_diffdays))
  }
  
  newtests <- mget(newtexts)
  
  # print both regular pairs.panels as well as drop_3sd implemented pairs.panels
  # pdf("data/outputs/full_full/all_testretest_noPRA_moreQC_365_220802.pdf",height=9,width=12)
  # for (i in 1:length(tests)) {
  #   pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  #   pairs.panels(newtests[[i]] %>% filter(drop_3sd != 1) %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  # }
  # dev.off()
  
  # print table
  # test_row %>%
  #   kbl(caption = "More Info for Scatters", align = rep("c", 8)) %>%
  #   kable_classic(full_width = F, html_font = "Cambria") %>%
  #   save_kable(file = "data/outputs/full_full/full_full_info_table_365_220803.pdf", self_contained = T)
  
}




# checking memory tests
# starting with CPF: 4  groups
# 1) people who got the same test form (< 2 months apart from each other)
# 2) people who got the same test form (> 2 months apart from each other)
# 3) people who got different test form (< 2 months apart from each other)
# 4) people who got different test form (> 2 months apart from each other)

{
  cpf_dat <- extralong_repeat %>% dplyr::select(matches("cpf")) %>% cbind(demos,.) %>% 
    filter(test_sessions_v.age %in% 18:35)
  cpf_dat <- subset(cpf_dat, !is.na(cpf_dat[,ncol(cpf_dat)]))
  
  # getting rid of invalid codes
  cpf_dat <- subset(cpf_dat, (cpf_dat[,grepl("_valid",colnames(cpf_dat))] != "N"))
  
  cpf_dat_a <- cpf_dat %>% filter(cpf_genus == "cpf_a")
  cpf_dat_b <- cpf_dat %>% filter(cpf_genus == "cpf_b")
  
  # working with CPF A first
  temp <- cpf_dat_a %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpf_dat_a_t1 <- cpf_dat_a %>% filter(combo %in% t.first$combo)  # first time point of anyone who has taken CPF A
  cpf_dat_a_tn <- cpf_dat_a %>% filter(combo %notin% t.first$combo)
  
  temp <- cpf_dat_a_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),] 
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpf_dat_a_t2 <- cpf_dat_a_tn %>% filter(combo %in% t.first$combo)  # second time point of anyone who has taken CPF A
  
  names(cpf_dat_a_t1)[21:ncol(cpf_dat_a_t1)] <- paste(names(cpf_dat_a_t1)[21:ncol(cpf_dat_a_t1)],"t1",sep="_")
  names(cpf_dat_a_t2)[c(9,19:ncol(cpf_dat_a_t2))] <- paste(names(cpf_dat_a_t2)[c(9,19:ncol(cpf_dat_a_t2))],"t2",sep="_")
  
  new_cpf_dat_a <- left_join(cpf_dat_a_t1,cpf_dat_a_t2[c(5,9,19:ncol(cpf_dat_a_t2))],by="bblid")
  new_cpf_dat_a <- subset(new_cpf_dat_a, !is.na(new_cpf_dat_a[,ncol(new_cpf_dat_a)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                                 dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                                 new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
  
  cpf_dat_a_2mon <- new_cpf_dat_a %>% filter(new_diffdays <= 60)   # group 1.a, n = 25
  cpf_dat_a_2more <- new_cpf_dat_a %>% filter(new_diffdays > 60)   # group 2.a, n = 155
  
  # print scatters
  # pdf("data/outputs/full_full/CPF_A_60d_220808.pdf",height=9,width=12)
  # pairs.panels(cpf_dat_a_2mon %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(cpf_dat_a_2more %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  
  # same thing for CPF B
  temp <- cpf_dat_b %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpf_dat_b_t1 <- cpf_dat_b %>% filter(combo %in% t.first$combo)  # first time point of anyone who has taken CPF A
  cpf_dat_b_tn <- cpf_dat_b %>% filter(combo %notin% t.first$combo)
  
  temp <- cpf_dat_b_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),] 
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpf_dat_b_t2 <- cpf_dat_b_tn %>% filter(combo %in% t.first$combo)  # second time point of anyone who has taken CPF A
  
  names(cpf_dat_b_t1)[21:ncol(cpf_dat_b_t1)] <- paste(names(cpf_dat_b_t1)[21:ncol(cpf_dat_b_t1)],"t1",sep="_")
  names(cpf_dat_b_t2)[c(9,19:ncol(cpf_dat_b_t2))] <- paste(names(cpf_dat_b_t2)[c(9,19:ncol(cpf_dat_b_t2))],"t2",sep="_")
  
  new_cpf_dat_b <- left_join(cpf_dat_b_t1,cpf_dat_b_t2[c(5,9,19:ncol(cpf_dat_b_t2))],by="bblid")
  new_cpf_dat_b <- subset(new_cpf_dat_b, !is.na(new_cpf_dat_b[,ncol(new_cpf_dat_b)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                                 dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                                 new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
  
  cpf_dat_b_2mon <- new_cpf_dat_b %>% filter(new_diffdays <= 60)   # group 1.b, n = 11
  cpf_dat_b_2more <- new_cpf_dat_b %>% filter(new_diffdays > 60)   # group 2.b, n = 355
  
  # print scatters
  # pdf("data/outputs/full_full/CPF_B_60d_220808.pdf",height=9,width=12)
  # pairs.panels(cpf_dat_b_2mon %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(cpf_dat_b_2more %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()

  # combining forms
  cpf_dat_same_2mon <- rbind(cpf_dat_a_2mon,cpf_dat_b_2mon)
  cpf_dat_same_2more <- rbind(cpf_dat_a_2more,cpf_dat_b_2more)
  
  # print scatters
  # pdf("data/outputs/full_full/CPF_same_60d_220808.pdf",height=9,width=12)
  # pairs.panels(cpf_dat_same_2mon %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(cpf_dat_same_2more %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  # different forms
  temp <- cpf_dat %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpf_dat_t1 <- cpf_dat %>% filter(combo %in% t.first$combo)  # first time point of anyone who has taken CPF A
  cpf_dat_tn <- cpf_dat %>% filter(combo %notin% t.first$combo)
  
  temp <- cpf_dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),] 
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpf_dat_t2 <- cpf_dat_tn %>% filter(combo %in% t.first$combo)  # second time point of anyone who has taken CPF A
  
  names(cpf_dat_t1)[21:ncol(cpf_dat_t1)] <- paste(names(cpf_dat_t1)[21:ncol(cpf_dat_t1)],"t1",sep="_")
  names(cpf_dat_t2)[c(9,19:ncol(cpf_dat_t2))] <- paste(names(cpf_dat_t2)[c(9,19:ncol(cpf_dat_t2))],"t2",sep="_")
  
  new_cpf_dat <- left_join(cpf_dat_t1,cpf_dat_t2[c(5,9,19:ncol(cpf_dat_t2))],by="bblid")
  new_cpf_dat <- subset(new_cpf_dat, !is.na(new_cpf_dat[,ncol(new_cpf_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                         dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                         new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
  
  new_cpf_dat_dif <- new_cpf_dat %>% filter(cpf_genus_t1 != cpf_genus_t2)
  
  new_cpf_dat_dif_2mon <- new_cpf_dat_dif %>% filter(new_diffdays <= 60)   # group 3, n = 30
  new_cpf_dat_dif_2more <- new_cpf_dat_dif %>% filter(new_diffdays > 60)   # group 4, n = 236
  
  # print scatters
  # pdf("data/outputs/full_full/CPF_dif_60d_220808.pdf",height=9,width=12)
  # pairs.panels(new_cpf_dat_dif_2mon %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(new_cpf_dat_dif_2more %>% dplyr::select(matches("cpf_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
}


# same thing but with CPW
{
  cpw_dat <- extralong_repeat %>% dplyr::select(matches("cpw")) %>% cbind(demos,.) %>% 
    filter(test_sessions_v.age %in% 18:35)
  cpw_dat <- subset(cpw_dat, !is.na(cpw_dat[,ncol(cpw_dat)]))
  
  # getting rid of invalid codes
  cpw_dat <- subset(cpw_dat, (cpw_dat[,grepl("_valid",colnames(cpw_dat))] != "N"))
  
  cpw_dat_a <- cpw_dat %>% filter(cpw_genus == "cpw_a")
  cpw_dat_ka <- cpw_dat %>% filter(cpw_genus == "kcpw_a")
  
  # working with CPW A first
  temp <- cpw_dat_a %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpw_dat_a_t1 <- cpw_dat_a %>% filter(combo %in% t.first$combo)  # first time point of anyone who has taken CPW A
  cpw_dat_a_tn <- cpw_dat_a %>% filter(combo %notin% t.first$combo)
  
  temp <- cpw_dat_a_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),] 
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpw_dat_a_t2 <- cpw_dat_a_tn %>% filter(combo %in% t.first$combo)  # second time point of anyone who has taken CPW A
  
  names(cpw_dat_a_t1)[21:ncol(cpw_dat_a_t1)] <- paste(names(cpw_dat_a_t1)[21:ncol(cpw_dat_a_t1)],"t1",sep="_")
  names(cpw_dat_a_t2)[c(9,19:ncol(cpw_dat_a_t2))] <- paste(names(cpw_dat_a_t2)[c(9,19:ncol(cpw_dat_a_t2))],"t2",sep="_")
  
  new_cpw_dat_a <- left_join(cpw_dat_a_t1,cpw_dat_a_t2[c(5,9,19:ncol(cpw_dat_a_t2))],by="bblid")
  new_cpw_dat_a <- subset(new_cpw_dat_a, !is.na(new_cpw_dat_a[,ncol(new_cpw_dat_a)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                                 dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                                 new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
  
  cpw_dat_a_2mon <- new_cpw_dat_a %>% filter(new_diffdays <= 60)   # group 1.a, n = 21
  cpw_dat_a_2more <- new_cpw_dat_a %>% filter(new_diffdays > 60)   # group 2.a, n = 147
  
  # print scatters
  # pdf("data/outputs/full_full/CPW_A_60d_220808.pdf",height=9,width=12)
  # pairs.panels(cpw_dat_a_2mon %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(cpw_dat_a_2more %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  
  # same thing for CPW KA
  temp <- cpw_dat_ka %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpw_dat_ka_t1 <- cpw_dat_ka %>% filter(combo %in% t.first$combo)  # first time point of anyone who has taken CPW A
  cpw_dat_ka_tn <- cpw_dat_ka %>% filter(combo %notin% t.first$combo)
  
  temp <- cpw_dat_ka_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),] 
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpw_dat_ka_t2 <- cpw_dat_ka_tn %>% filter(combo %in% t.first$combo)  # second time point of anyone who has taken CPW A
  
  names(cpw_dat_ka_t1)[21:ncol(cpw_dat_ka_t1)] <- paste(names(cpw_dat_ka_t1)[21:ncol(cpw_dat_ka_t1)],"t1",sep="_")
  names(cpw_dat_ka_t2)[c(9,19:ncol(cpw_dat_ka_t2))] <- paste(names(cpw_dat_ka_t2)[c(9,19:ncol(cpw_dat_ka_t2))],"t2",sep="_")
  
  new_cpw_dat_ka <- left_join(cpw_dat_ka_t1,cpw_dat_ka_t2[c(5,9,19:ncol(cpw_dat_ka_t2))],by="bblid")
  new_cpw_dat_ka <- subset(new_cpw_dat_ka, !is.na(new_cpw_dat_ka[,ncol(new_cpw_dat_ka)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                                 dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                                 new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
  
  cpw_dat_ka_2mon <- new_cpw_dat_ka %>% filter(new_diffdays <= 60)   # group 1.b, n = 0
  cpw_dat_ka_2more <- new_cpw_dat_ka %>% filter(new_diffdays > 60)   # group 2.b, n = 398
  
  # print scatters
  # pdf("data/outputs/full_full/CPW_KA_60d_220808.pdf",height=9,width=12)
  # # pairs.panels(cpw_dat_ka_2mon %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE) # no obs.
  # pairs.panels(cpw_dat_ka_2more %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  # combining forms
  cpw_dat_same_2mon <- rbind(cpw_dat_a_2mon,cpw_dat_ka_2mon)
  cpw_dat_same_2more <- rbind(cpw_dat_a_2more,cpw_dat_ka_2more)
  
  # print scatters
  # pdf("data/outputs/full_full/CPW_same_60d_220808.pdf",height=9,width=12)
  # pairs.panels(cpw_dat_same_2mon %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(cpw_dat_same_2more %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  
  # different forms
  temp <- cpw_dat %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpw_dat_t1 <- cpw_dat %>% filter(combo %in% t.first$combo)  # first time point of anyone who has taken CPW A
  cpw_dat_tn <- cpw_dat %>% filter(combo %notin% t.first$combo)
  
  temp <- cpw_dat_tn %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),] 
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  cpw_dat_t2 <- cpw_dat_tn %>% filter(combo %in% t.first$combo)  # second time point of anyone who has taken CPW A
  
  names(cpw_dat_t1)[21:ncol(cpw_dat_t1)] <- paste(names(cpw_dat_t1)[21:ncol(cpw_dat_t1)],"t1",sep="_")
  names(cpw_dat_t2)[c(9,19:ncol(cpw_dat_t2))] <- paste(names(cpw_dat_t2)[c(9,19:ncol(cpw_dat_t2))],"t2",sep="_")
  
  new_cpw_dat <- left_join(cpw_dat_t1,cpw_dat_t2[c(5,9,19:ncol(cpw_dat_t2))],by="bblid")
  new_cpw_dat <- subset(new_cpw_dat, !is.na(new_cpw_dat[,ncol(new_cpw_dat)])) %>% mutate(dotest_t1 = as.numeric(ymd(test_sessions_v.dotest)),
                                                                                         dotest_t2 = as.numeric(ymd(test_sessions_v.dotest_t2)),
                                                                                         new_diffdays = dotest_t2 - dotest_t1) %>% arrange(new_diffdays)
  
  new_cpw_dat_dif <- new_cpw_dat %>% filter(cpw_genus_t1 != cpw_genus_t2)
  
  new_cpw_dat_dif_2mon <- new_cpw_dat_dif %>% filter(new_diffdays <= 60)   # group 3, n = 15
  new_cpw_dat_dif_2more <- new_cpw_dat_dif %>% filter(new_diffdays > 60)   # group 4, n = 91
  
  # print scatters
  # pdf("data/outputs/full_full/CPW_dif_60d_220808.pdf",height=9,width=12)
  # pairs.panels(new_cpw_dat_dif_2mon %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(new_cpw_dat_dif_2more %>% dplyr::select(matches("cpw_cr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
}







