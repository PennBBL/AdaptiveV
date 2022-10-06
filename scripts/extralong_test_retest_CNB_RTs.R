# Extralong test-retest (CNB RTs) as part of AdaptiveV

# Akira Di Sandro, 10.05.22


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



# load CSVs and create dataset ----

XL <- read.csv("data/inputs/cnb/cnb_merged_webcnp_surveys_smryscores_allbbl_longform.csv")  #30020 rows , 10/3/22
extralong <- XL %>% filter(test_sessions.bblid.clean>9999) %>% rename(bblid = test_sessions.bblid.clean)  # 16181 rows
# extralong <- XL %>% filter(test_sessions.bblid.clean>9999,test_sessions.siteid != "adaptive_v") %>% rename(bblid = test_sessions.bblid.clean)

# cpt_acc generation
extralong <- extralong %>% mutate(cpt_acc = cpt_ptp - cpt_pfp) %>% dplyr::select(datasetid_platform:bblid,test_sessions_v.age:cpt_pfp,cpt_acc,cpt_fprt:EDISC.test)

# disc scoring
{
  ddisc_qs <- extralong %>% dplyr::select((names(extralong)[grep('DDISC.q',colnames(extralong))]))
  extralong$ddisc_sum <- rowSums(ddisc_qs-1)
  ddisc_ttrs <- extralong %>% dplyr::select((names(extralong)[grep('DDISC.trr',colnames(extralong))]))
  extralong$ddisc_mcr <- rowMedians(as.matrix(ddisc_ttrs))
  
  rdisc_qs <- extralong %>% dplyr::select((names(extralong)[grep('RDISC.q',colnames(extralong))]))
  extralong$rdisc_sum <- rowSums(rdisc_qs-1)
  rdisc_ttrs <- extralong %>% dplyr::select((names(extralong)[grep('RDISC.trr',colnames(extralong))]))
  extralong$rdisc_mcr <- rowMedians(as.matrix(rdisc_ttrs))
  
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
}


# * creating extralong_repeat: keep only those with multiple timepoints ----
{
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
}


# No PRA necessary because we don't have RT for those anyways!


# Creating datasets for each test with RT ----
tests <- c("adt","aim","cpf","cpt","cpw","ddisc","digsym","edisc","er40","gng","medf","plot","pmat","pvrt","rdisc","volt")  # no pra for now
test_RTs <- c("adt_rtcr","aim_totrt","cpf_w_rtcr","cpt_tprt","cpw_w_rtcr","ddisc_mcr","dscorrt","edisc_mcr","er40_rtcr","gng_rtcr","medf_rtcr","plot_rtcr","pmat_rtcr","pvrt_rtcr","rdisc_mcr","volt_w_rtcr")

demos <- extralong_repeat %>% dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo)


# * 18-35yo, <= 365 ----
# 2 TP for most, except expand to 3 TP for those tests that don't have N >= 50 with 2 TP only
# keep only one row per person, limit data to max 365 diffdays, add more rows for tests that need it.
{
  for (i in 1:length(tests)) {
    test <- tests[i]
    test_RT <- test_RTs[i]
    
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
    # remove extreme outlier for CPW
    if (test == "cpw"){
      new_dat <- new_dat %>% filter(cpw_cr_t1 > 15)
    }
    
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
  
  # tests that need more rows: DIGSYM (i=7,n=39), EDISC (i=8,n=31), RDISC (i=15,n=30)
  for (i in c(7,8,15)) {
    test <- tests[i]
    test_RT <- test_RTs[i]
    
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
    temp <- test_dat[,grepl(test_RTs[i],colnames(test_dat))]
    
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
  
  
  # print both regular pairs.panels as well as drop_3sd implemented pairs.panels
  # pdf("data/outputs/full_full/all_testretest_someTP3_365_RT_221005.pdf",height=9,width=12)
  # for (i in 1:length(tests)) {
  #   pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_RTs[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  #   # pairs.panels(newtests[[i]] %>% filter(drop_3sd != 1) %>% dplyr::select(matches(test_RTs[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  # }
  # pairs.panels(new_digsym_123_365 %>% dplyr::select(matches("dsmemcr_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  # print table (with or without QC column)
  test_row %>%
    dplyr::select(!matches("QC")) %>%  # add this to ignoore the QC row
    kbl(caption = "More Info for Scatters", align = rep("c", 8)) %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    # save_kable(file = "data/outputs/full_full/full_full_info_table_someTP3_365_221005.pdf", self_contained = T)
    save_kable(file = "data/outputs/full_full/full_full_info_table_someTP3_365_RT_noQCcol_221005.pdf", self_contained = T)
  
}






