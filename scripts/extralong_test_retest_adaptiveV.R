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




# load CSVs and create dataset ----

extralong <- read.csv("data/inputs/cnb_merged_webcnp_surveys_allbblprjcts_longform.csv")  #29851 rows , 7/21/22
extralong <- extralong %>% filter(test_sessions.bblid.clean>9999) %>% rename(bblid = test_sessions.bblid.clean)  # 16093 rows

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
  
  pdf("data/outputs/extralong_repeat_date_dist_age_220721.pdf",height = 7,width = 10)
  my_plot
  dev.off()
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
pdf("data/outputs/full_full/CPW_testretest_220725.pdf",height=9,width=12)
pairs.panels(combo_cpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(combo_kcpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()





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
pdf("data/outputs/full_full/CPW_testretest_220725.pdf",height=9,width=12)
pairs.panels(combo_cpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(combo_kcpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()







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
pdf("data/outputs/full_full/CPF_testretest_220725.pdf",height=9,width=12)
pairs.panels(combo_cpfa %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(combo_cpfb %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()




# adjust diffdays to get n >=50 for each test ----
# need to get at least n = 50 (preferably more like 100, but don't want to widen 
# diffdays range too much) per test for good comparison
diffday_order <- repeats_only %>% arrange(diffdays) %>% filter(test_sessions_v.age %in% 18:35) %>% 
  dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,mpraxis_genus:EDISC.test) 
# dim(diffday_order %>% filter(diffdays <= 100)) # 258 rows



# creating new datasets with a loop

# pra from itemwise??
tests <- c("adt","aim","cpf","cpt","cpw","ddisc","digsym","edisc","er40","gng","medf","plot","pmat","pvrt","rdisc","volt")  # no pra for now
test_sums <- c("adt_pc","aim_tot_","cpf_cr","cpt_acc","cpw_cr","ddisc_sum","dscor_","edisc_sum","er40_cr","gng_cr","medf_pc","plot_pc","pmat_pc","pvrt_cr","rdisc_sum","volt_cr")

demos <- extralong_repeat %>% dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo)

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
  
  # getting rid of invalid codes, excluding disc tasks because they don't have any valid codes
  if (test %notin% c("ddisc","edisc","rdisc")){
    dat <- subset(dat, (dat[,grepl("_valid",colnames(dat))] %notin% c("F","N","V3")))
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


# printing the pairs.panels produced above
newtexts <- paste0("new_",tests)
newtests <- mget(newtexts)


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
}

newtests <- mget(newtexts)

pdf("data/outputs/full_full/all_testretest_noPRA_moreQC_220802.pdf",height=9,width=12)
for (i in 1:length(tests)) {
  pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(newtests[[i]] %>% filter(drop_3sd != 1) %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
}
dev.off()


















