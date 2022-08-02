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
  left_join(extralong,.,by="bblid") %>% filter(n>1,test_sessions.datasetid %notin% c(48505,42470)) %>% arrange(bblid,test_sessions_v.dotest)

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




# trying something new

# while loop? 

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
  
  if (test == "digsym"){
    dat <- extralong_repeat %>% dplyr::select(matches("digsym|ds")) %>% cbind(demos,.) %>% 
      filter(test_sessions_v.age %in% 18:35)
  }
  
  if (test %in% c("digsym","gng","plot")){
    dat <- subset(dat, !is.na(dat[,ncol(dat)-1]))
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
  
  # need to add QA method of getting rid of people w/ > 3 SD of a difference
  
  assign(paste0(test,"_new"),new_dat)
  pairs.panels(new_dat %>% dplyr::select(matches(test_sum)),lm=TRUE,scale=TRUE,ci=TRUE)
}


# printing the pairs.panels produced above
newtexts <- paste0(tests,"_new")
newtests <- mget(newtexts)

pdf("data/outputs/full_full/all_testretest_noPRA_220802.pdf",height=9,width=12)
for (i in 1:length(tests)) {
  pairs.panels(newtests[[i]] %>% dplyr::select(matches(test_sums[i])),lm=TRUE,scale=TRUE,ci=TRUE)
}
dev.off()






# old code (only comparing t1 to t2)
{
  bbl_keep100 <- repeats_only %>% filter(diffdays <= 100, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep100 <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep100) %>% pull(bblid) %>% unique()  # n = 146
  
  dat100 <- extralong_repeat %>% filter(diffdays <= 100, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep100) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,mpraxis_genus:edisc_mcr)
  
  demos100 <- dat100 %>% dplyr::select(datasetid_platform:combo)
  
  adt100 <-  dat100 %>% dplyr::select(matches("adt")) %>% cbind(demos100,.) %>% filter(!is.na(adt_pc))  # 116 rows
  aim100 <-  dat100 %>% dplyr::select(matches("aim")) %>% cbind(demos100,.) %>% filter(!is.na(aim_tot))  # 77 rows
  cpf100 <-  dat100 %>% dplyr::select(matches("cpf")) %>% cbind(demos100,.) %>% filter(!is.na(cpf_cr))  # 256 rows
  cpt100 <-  dat100 %>% dplyr::select(matches("cpt")) %>% cbind(demos100,.) %>% filter(!is.na(cpt_ptp))  # 129 rows
  cpw100 <-  dat100 %>% dplyr::select(matches("cpw")) %>% cbind(demos100,.) %>% filter(!is.na(cpw_cr))  # 221 rows
  ds100 <-  dat100 %>% dplyr::select(matches("ds|digsym")) %>% cbind(demos100,.) %>% filter(!is.na(dscor))  # 70 rows
  ddisc100 <-  dat100 %>% dplyr::select(matches("ddisc")) %>% cbind(demos100,.) %>% filter(!is.na(ddisc_sum))  # 79 rows
  edisc100 <-  dat100 %>% dplyr::select(matches("edisc")) %>% cbind(demos100,.) %>% filter(!is.na(edisc_sum))  # 64 rows
  er40100 <-  dat100 %>% dplyr::select(matches("er40")) %>% cbind(demos100,.) %>% filter(!is.na(er40_cr))  # 274 rows
  gng100 <-  dat100 %>% dplyr::select(matches("gng")) %>% cbind(demos100,.) %>% filter(!is.na(gng_cr))  # 64 rows
  medf100 <-  dat100 %>% dplyr::select(matches("medf")) %>% cbind(demos100,.) %>% filter(!is.na(medf_pc))  # 119 rows
  plot100 <-  dat100 %>% dplyr::select(matches("plot")) %>% cbind(demos100,.) %>% filter(!is.na(plot_pc))  # 100 rows
  pmat100 <-  dat100 %>% dplyr::select(matches("pmat")) %>% cbind(demos100,.) %>% filter(!is.na(pmat_pc))  # 103 rows
  # PRA is not in this data set
  # pra100 <-  dat100 %>% dplyr::select(matches("pra")) %>% cbind(demos100,.) %>% filter(!is.na(pra_pc))  # 64 rows
  pvrt100 <-  dat100 %>% dplyr::select(matches("pvrt")) %>% cbind(demos100,.) %>% filter(!is.na(pvrt_cr))  # 189 rows
  rdisc100 <-  dat100 %>% dplyr::select(matches("rdisc")) %>% cbind(demos100,.) %>% filter(!is.na(rdisc_sum))  # 75 rows
  volt100 <-  dat100 %>% dplyr::select(matches("volt")) %>% cbind(demos100,.) %>% filter(!is.na(volt_cr))  # 113 rows
  
  
  dat100_t1 <- dat100 %>% filter(diffdays == 0) %>% arrange(bblid,test_sessions_v.dotest)  # n = 146
  dat100_t2 <- dat100 %>% filter(diffdays > 0) %>% arrange(bblid,test_sessions_v.dotest)   # n = 208
  
  demos100_t1 <- dat100_t1 %>% dplyr::select(datasetid_platform:combo)
  demos100_t2 <- dat100_t2 %>% dplyr::select(datasetid_platform:combo)
  
  
  # ADT -- 54 rows
  # more data
  bbl_keep <- repeats_only %>% filter(diffdays <= 300, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 300, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(adt_pc)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,adt_genus:adt_rtcr)
  
  t1_ADT <- dat %>% filter(diffdays == 0) # n = 95
  t2_ADT <- dat %>% filter(diffdays > 0)  # n = 137
  
  # t1_ADT <- dat100_t1 %>% dplyr::select(matches("adt")) %>% cbind(demos100_t1,.) %>% filter(!is.na(adt_pc)) # 52 rows
  # t2_ADT <- dat100_t2 %>% dplyr::select(matches("adt")) %>% cbind(demos100_t2,.) %>% filter(!is.na(adt_pc)) # 64 rows
  
  temp <- t2_ADT %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_ADT$combo <- paste(as.character(t2_ADT$bblid),as.character(t2_ADT$test_sessions_v.dotest),sep="_")
  t2_ADT2 <- t2_ADT %>% filter(combo %in% t.first$combo)  # 110 rows
  
  
  combo_ADT <- left_join(t1_ADT,t2_ADT2 %>% dplyr::select(bblid,diffdays,adt_genus:adt_rtcr), by="bblid") %>% filter(!is.na(adt_pc.y)) %>%  # 54 rows
    rename(adt_pc_t1 = adt_pc.x,adt_pc_t2 = adt_pc.y)
  combo_adt60a <- combo_ADT %>% filter(adt_genus.x == "adt60_a",adt_genus.y == "adt60_a")   # 2 rows
  combo_adt36a <- combo_ADT %>% filter(adt_genus.x == "adt36_a",adt_genus.y == "adt36_a")   # 33 rows
  combo_adt36b <- combo_ADT %>% filter(adt_genus.x == "adt36_b",adt_genus.y == "adt36_b")   # 1 rows
  
  # scatters
  pdf("data/outputs/full_full/ADT_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_ADT %>% dplyr::select(matches("adt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_adt36a %>% dplyr::select(matches("adt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_adt36b %>% dplyr::select(matches("adt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_adt6a0 %>% dplyr::select(matches("adt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # AIM -- 61 rows
  bbl_keep <- repeats_only %>% filter(diffdays <= 200, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 200, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(aim_tot)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,aim_genus:aim_mcrrt)
  
  t1_AIM <- dat %>% filter(diffdays == 0) # n = 66
  t2_AIM <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_AIM <- dat100_t1 %>% dplyr::select(matches("aim")) %>% cbind(demos100_t1,.) %>% filter(!is.na(aim_tot)) # 38 rows
  # t2_AIM <- dat100_t2 %>% dplyr::select(matches("aim")) %>% cbind(demos100_t2,.) %>% filter(!is.na(aim_tot)) # 39 rows
  
  temp <- t2_AIM %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_AIM$combo <- paste(as.character(t2_AIM$bblid),as.character(t2_AIM$test_sessions_v.dotest),sep="_")
  t2_AIM2 <- t2_AIM %>% filter(combo %in% t.first$combo)  # 79 rows
  
  
  combo_AIM <- left_join(t1_AIM,t2_AIM2 %>% dplyr::select(bblid,diffdays,aim_genus:aim_mcrrt), by="bblid") %>% filter(!is.na(aim_tot.y)) %>%  # 61 rows
    rename(aim_tot_t1 = aim_tot.x,aim_tot_t2 = aim_tot.y)
  
  # scatters
  pdf("data/outputs/full_full/AIM_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_AIM %>% dplyr::select(matches("aim_tot_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # CPF -- 74, good
  t1_CPF <- dat100_t1 %>% dplyr::select(matches("cpf")) %>% cbind(demos100_t1,.) %>% filter(!is.na(cpf_cr)) # 117 rows
  t2_CPF <- dat100_t2 %>% dplyr::select(matches("cpf")) %>% cbind(demos100_t2,.) %>% filter(!is.na(cpf_cr)) # 139 rows
  
  temp <- t2_CPF %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_CPF$combo <- paste(as.character(t2_CPF$bblid),as.character(t2_CPF$test_sessions_v.dotest),sep="_")
  t2_CPF2 <- t2_CPF %>% filter(combo %in% t.first$combo)  # 98 rows
  
  
  combo_CPF <- left_join(t1_CPF,t2_CPF2 %>% dplyr::select(bblid,diffdays,cpf_genus:cpf_w_rtcr), by="bblid") %>% filter(!is.na(cpf_cr.y)) %>%  # 74 rows
    rename(cpf_cr_t1 = cpf_cr.x,cpf_cr_t2 = cpf_cr.y)
  combo_cpfa <- combo_CPF %>% filter(cpf_genus.x == "cpf_a",cpf_genus.y == "cpf_a")   # 28 rows
  combo_cpfb <- combo_CPF %>% filter(cpf_genus.x == "cpf_b",cpf_genus.y == "cpf_b")   # 9 rows
  
  # scatters
  pdf("data/outputs/full_full/CPF_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_CPF %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_cpfa %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_cpfb %>% dplyr::select(matches("cpf_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # CPT -- 29, need to expand
  bbl_keep <- repeats_only %>% filter(diffdays <= 250, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 250, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(cpt_ptp)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,cpt_genus:cpt_n_tprt)
  
  t1_CPT <- dat %>% filter(diffdays == 0) %>% mutate(cpt_acc = cpt_ptp - cpt_pfp) # n = 73
  t2_CPT <- dat %>% filter(diffdays > 0) %>% mutate(cpt_acc = cpt_ptp - cpt_pfp)  # n = 120
  
  # t1_CPT <- dat100_t1 %>% dplyr::select(matches("cpt")) %>% cbind(demos100_t1,.) %>% filter(!is.na(cpt_ptp)) %>% 
  #   mutate(cpt_acc = cpt_ptp - cpt_pfp) # 54 rows
  # t2_CPT <- dat100_t2 %>% dplyr::select(matches("cpt")) %>% cbind(demos100_t2,.) %>% filter(!is.na(cpt_ptp)) %>% 
  #   mutate(cpt_acc = cpt_ptp - cpt_pfp) # 75 rows
  
  temp <- t2_CPT %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_CPT$combo <- paste(as.character(t2_CPT$bblid),as.character(t2_CPT$test_sessions_v.dotest),sep="_")
  t2_CPT2 <- t2_CPT %>% filter(combo %in% t.first$combo)  # 99 rows
  
  
  combo_CPT <- left_join(t1_CPT,t2_CPT2 %>% dplyr::select(bblid,diffdays,cpt_genus:cpt_acc), by="bblid") %>% filter(!is.na(cpt_acc.y)) %>%  # 62 rows
    rename(cpt_acc_t1 = cpt_acc.x,cpt_acc_t2 = cpt_acc.y)
  
  # scatters
  pdf("data/outputs/full_full/CPT_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_CPT %>% dplyr::select(matches("cpt_acc")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # CPW
  t1_CPW <- dat100_t1 %>% dplyr::select(matches("cpw")) %>% cbind(demos100_t1,.) %>% filter(!is.na(cpw_cr)) # 108 rows
  t2_CPW <- dat100_t2 %>% dplyr::select(matches("cpw")) %>% cbind(demos100_t2,.) %>% filter(!is.na(cpw_cr)) # 113 rows
  
  temp <- t2_CPW %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_CPW$combo <- paste(as.character(t2_CPW$bblid),as.character(t2_CPW$test_sessions_v.dotest),sep="_")
  t2_CPW2 <- t2_CPW %>% filter(combo %in% t.first$combo)  # 89 rows
  
  
  combo_CPW <- left_join(t1_CPW,t2_CPW2 %>% dplyr::select(bblid,diffdays,cpw_genus:cpw_w_rtcr), by="bblid") %>% filter(!is.na(cpw_cr.y)) %>%  # 58 rows
    rename(cpw_cr_t1 = cpw_cr.x,cpw_cr_t2 = cpw_cr.y) %>% filter(cpw_cr_t1 > 7)  # removing extreme outlier
  combo_cpwa <- combo_CPW %>% filter(cpw_genus.x == "cpw_a",cpw_genus.y == "cpw_a")  # 22 rows
  combo_kcpwa <- combo_CPW %>% filter(cpw_genus.x == "kcpw_a",cpw_genus.y == "kcpw_a")  # 7 rows
  
  # scatters
  pdf("data/outputs/full_full/CPW_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_CPW %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_cpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_kcpwa %>% dplyr::select(matches("cpw_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # DDISC -- 3 rows, need to expand
  bbl_keep <- repeats_only %>% filter(diffdays <= 2000, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 2000, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(ddisc_sum)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,KDDISC.valid_code:ddisc_mcr)
  
  t1_DDISC <- dat %>% filter(diffdays == 0) # n = 66
  t2_DDISC <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_DDISC <- dat100_t1 %>% dplyr::select(matches("ddisc")) %>% cbind(demos100_t1,.) %>% filter(!is.na(ddisc_sum)) # 8 rows
  # t2_DDISC <- dat100_t2 %>% dplyr::select(matches("ddisc")) %>% cbind(demos100_t2,.) %>% filter(!is.na(ddisc_sum)) # 71 rows
  
  temp <- t2_DDISC %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_DDISC$combo <- paste(as.character(t2_DDISC$bblid),as.character(t2_DDISC$test_sessions_v.dotest),sep="_")
  t2_DDISC2 <- t2_DDISC %>% filter(combo %in% t.first$combo)  # 36 rows
  
  
  combo_DDISC <- left_join(t1_DDISC,t2_DDISC2 %>% dplyr::select(bblid,diffdays,KDDISC.valid_code:ddisc_mcr), by="bblid") %>% filter(!is.na(ddisc_sum.y)) %>%  # 3 rows
    rename(ddisc_sum_t1 = ddisc_sum.x,ddisc_sum_t2 = ddisc_sum.y)
  
  # scatters
  pdf("data/outputs/full_full/DDISC_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_DDISC %>% dplyr::select(matches("ddisc_sum")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # DIGSYM -- 4 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 4000, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 4000, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(dscor)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,digsym_genus:digsym_valid)
  
  t1_DIGSYM <- dat %>% filter(diffdays == 0) # n = 66
  t2_DIGSYM <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_DIGSYM <- dat100_t1 %>% dplyr::select(matches("digsym|^ds")) %>% cbind(demos100_t1,.) %>% filter(!is.na(dscor)) # 5 rows
  # t2_DIGSYM <- dat100_t2 %>% dplyr::select(matches("digsym|^ds")) %>% cbind(demos100_t2,.) %>% filter(!is.na(dscor)) # 65 rows
  
  temp <- t2_DIGSYM %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_DIGSYM$combo <- paste(as.character(t2_DIGSYM$bblid),as.character(t2_DIGSYM$test_sessions_v.dotest),sep="_")
  t2_DIGSYM2 <- t2_DIGSYM %>% filter(combo %in% t.first$combo)  # 89 rows
  
  
  combo_DIGSYM <- left_join(t1_DIGSYM,t2_DIGSYM2 %>% dplyr::select(bblid,diffdays,digsym_genus:digsym_valid), by="bblid") %>% filter(!is.na(dscor.y)) %>%  # 3 rows
    rename(dscor_t1 = dscor.x,dscor_t2 = dscor.y)
  
  # scatters
  pdf("data/outputs/full_full/DIGSYM_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_DIGSYM %>% dplyr::select(matches("dscor_t")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # EDISC -- 7 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 4000, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 4000, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(edisc_sum)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,EDISC.valid_code:edisc_mcr)
  
  t1_EDISC <- dat %>% filter(diffdays == 0) # n = 66
  t2_EDISC <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_EDISC <- dat100_t1 %>% dplyr::select(matches("edisc")) %>% cbind(demos100_t1,.) %>% filter(!is.na(edisc_sum)) # 4 rows
  # t2_EDISC <- dat100_t2 %>% dplyr::select(matches("edisc")) %>% cbind(demos100_t2,.) %>% filter(!is.na(edisc_sum)) # 60 rows
  
  temp <- t2_EDISC %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_EDISC$combo <- paste(as.character(t2_EDISC$bblid),as.character(t2_EDISC$test_sessions_v.dotest),sep="_")
  t2_EDISC2 <- t2_EDISC %>% filter(combo %in% t.first$combo)  # 25 rows
  
  
  combo_EDISC <- left_join(t1_EDISC,t2_EDISC2 %>% dplyr::select(bblid,diffdays,EDISC.valid_code:edisc_mcr), by="bblid") %>% filter(!is.na(edisc_sum.y)) %>%  # 2 rows
    rename(edisc_sum_t1 = edisc_sum.x,edisc_sum_t2 = edisc_sum.y)
  
  # scatters
  pdf("data/outputs/full_full/EDISC_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_EDISC %>% dplyr::select(matches("edisc_sum")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # ER40
  t1_ER40 <- dat100_t1 %>% dplyr::select(matches("er40")) %>% cbind(demos100_t1,.) %>% filter(!is.na(er40_cr)) # 132 rows
  t2_ER40 <- dat100_t2 %>% dplyr::select(matches("er40")) %>% cbind(demos100_t2,.) %>% filter(!is.na(er40_cr)) # 142 rows
  
  temp <- t2_ER40 %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_ER40$combo <- paste(as.character(t2_ER40$bblid),as.character(t2_ER40$test_sessions_v.dotest),sep="_")
  t2_ER402 <- t2_ER40 %>% filter(combo %in% t.first$combo)  # 120 rows
  
  
  combo_ER40 <- left_join(t1_ER40,t2_ER402 %>% dplyr::select(bblid,diffdays,er40_genus:er40_rtcr), by="bblid") %>% filter(!is.na(er40_cr.y)) %>%  # 109 rows
    rename(er40_cr_t1 = er40_cr.x,er40_cr_t2 = er40_cr.y)
  combo_er40a <- combo_ER40 %>% filter(er40_genus.x == "ER40_A",er40_genus.y == "ER40_A")  # 43 rows
  combo_er40c <- combo_ER40 %>% filter(er40_genus.x == "ER40_C",er40_genus.y == "ER40_C")  # 1 rows
  combo_er40d <- combo_ER40 %>% filter(er40_genus.x == "ER40_D",er40_genus.y == "ER40_D")  # 19 rows
  
  # scatters
  pdf("data/outputs/full_full/ER40_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_ER40 %>% dplyr::select(matches("er40_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_er40a %>% dplyr::select(matches("er40_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_er40c %>% dplyr::select(matches("er40_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_er40d %>% dplyr::select(matches("er40_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # GNG -- 28 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 4000, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(test_sessions_v.age %in% 18:35) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(gng_cr)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,gng_genus:gng_valid)
  
  t1_GNG <- dat %>% filter(diffdays == 0) # n = 66
  t2_GNG <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_GNG <- dat100_t1 %>% dplyr::select(matches("gng")) %>% cbind(demos100_t1,.) %>% filter(!is.na(gng_cr)) # 20 rows
  # t2_GNG <- dat100_t2 %>% dplyr::select(matches("gng")) %>% cbind(demos100_t2,.) %>% filter(!is.na(gng_cr)) # 44 rows
  
  temp <- t2_GNG %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_GNG$combo <- paste(as.character(t2_GNG$bblid),as.character(t2_GNG$test_sessions_v.dotest),sep="_")
  t2_GNG2 <- t2_GNG %>% filter(combo %in% t.first$combo)  # 37 rows
  
  
  combo_GNG <- left_join(t1_GNG,t2_GNG2 %>% dplyr::select(bblid,diffdays,gng_genus:gng_valid), by="bblid") %>% filter(!is.na(gng_cr.y)) %>%  # 5 rows
    rename(gng_cr_t1 = gng_cr.x,gng_cr_t2 = gng_cr.y) %>% filter(gng_cr_t2 > 90) # remove extreme outlier of gng_cr_t2 < 90
  
  # scatters
  pdf("data/outputs/full_full/GNG_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_GNG %>% dplyr::select(matches("gng_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # MEDF -- 56 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 300, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 300, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(medf_pc)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,medf_genus:medf_rtcr)
  
  t1_MEDF <- dat %>% filter(diffdays == 0) # n = 66
  t2_MEDF <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_MEDF <- dat100_t1 %>% dplyr::select(matches("medf")) %>% cbind(demos100_t1,.) %>% filter(!is.na(medf_pc)) # 53 rows
  # t2_MEDF <- dat100_t2 %>% dplyr::select(matches("medf")) %>% cbind(demos100_t2,.) %>% filter(!is.na(medf_pc)) # 66 rows
  
  temp <- t2_MEDF %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_MEDF$combo <- paste(as.character(t2_MEDF$bblid),as.character(t2_MEDF$test_sessions_v.dotest),sep="_")
  t2_MEDF2 <- t2_MEDF %>% filter(combo %in% t.first$combo)  # 56 rows
  
  
  combo_MEDF <- left_join(t1_MEDF,t2_MEDF2 %>% dplyr::select(bblid,diffdays,medf_genus:medf_rtcr), by="bblid") %>% filter(!is.na(medf_pc.y)) %>%  # 21 rows
    rename(medf_pc_t1 = medf_pc.x,medf_pc_t2 = medf_pc.y)
  combo_medf36a <- combo_MEDF %>% filter(medf_genus.x == "medf36_a",medf_genus.y == "medf36_a")   # 34 rows
  combo_medf60a <- combo_MEDF %>% filter(medf_genus.x == "medf60_a",medf_genus.y == "medf60_a")   # 2 rows
  
  # scatters
  pdf("data/outputs/full_full/MEDF_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_MEDF %>% dplyr::select(matches("medf_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_medf36a %>% dplyr::select(matches("medf_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_medf60a %>% dplyr::select(matches("medf_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # PLOT -- 14 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 550, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 550, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(plot_pc)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,plot_genus:plot_rtcr)
  
  t1_PLOT <- dat %>% filter(diffdays == 0) # n = 66
  t2_PLOT <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_PLOT <- dat100_t1 %>% dplyr::select(matches("plot")) %>% cbind(demos100_t1,.) %>% filter(!is.na(plot_pc)) # 43 rows
  # t2_PLOT <- dat100_t2 %>% dplyr::select(matches("plot")) %>% cbind(demos100_t2,.) %>% filter(!is.na(plot_pc)) # 57 rows
  
  temp <- t2_PLOT %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_PLOT$combo <- paste(as.character(t2_PLOT$bblid),as.character(t2_PLOT$test_sessions_v.dotest),sep="_")
  t2_PLOT2 <- t2_PLOT %>% filter(combo %in% t.first$combo)  # 45 rows
  
  
  combo_PLOT <- left_join(t1_PLOT,t2_PLOT2 %>% dplyr::select(bblid,diffdays,plot_genus:plot_rtcr), by="bblid") %>% filter(!is.na(plot_pc.y)) %>%  # 14 rows
    rename(plot_pc_t1 = plot_pc.x,plot_pc_t2 = plot_pc.y)
  combo_splot12 <- combo_PLOT %>% filter(plot_genus.x == "splot12",plot_genus.y == "splot12")   # 3 rows
  combo_vsplot15 <- combo_PLOT %>% filter(plot_genus.x == "vsplot15",plot_genus.y == "vsplot15")   # 21 rows
  combo_vsplot24 <- combo_PLOT %>% filter(plot_genus.x == "vsplot24",plot_genus.y == "vsplot24")   # 2 rows
  
  # scatters
  pdf("data/outputs/full_full/PLOT_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_PLOT %>% dplyr::select(matches("plot_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_splot12 %>% dplyr::select(matches("plot_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_vsplot15 %>% dplyr::select(matches("plot_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_vsplot24 %>% dplyr::select(matches("plot_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # PMAT -- 53 rows
  bbl_keep <- repeats_only %>% filter(diffdays <= 400, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 400, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(pmat_pc)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,pmat_genus:pmat_rtcr)
  
  t1_PMAT <- dat %>% filter(diffdays == 0) # n = 66
  t2_PMAT <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_PMAT <- dat100_t1 %>% dplyr::select(matches("pmat")) %>% cbind(demos100_t1,.) %>% filter(!is.na(pmat_pc)) # 49 rows
  # t2_PMAT <- dat100_t2 %>% dplyr::select(matches("pmat")) %>% cbind(demos100_t2,.) %>% filter(!is.na(pmat_pc)) # 54 rows
  
  temp <- t2_PMAT %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_PMAT$combo <- paste(as.character(t2_PMAT$bblid),as.character(t2_PMAT$test_sessions_v.dotest),sep="_")
  t2_PMAT2 <- t2_PMAT %>% filter(combo %in% t.first$combo)  # 44 rows
  
  
  combo_PMAT <- left_join(t1_PMAT,t2_PMAT2 %>% dplyr::select(bblid,diffdays,pmat_genus:pmat_rtcr), by="bblid") %>% filter(!is.na(pmat_pc.y)) %>%  # 17 rows
    rename(pmat_pc_t1 = pmat_pc.x,pmat_pc_t2 = pmat_pc.y)
  combo_pmat18_a <- combo_PMAT %>% filter(pmat_genus.x == "pmat18_a",pmat_genus.y == "pmat18_a")   # 1 rows
  combo_pmat24_a <- combo_PMAT %>% filter(pmat_genus.x == "pmat24_a",pmat_genus.y == "pmat24_a")   # 23 rows
  combo_pmat24_b <- combo_PMAT %>% filter(pmat_genus.x == "pmat24_b",pmat_genus.y == "pmat24_b")   # 6 rows
  
  # scatters
  pdf("data/outputs/full_full/PMAT_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_PMAT %>% dplyr::select(matches("pmat_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_pmat18_a %>% dplyr::select(matches("pmat_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_pmat24_a %>% dplyr::select(matches("pmat_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_pmat24_b %>% dplyr::select(matches("pmat_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # PVRT -- 56 rows (diffdays <= 150; 78 rows diffdays <= 175), need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 175, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 175, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(pvrt_pc)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,pvrt_genus:pvrt_pc)
  
  t1_PVRT <- dat %>% filter(diffdays == 0) # n = 66
  t2_PVRT <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_PVRT <- dat100_t1 %>% dplyr::select(matches("pvrt")) %>% cbind(demos100_t1,.) %>% filter(!is.na(pvrt_pc)) # 90 rows
  # t2_PVRT <- dat100_t2 %>% dplyr::select(matches("pvrt")) %>% cbind(demos100_t2,.) %>% filter(!is.na(pvrt_pc)) # 99 rows
  
  temp <- t2_PVRT %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_PVRT$combo <- paste(as.character(t2_PVRT$bblid),as.character(t2_PVRT$test_sessions_v.dotest),sep="_")
  t2_PVRT2 <- t2_PVRT %>% filter(combo %in% t.first$combo)  # 86 rows
  
  
  combo_PVRT <- left_join(t1_PVRT,t2_PVRT2 %>% dplyr::select(bblid,diffdays,pvrt_genus:pvrt_pc), by="bblid") %>% filter(!is.na(pvrt_pc.y)) %>%  # 41 rows
    rename(pvrt_pc_t1 = pvrt_pc.x,pvrt_pc_t2 = pvrt_pc.y) %>% filter(abs(pvrt_pc_t1-pvrt_pc_t2) < 60) # remove outliers whose pc differ by > 60%
  combo_KSPVRT_A <- combo_PVRT %>% filter(pvrt_genus.x == "KSPVRT_A",pvrt_genus.y == "KSPVRT_A")   # 0 rows
  combo_KSPVRT_B <- combo_PVRT %>% filter(pvrt_genus.x == "KSPVRT_B",pvrt_genus.y == "KSPVRT_B")   # 0 rows
  combo_KSPVRT_D <- combo_PVRT %>% filter(pvrt_genus.x == "KSPVRT_D",pvrt_genus.y == "KSPVRT_D")   # 9 rows
  combo_SPVRT_A <- combo_PVRT %>% filter(pvrt_genus.x == "SPVRT_A",pvrt_genus.y == "SPVRT_A")   # 34 rows
  
  # scatters
  pdf("data/outputs/full_full/PVRT_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_PVRT %>% dplyr::select(matches("pvrt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_KSPVRT_A %>% dplyr::select(matches("pvrt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(combo_KSPVRT_B %>% dplyr::select(matches("pvrt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_KSPVRT_D %>% dplyr::select(matches("pvrt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(combo_SPVRT_A %>% dplyr::select(matches("pvrt_pc")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # RDISC -- 2 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 1000, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(test_sessions_v.age %in% 18:35) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(rdisc_sum)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,KRDISC.valid_code:rdisc_mcr)
  
  t1_RDISC <- dat %>% filter(diffdays == 0) # n = 66
  t2_RDISC <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_RDISC <- dat100_t1 %>% dplyr::select(matches("rdisc")) %>% cbind(demos100_t1,.) %>% filter(!is.na(rdisc_sum)) # 7 rows
  # t2_RDISC <- dat100_t2 %>% dplyr::select(matches("rdisc")) %>% cbind(demos100_t2,.) %>% filter(!is.na(rdisc_sum)) # 68 rows
  
  temp <- t2_RDISC %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_RDISC$combo <- paste(as.character(t2_RDISC$bblid),as.character(t2_RDISC$test_sessions_v.dotest),sep="_")
  t2_RDISC2 <- t2_RDISC %>% filter(combo %in% t.first$combo)  # 33 rows
  
  
  combo_RDISC <- left_join(t1_RDISC,t2_RDISC2 %>% dplyr::select(bblid,diffdays,KRDISC.valid_code:rdisc_mcr), by="bblid") %>% filter(!is.na(rdisc_sum.y)) %>%  # 2 rows
    rename(rdisc_sum_t1 = rdisc_sum.x,rdisc_sum_t2 = rdisc_sum.y)
  
  # scatters
  pdf("data/outputs/full_full/RDISC_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_RDISC %>% dplyr::select(matches("rdisc_sum")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  
  
  # VOLT -- 52 rows, need more
  bbl_keep <- repeats_only %>% filter(diffdays <= 350, test_sessions_v.age %in% 18:35) %>% pull(bblid) %>% unique()  # n = 176, realized this includes people who may have t1's < 18  yo
  bbl_keep <- extralong_repeat %>% filter(diffdays == 0, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% pull(bblid) %>% unique()  # n = 146
  
  dat <- extralong_repeat %>% filter(diffdays <= 350, test_sessions_v.age %in% 18:35, bblid %in% bbl_keep) %>% 
    arrange(bblid,test_sessions_v.dotest) %>% filter(!is.na(volt_cr)) %>%    # 172 rows
    dplyr::select(datasetid_platform:test_sessions.famid,bblid,test_sessions_v.age:platform,diffdays:combo,volt_genus:volt_w_rtcr)
  
  t1_VOLT <- dat %>% filter(diffdays == 0) # n = 66
  t2_VOLT <- dat %>% filter(diffdays > 0)  # n = 92
  
  # t1_VOLT <- dat100_t1 %>% dplyr::select(matches("volt")) %>% cbind(demos100_t1,.) %>% filter(!is.na(volt_cr)) # 51 rows
  # t2_VOLT <- dat100_t2 %>% dplyr::select(matches("volt")) %>% cbind(demos100_t2,.) %>% filter(!is.na(volt_cr)) # 62 rows
  
  temp <- t2_VOLT %>% dplyr::select("bblid","test_sessions_v.dotest")
  t.first <- temp[match(unique(temp$bblid), temp$bblid),]    # this has the first instances of each unique bblid
  t.first$combo <- paste(as.character(t.first$bblid),as.character(t.first$test_sessions_v.dotest),sep="_")
  
  t2_VOLT$combo <- paste(as.character(t2_VOLT$bblid),as.character(t2_VOLT$test_sessions_v.dotest),sep="_")
  t2_VOLT2 <- t2_VOLT %>% filter(combo %in% t.first$combo)  # 43 rows
  
  
  combo_VOLT <- left_join(t1_VOLT,t2_VOLT2 %>% dplyr::select(bblid,diffdays,volt_genus:volt_w_rtcr), by="bblid") %>% filter(!is.na(volt_cr.y)) %>%  # 17 rows
    rename(volt_cr_t1 = volt_cr.x,volt_cr_t2 = volt_cr.y)
  
  # scatters
  pdf("data/outputs/full_full/VOLT_100days_testretest_220728.pdf",height=9,width=12)
  pairs.panels(combo_VOLT %>% dplyr::select(matches("volt_cr")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
  
  }














