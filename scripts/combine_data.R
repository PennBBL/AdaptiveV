# initial data analysis for CNB (started 02.04.22)

# Aki and Kat


# load packages ----
library(readr)
library(tidyr)
library(dplyr)
library(matrixStats)



# functions ----
colClean <- function(x){ colnames(x) <- gsub(".1$", "", colnames(x)); x } 


# Group data long to wide ----
{
  dat1 <- read_csv("data/inputs/adaptive_study_data_group_order_oracle_030522.csv")
  
  dat1$ps <- ifelse(dat1$GRP=="Psychosis",1,0)
  dat1$td <- ifelse(dat1$GRP=="Healthy Controls",1,0)
  dat1$mood <- ifelse(dat1$GRP=="Mood-Anx-BP",1,0)
  
  temp <- dat1[,c(1:2,10:12)]
  
  bblid <- unique(dat1$BBLID)
  
  newdat1 <- data.frame(matrix(rep(NA,5),nrow=1))
  names(newdat1) <- names(temp)
  for (id in bblid) {
    newrow <- temp[which(temp$BBLID==id),]
    if (!(id %in% newdat1$BBLID)){
      newdat1 <- rbind(newdat1,newrow[1,])
    }
  }
  newdat1 <- newdat1 %>% filter(BBLID>9999)
  
  # getting protocol order by a loop right now, but maybe make this into a function?
  proto_order <- data.frame(bblid,
                            "proto_1"=NA,
                            "proto_2"=NA,
                            "proto_3"=NA,
                            "proto_4"=NA)
  for (bbl in bblid) {
    proto_order[proto_order$bblid==bbl,2:5] <- dat1 %>% 
      filter(BBLID==bbl) %>% 
      distinct(ORD, .keep_all = T) %>% 
      dplyr::select(ORD, ORDER_NAME) %>% 
      arrange(ORD) %>% 
      pull(ORDER_NAME)
  }
  
  grp_order <- left_join(newdat1,proto_order,by=c("BBLID"="bblid"))
}





# CAT CNB long to wide ----

# dat2 <- read_csv("data/inputs/cnb_cat_out.csv")
# dat2 <- data.frame(dat2)
# dat2$uniqtest <- paste(dat2$test_code,dat2$Decision_Tree_Collection_Name,sep="_")
# 
# newdat2 <- pivot_wider(dat2[,c(1,2,9,10,11)], 
#                         names_from = "uniqtest",
#                         values_from = "accuracy_score")

# for newest cat_out:
{
  cat_pra <- read_csv("data/inputs/CNB_CAT_session_pra-d_20220428_125704.csv")
  cat_v2_cpf_er40 <- read_csv("data/inputs/CNB_CAT_session_adaptive_v_cpfv2_er40v2_20220428_125641.csv")
  cat_v <- read_csv("data/inputs/CNB_CAT_session_adaptive_v_20220428_125332.csv")
  cat_v2 <- read_csv("data/inputs/CNB_CAT_session_adaptive_v2_20220428_125604.csv")
  
  cat_pra <- cat_pra %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_pra=`Dataset ID`,date_pra=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  cat_pra$uniqtest <- paste(cat_pra$testcode,cat_pra$decisiontreename,sep="_")
  cat_pra <- cat_pra %>% dplyr::select(BBLID,datasetid_pra,date_pra,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID) %>% mutate(`pra-1.00-d-cat_default` = as.numeric(`pra-1.00-d-cat_default`))
  
  cat_v2_cpf_er40 <- cat_v2_cpf_er40 %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v2CE=`Dataset ID`,date=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  cat_v2_cpf_er40$uniqtest <- paste(cat_v2_cpf_er40$testcode,cat_v2_cpf_er40$decisiontreename,sep="_")
  cat_v2_cpf_er40 <- cat_v2_cpf_er40 %>% dplyr::select(BBLID,datasetid_v2CE,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  
  cat_v <- cat_v %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v=`Dataset ID`,date=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  cat_v$uniqtest <- paste(cat_v$testcode,cat_v$decisiontreename,sep="_")
  cat_v <- cat_v %>% dplyr::select(BBLID,datasetid_v,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  
  cat_v2 <- cat_v2 %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v2=`Dataset ID`,date=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  cat_v2$uniqtest <- paste(cat_v2$testcode,cat_v2$decisiontreename,sep="_")
  cat_v2 <- cat_v2 %>% dplyr::select(BBLID,datasetid_v2,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  # replace temp with cat_v2[,4:ncol(cat_v2)]
  cat_v2[,4:ncol(cat_v2)] <- cat_v2 %>% select(matches(".00")) %>% data.frame(sapply(., as.numeric)) %>% select(matches(".1$")) %>% colClean(.)
  
  
  cat_dat <- full_join(cat_pra,cat_v,by=c("BBLID","date_pra" = "date")) %>% arrange(BBLID)
  cat_dat_nodates <- full_join(cat_dat,cat_v2_cpf_er40,by=c("BBLID")) %>% arrange(BBLID) %>%    #ignoring dates to combine cat_pra, cat_v, and cat_v2cpf_er40
    mutate(`cpf2-1.00-v1-cat_target` = NA, `cpf2-1.00-v1-cat_foil` = NA, datasetid_v2 = NA)
  
  # 199 unique BBLIDs in cat_dat_nodates:
  #   1) BBLID 90922, datasetid_pra 489 only has pra_d
  #   2) BBLID 90158, has one pra_d but two records of every other test on adaptive_v, both incomplete
  #   3) BBLID 22668, datasetid_pra 432 only has pra_d
  #   4) BBLID 22454, has one incomplete adaptive_v battery (datasetid_v 322) and one complete adaptive_v battery (datasetid_v 323) with the same pra_d info
  
  cnb_all <- full_join(cat_dat_nodates,cat_v2,by=c("BBLID",colnames(cat_v2)[4:21])) %>% arrange(BBLID)
  
  duplicates <- data.frame(table(cnb_all$BBLID)) %>% filter(Freq == 2) %>% pull(Var1)
  cnb_all_dup <- cnb_all %>% filter(BBLID %in% duplicates)

}




# Full CNB data, merge with other data ----


# dat3 <- read_csv("data/cnb_data_20220126.csv")
# dat3 <- dat3[,-6]
# names(dat3)[6] <- "bblid"

# full_dat3 <- dat3 %>% 
#   subset(test_sessions_v.battery!="cat_adaptive_v_battery_link3")
# 
# cat_dat3 <- dat3 %>% 
#   subset(test_sessions_v.battery=="cat_adaptive_v_battery_link3") %>% 
#   rename(cat_digsym_genus = digsym_genus) %>% 
#   rename(cat_dscor = dscor) %>% 
#   rename(cat_dsmemcr = dsmemcr) %>% 
#   rename(cat_dscorrt = dscorrt) %>% 
#   rename(cat_dsmcrrt = dsmcrrt) %>% 
#   rename(cat_digsym_valid = digsym_valid) %>% 
#   select(bblid, test_sessions.subid,matches("cat_"))
# 
# newdat3 <- full_join(full_dat3,cat_dat3,by="bblid")

# everything related to dat3 above this line was for the old CSV, cnb_data_20220126.csv

# dat3 <- read_csv("data/inputs/cnb_merged_20220311.csv")
# dat3 <- dat3[,-7]
# names(dat3)[6] <- "bblid"
# dat3 <- dat3 %>%
#   filter(test_sessions.siteid == "adaptive_v") %>% 
#   filter(bblid > 9999) %>% 
#   dplyr::select(!(matches("DISC") | matches("mpraxis") | matches("pcet") | matches("lnb") | matches("tap")))

{
dat05 <- read_csv("data/inputs/cnb_merged_webcnp_surveys_allbblprjcts_longform.csv")
dat05 <- dat05[,-7]
names(dat05)[6] <- "bblid"
dat05 <- dat05 %>%
  filter(test_sessions.siteid == "adaptive_v") %>% 
  filter(bblid > 9999) %>% 
  dplyr::select(!(matches("mpraxis") | matches("pcet") | matches("lnb") | matches("tap") | (names(dat05)[grep('^KRDISC.*\\.1$',colnames(dat05))])))


# * DISC task sum scores ----
ddisc_qs <- dat05 %>% dplyr::select((names(dat05)[grep('DDISC.q',colnames(dat05))]))
dat05$ddisc_sum <- rowSums(ddisc_qs-1)
ddisc_ttrs <- dat05 %>% dplyr::select((names(dat05)[grep('DDISC.trr',colnames(dat05))]))
dat05$ddisc_mcr <- rowMedians(as.matrix(ddisc_ttrs))

rdisc_qs <- dat05 %>% dplyr::select((names(dat05)[grep('RDISC.q',colnames(dat05))]))
dat05$rdisc_sum <- rowSums(rdisc_qs-1)
rdisc_ttrs <- dat05 %>% dplyr::select((names(dat05)[grep('RDISC.trr',colnames(dat05))]))
dat05$rdisc_mcr <- rowMedians(as.matrix(rdisc_ttrs))

edisc_qs <- dat05 %>% dplyr::select((names(dat05)[grep('EDISC',colnames(dat05))])) %>% 
  dplyr::select((names(dat05)[grep('_resp',colnames(dat05))])) %>% 
  dplyr::select(EDISC.q_101_resp:EDISC.q_134_resp)
dat05$edisc_sum <- rowSums(edisc_qs-1)
edisc_ttrs <- dat05 %>% dplyr::select((names(dat05)[grep('EDISC',colnames(dat05))])) %>% 
  dplyr::select((names(dat05)[grep('_ttr',colnames(dat05))])) %>% 
  dplyr::select(EDISC.q_101_ttr:EDISC.q_134_ttr)
dat05$edisc_mcr <- rowMedians(as.matrix(edisc_ttrs))
fullcnb <- dat05 %>% select(!(matches("KDDISC.") | matches("KRDISC.") | matches("EDISC.q") | matches("EDISC.t") | matches("EDISC.v")))
}


# combine data ----


bbl_grp <- unique(grp_order$BBLID)
bbl_cat <- unique(cat_dat_nodates$BBLID)
bbl_full <- unique(fullcnb$bblid)

com_bbl <- intersect(bbl_grp,bbl_cat)
com_bbl <- intersect(bbl_full,com_bbl)    # 185 BBLIDs, 5.4.22

dat <- left_join(fullcnb,grp_order,by=c("bblid"="BBLID"))
dat <- left_join(dat,cat_dat_nodates,by=c("bblid"="BBLID"))

# write_csv(dat, "data/full_CAT_cnb.csv")

# adding short tests ----
dat4 <- read_csv("data/inputs/cnb_rawscores_short_tests.csv") %>%      # there's one BBLID that has two rows, getting rid of one arbitrarily for now
  subset(subid!="111789_4")

# temp <- dat4[dat4$bblid==111789,]
# write.csv(temp,"data/duplicate_bblid.csv",row.names = F)

dat <- left_join(dat,dat4[,-c(2:5)],by="bblid")
write.csv(dat,"data/outputs/full_CAT_cnb_4May22.csv",row.names = F,na="")


