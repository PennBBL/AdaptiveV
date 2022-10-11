# this script comes from Mrugank's adaptive_tracker.R I've been editing it to fit my needs for AdaptiveV analysis

# Akira Di Sandro, 5.11.22

#libraries
library(data.table)
library(dplyr)
library(lubridate)
library(readr)
library(gtools)
library(tidyr)
library(matrixStats)

# function for datediff in months
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01"));lt$year*12 + lt$mon }
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }
`%notin%` <- Negate(`%in%`)
colClean <- function(x){ colnames(x) <- gsub(".1$", "", colnames(x)); x } 

# Study Enroll ----
{
  studyAll <- read.csv("data/inputs/cnb/bbl_study_all.csv", comment.char="#")  # 315 rows of unique bblid as of 9/30/22
 
  #pull out adaptive enrolled subjects
  studyEnroll<-studyAll%>%
    filter(PROTOCOL %in% c("834552 - Adaptive-V"))
  studyEnroll$BBLID = as.numeric(studyEnroll$BBLID)
  
  studyEnroll <- studyEnroll%>%
    filter( BBLID>1000, BBLID != 18026) # 18026 is an old test record according to Mrugank
  studyEnroll$formatted_doenroll<- as.Date(as.character(studyEnroll$DOENROLL), format = "%d-%b-%y")
  
  names(studyEnroll)<-tolower(names(studyEnroll))
  studyEnroll2<-studyEnroll%>%
    dplyr::select(bblid, study_status, study_group, timepoints, formatted_doenroll)   # no date for bblid 22139, says it's cross-listed as 20871, but 20871 is not on this csv, only 22867 
  
  
  # adding protocol order from what used to be adaptive_study_data_group_order_oracle.csv, now adaptive.csv
  group_order <- read.csv("data/inputs/cnb/adaptive.csv")
  # remove 22015, 22412, 22680, 23217, 91919, 94144, 112083, 119527 because they've been removed from study; 
  # 22139 is a duplicate for 20871; 
  # 23064 should exist, problem unresolved (apparently should be the same person as 23179).
  
  bblid <- group_order %>% dplyr::select(BBLID) %>% filter(BBLID > 9999) %>% unique()
  
  proto_order <- data.frame(bblid,
                            "proto_1"=NA,
                            "proto_2"=NA,
                            "proto_3"=NA,
                            "proto_4"=NA)
  for (i in 1:nrow(bblid)) {
    proto_order[proto_order$BBLID==bblid[i,1],2:5] <- group_order %>% 
      filter(BBLID==bblid[i,1]) %>% 
      distinct(ORD, .keep_all = T) %>% 
      dplyr::select(ORD, ORDER_NAME) %>% 
      arrange(ORD) %>% 
      pull(ORDER_NAME)
  }
  
  studyEnroll3 <- left_join(studyEnroll2,proto_order,by=c("bblid"="BBLID")) %>%      # 315 unique bblid rows, 9/30/22 (13 NA rows from group_order)
    filter(study_status %notin% c("excluded","dropout","not enrolled"))              # 288 rows left after getting rid of "excluded","dropout","not enrolled", 9/30/22
  
  # all rows that have data from bbl_study_all.csv, but not adaptive.csv
  bblstudyall_no_adapive <- studyEnroll3 %>% filter(is.na(proto_1))    # 13 after excluding non-active/complete ppl
}


#Demographics ----
{
  demo <- read.csv("data/inputs/cnb/subjectvisitsall_v.csv")   # still the same as most updated, 7/27/22

  names(demo)<-tolower(names(demo))
  
  #all demographics for bblids matching those that got ENROLLED
  
  demo2<-demo[demo$bblid %in% studyEnroll$bblid,]
  
  #adaptive specific demo collected
  demo3<-demo2[which(demo2$protocol=='834552 - Adaptive-V'),]
  
  
  demo3<-demo3%>%
    mutate(new_dovisit = lubridate::dmy(dovisit))%>%
    group_by(bblid)%>%
    mutate(rank = rank(desc(new_dovisit)))%>%
    filter(rank == 1)  
  
  
  ## BBLID that are in study enroll but not marked as Adaptive
  
  nodemo<-demo2[demo2$bblid %notin%demo3$bblid,] # 13, 10/5/22
  
  ## Descending Rank filters to the latest DOVISIT
  nodemo<-nodemo %>%
    mutate(new_dovisit = lubridate::dmy(dovisit))%>%
    group_by(bblid)%>%
    mutate(rank = rank(desc(new_dovisit)))%>%
    filter(rank == 1)
  
  #combine
  demofinal<-bind_rows(demo3,nodemo)
  names(demofinal)<-tolower(names(demofinal))
  #merge with study enroll
  
  # The datediff gives you the demographics collected closest to the doenroll
  
  demographics_adaptive <- left_join(studyEnroll3,demofinal, by = "bblid")      
  demographics_adaptive$age_enroll<- (mapply(mondf,as.Date(as.character(demographics_adaptive$dobirth), format = "%d-%b-%y"),demographics_adaptive$formatted_doenroll)/12)
  demographics_adaptive <- demographics_adaptive%>%
    dplyr::select(bblid, study_status, study_group, sex, race, ethnic, educ, age_enroll,proto_1,proto_2,proto_3,proto_4)

  # all rows that have data from bbl_study_all.csv, but not subjectvisitsall_v.csv
  no_subjectvisitsall <- demographics_adaptive %>% filter(is.na(sex))   # 0 missing as of 10/05/22
}


#CNB ----

{
  # cnb <- read.csv("data/inputs/cnb/cnb_merged_webcnp_surveys_allbblprjcts_longform.csv")   # old file name
  cnb <- read.csv("data/inputs/cnb/cnb_merged_webcnp_surveys_smryscores_allbbl_longform.csv") # new file name
  
  cnb1 <- cnb %>%
    dplyr::select(matches("datasetid_platform|test_sessions.datasetid|siteid|bblid.clean|age|battery|dotest|gng|cpf|medf|pvrt|er40|cpt|cpw|aim|adt|plot|volt|pmat|digsym|dscor|dsmemcr|dscorrt|dsmcrrt|KDDISC|KRDISC|EDISC"))
  
  colnames(cnb1)[colnames(cnb1) %in% c("test_sessions.datasetid", "test_sessions.siteid",
                                       "test_sessions.bblid.clean", "test_sessions_v.age",            
                                       "test_sessions_v.battery",           
                                       "test_sessions_v.dotest" )] <-
    c("datasetid", "siteid",       
      "bblid", "age","battery","dotest")
  
  # need to make sure bblid 23064 --> 23179
  cnb1[which(cnb1$bblid == 23064),"bblid"] <- 23179
  
  cnb2 <- cnb1[cnb1$bblid %in% studyEnroll3$bblid, ] %>% filter(siteid == "adaptive_v")
  
  # inserted DISC summary score script
  ddisc_qs <- cnb2 %>% dplyr::select((names(cnb2)[grep('DDISC.q',colnames(cnb2))]))
  cnb2$ddisc_sum <- rowSums(ddisc_qs-1)
  ddisc_ttrs <- cnb2 %>% dplyr::select((names(cnb2)[grep('DDISC.trr',colnames(cnb2))]))
  cnb2$ddisc_mcr <- rowMedians(as.matrix(ddisc_ttrs))
  
  rdisc_qs <- cnb2 %>% dplyr::select((names(cnb2)[grep('RDISC.q',colnames(cnb2))]))
  cnb2$rdisc_sum <- rowSums(rdisc_qs[,1:41]-1)
  rdisc_ttrs <- cnb2 %>% dplyr::select((names(cnb2)[grep('RDISC.trr',colnames(cnb2))]))
  cnb2$rdisc_mcr <- rowMedians(as.matrix(rdisc_ttrs[,1:41]))
  
  edisc_qs <- cnb2 %>% dplyr::select((names(cnb2)[grep('EDISC',colnames(cnb2))])) %>% 
    dplyr::select((names(cnb2)[grep('_resp',colnames(cnb2))])) %>% 
    dplyr::select(EDISC.q_101_resp:EDISC.q_134_resp)
  cnb2$edisc_sum <- rowSums(edisc_qs-1)
  edisc_ttrs <- cnb2 %>% dplyr::select((names(cnb2)[grep('EDISC',colnames(cnb2))])) %>% 
    dplyr::select((names(cnb2)[grep('_ttr',colnames(cnb2))])) %>% 
    dplyr::select(EDISC.q_101_ttr:EDISC.q_134_ttr)
  cnb2$edisc_mcr <- rowMedians(as.matrix(edisc_ttrs))
  
  # keep all tasks, get rid of itemwise DISC cols
  fullcnb <- cnb2 %>% dplyr::select(datasetid_platform:aim_mcrrt,matches("ddisc_|rdisc_|edisc_")) %>%  # 278 rows as of 10/03/22
    filter(datasetid %notin% c(47859,48505,52277))     # get rid of doubled up bblid: datasetid = 47859, datasetid = 48505, datasetid = 52277
  
  
  # adding PRA from iw
  for_pra <- read.csv("data/inputs/cnb/athena_195_360.csv",na.strings=c(""," ","NA"))
  
  PRA_iw <- for_pra %>% mutate(bblid = as.numeric(test_sessions_v.bblid)) %>% arrange(bblid) %>% 
    filter(test_sessions.siteid == "adaptive_v",test_sessions_v.battery == "PRA_D", bblid > 9999,!is.na(test_sessions_v.battery_complete),test_sessions.datasetid!=48671) %>% 
    dplyr::select(matches("test_session|^bblid|PRA_D")) %>% dplyr::select(test_sessions.datasetid:test_sessions.famid,bblid,test_sessions_v.admin_comments:test_sessions_v.batt_consent,test_sessions_v.deleted_flag:PRA_D.AGE)
  
  PRA_resp <- PRA_iw %>% dplyr::select(matches("RESP"))
  PRA_resp[PRA_resp == 2] <- 0
  PRA_iw$pra_cr <- rowSums(PRA_resp,na.rm = T)
  
  fullcnb2 <- left_join(fullcnb,PRA_iw %>% dplyr::select(bblid,pra_cr),by="bblid") # 11 NA as of 10/03/22, duplicate of 20189 and we should only keep the 2022 record 9/6/22
  
  # subset (old, from Mrugank's original script)
  cnb2_subset <- cnb2 %>%
    dplyr::select(bblid, dotest, battery) %>%
    distinct(bblid, dotest, battery) %>%
    mutate(status = 1) %>%
    pivot_wider(names_from = battery, values_from = status)

  # all rows that have data from bbl_study_all.csv, but not cnb_merged_webcnp_surveys_allbblprjcts_longform.csv 
  missing <- studyEnroll3[studyEnroll3$bblid %notin% cnb2$bblid,]  # 10 as of 10/5/22
  
  tracker <- left_join(demographics_adaptive,cnb2_subset, by = c("bblid"))
  dat_combined <- left_join(demographics_adaptive,fullcnb2, by = c("bblid"))  
}




#Adaptive CNB ----

{
  adaptive_v2 <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v2_20220926_095354.csv")
  adaptive_v <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v_20220909_221240.csv")
  adaptive_cpfv2_er40v2 <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v_cpfv2_er40v2_20220926_095401.csv")
  adaptive_prad <- read_csv("data/inputs/cnb/CNB_CAT_session_pra-d_20220926_095411.csv")
  
  adaptive_v2 <- rename(adaptive_v2, Dataset.ID = 1) %>% 
    mutate(version = "adaptive_cnb_V2")
  adaptive_v <- rename(adaptive_v, Dataset.ID = 1) %>% 
    mutate(version = "adaptive_cnb_V")
  adaptive_cpfv2_er40v2 <- rename(adaptive_cpfv2_er40v2, Dataset.ID = 1) %>% 
    mutate(version = "adaptive_cnb_V2_cpf_er40")
  adaptive_prad <- rename(adaptive_prad, Dataset.ID = 1) %>% 
    mutate(version = "adaptive_cnb_pra_d")
  
  # clean and get ready to merge
  adaptive_prad <- adaptive_prad %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_pra=Dataset.ID,date_pra=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  adaptive_prad$uniqtest <- paste(adaptive_prad$testcode,adaptive_prad$decisiontreename,sep="_")
  adaptive_prad <- adaptive_prad %>% dplyr::select(BBLID,datasetid_pra,date_pra,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID) %>% mutate(`pra-1.00-d-cat_default` = as.numeric(`pra-1.00-d-cat_default`))
  
  # changing BBLID 22668 -> 13446 for datsetid == 432, 6.7.22
  # adaptive_prad[which(adaptive_prad$BBLID==22668 & adaptive_prad$datasetid_pra == 432),"BBLID"] <- 13446
  
  # BBLID 90922 has 4 different data entries for pra_d, only keep datasetid == 291 for now because that matches adaptive_v data, all other datasetids match up (with dotest as well) with some people in adpative_v2
  # BBLID 90922, datasetid == 291 --> keep the same
  # BBLID 90922, datasetid == 489 --> ignore for now
  # BBLID 90922, datasetid == 686 --> 95116 (according to matching dotest + almost consecutive datasetid, 684)
  # BBLID 90922, datasetid == 693 --> ignore for now
  adaptive_prad[which(adaptive_prad$datasetid_pra == 686),"BBLID"] <- 95116
  
  # also remove BBLID 94144 (pt no-showed 3x, excluded from study)
  adaptive_prad <- adaptive_prad %>% filter(datasetid_pra %notin% c(489,693,682))
  
  
  adaptive_cpfv2_er40v2 <- adaptive_cpfv2_er40v2 %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v2CE=Dataset.ID,date=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999,datasetid_v2CE!=893) # still deciding what to do with double data for 15507, for now remove second, aka datasetid_v2CE==893 
  adaptive_cpfv2_er40v2$uniqtest <- paste(adaptive_cpfv2_er40v2$testcode,adaptive_cpfv2_er40v2$decisiontreename,sep="_")
  adaptive_cpfv2_er40v2 <- adaptive_cpfv2_er40v2 %>% dplyr::select(BBLID,datasetid_v2CE,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  
  
  adaptive_v <- adaptive_v %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v=Dataset.ID,date=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  adaptive_v$uniqtest <- paste(adaptive_v$testcode,adaptive_v$decisiontreename,sep="_")
  adaptive_v <- adaptive_v %>% dplyr::select(BBLID,datasetid_v,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  
  # get rid of datasetid == 322 because it was only partially completed and BBLID 22454 has complete record in datasetid == 323
  # get rid of datasetid == 543 because it was only partially completed and BBLID 90158 has a more complete record in datasetid == 542
  adaptive_v <- adaptive_v %>% filter(datasetid_v %notin% c(322,543))
  
  
  adaptive_v2 <- adaptive_v2 %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v2=Dataset.ID,date=`Response Logged Date`,score=`Running Score`) %>% 
    filter(`Final Score`=="Yes",BBLID > 9999) 
  adaptive_v2$uniqtest <- paste(adaptive_v2$testcode,adaptive_v2$decisiontreename,sep="_")
  adaptive_v2 <- adaptive_v2 %>% dplyr::select(BBLID,datasetid_v2,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  
  # replace temp with adaptive_v2[,4:ncol(adaptive_v2)]
  adaptive_v2[,4:ncol(adaptive_v2)] <- adaptive_v2 %>% dplyr::select(matches(".00")) %>% data.frame(sapply(., as.numeric)) %>% dplyr::select(matches(".1$")) %>% colClean(.)
  
  
  
  
  
  # merging code
  # idea is to merge pra + link1 (then merge with cpf_2/er40_2), pra + link2 then combine these columns (adding NA columns for setdiff)
  pra_v <- inner_join(adaptive_prad,adaptive_v,by=c("BBLID","date_pra" = "date")) %>% arrange(BBLID)
  
  pra_v_new_cpf_er40 <- left_join(pra_v,adaptive_cpfv2_er40v2,by=c("BBLID")) %>% arrange(BBLID)
  
  pra_v1_tocombine <- data.frame(pra_v_new_cpf_er40 %>% dplyr::select(BBLID:datasetid_v),datasetid_v2 = NA,
                                 pra_v_new_cpf_er40 %>% dplyr::select(`medf-1.00-cat_same`:`er40-2.00-cat_neutral`),
                                 `cpf2.1.00.v1.cat_target` = NA, `cpf2.1.00.v1.cat_foil` = NA)
  
  pra_v2 <- inner_join(adaptive_prad,adaptive_v2,by=c("BBLID","date_pra" = "date")) %>% arrange(BBLID)
  
  pra_v2_tocombine <- data.frame(pra_v2 %>% dplyr::select(BBLID:`pra-1.00-d-cat_default`),
                                 datasetid_v = NA,pra_v2 %>% dplyr::select(datasetid_v2),
                                 pra_v2 %>% dplyr::select(`medf-1.00-cat_same`:`pvrt-1.00-cat_default`),
                                 `cpf-1.00-v1-cat_target` = NA, `cpf-1.00-v1-cat_foil` = NA,
                                 `er40-1.00-cat_emotive` = NA, `er40-1.00-cat_neutral` = NA,
                                 pra_v2 %>% dplyr::select(`cpw-1.00-v1-cat_target`:`edisc-1.00-cat_default`),
                                 datasetid_v2CE = NA, date = NA, 
                                 `cpf2-1.00-v2-cat_target` = NA, `cpf2-1.00-v2-cat_foil` = NA,
                                 pra_v2 %>% dplyr::select(`er40-2.00-cat_emotive`:`er40-2.00-cat_neutral`),
                                 pra_v2 %>% dplyr::select(`cpf2-1.00-v1-cat_target`:`cpf2-1.00-v1-cat_foil`))
  
  # some people who have PRA who don't have adaptive_v or adaptive_v2 data
  pra_no_other <- adaptive_prad %>% filter(BBLID %notin% pra_v1_tocombine$BBLID,BBLID %notin% pra_v2_tocombine$BBLID)
  
  adaptive_cnb <- rbind(pra_v1_tocombine,pra_v2_tocombine) # 277 as of 10/03/22 (with some tests having less)
  
  
  
  # testing out new CAT CNB short tests
  adaptive_short <- read_csv("data/inputs/cnb/athena_253_324.csv")
  
  # there's data now! 9/6/22
  adaptive_AIM_CPT <- adaptive_short %>% mutate(bblid = as.numeric(test_sessions.bblid)) %>% 
    filter(test_sessions.siteid == "adaptive_v",bblid>9999,!is.na(CPT108.valid_code) | !is.na(S_AIM.valid_code)) %>% arrange(bblid) %>% # BBLID 111789: (upon discussion w Kelly, keep datasetid 48932 and change sex to F since AFAB)
    dplyr::select(matches("test_sessions|^bblid|S_AIM|CPT108")) %>% rename(datasetid_AC = test_sessions.datasetid) %>% filter(datasetid_AC != 48930)  # had to switch CPT108.ScorVers for  CPT108.valid_code, same for AIM on 10/03/2022
  
  adaptive_DIGSYM_GNG <- adaptive_short %>% mutate(bblid = as.numeric(test_sessions.bblid)) %>% 
    filter(test_sessions.siteid == "adaptive_v",bblid>9999,!is.na(GNG60.valid_code) | !is.na(S_DIGSYM.valid_code)) %>% arrange(bblid) %>% 
    dplyr::select(matches("test_sessions|^bblid|S_DIGSYM|GNG60")) %>% rename(datasetid_DG = test_sessions.datasetid)
  adaptive_DIGSYM_GNG[which(adaptive_DIGSYM_GNG$datasetid_DG == 48900),20:33] <- adaptive_DIGSYM_GNG[which(adaptive_DIGSYM_GNG$datasetid_DG == 48898),20:33] # BBLID: 104265 (GNG and DIGSYM data separated into two datasetids because of technical issues, good to combine) from CAT CNB notes
  adaptive_DIGSYM_GNG <- adaptive_DIGSYM_GNG %>% filter(datasetid_DG != 48898)
  
  # combine AIM/CPT and DS/GNG people, select only relevant columns
  adaptive_short_tojoin <- left_join(adaptive_AIM_CPT %>% dplyr::select(matches("^bblid|datasetid_AC|CPT108.CATCPTT_TP|CPT108.CATCPTT_FP|CPT108.CATCPTT_TPRT$|CPT108.CATCPTT_FPRT$|S_AIM.AIM_NM|S_AIM.AIM_M|S_AIM.CRRT_NM|S_AIM.CRRT_M|S_AIM.AIMTOT|S_AIM.AIMTOTRT")),
                                     adaptive_DIGSYM_GNG %>% dplyr::select(matches("^bblid|datasetid_DG|S_DIGSYM.DSCOR|S_DIGSYM.DSCORRT|S_DIGSYM.DSMEMCR|S_DIGSYM.DSMCRRT|GNG60.GNG60_CR|GNG60.GNG60_RTCR")),by="bblid") %>% 
    dplyr::select(bblid,datasetid_AC,S_AIM.AIM_NM:S_AIM.AIMTOTRT,CPT108.CATCPTT_FP:CPT108.CATCPTT_TPRT,CPT108.CATCPTT_FPRT,datasetid_DG,S_DIGSYM.DSCOR:S_DIGSYM.DSMCRRT,GNG60.GNG60_CR:GNG60.GNG60_RTCR)
  
  # 20189 23402 23486 23536 23584 96149 in adaptive_short_tojoin, but not merged because other non-short test data doesn't exist
  # 20189,23536,23584,96149 not in adpative_pra either
  
  adaptive_cnb1 <- left_join(adaptive_cnb,adaptive_short_tojoin,by=c("BBLID"="bblid")) # 277 rows as of 10/3/22
  
  # 199 unique BBLIDs in cat_dat_nodates:
  #   1) BBLID 90922, datasetid_pra 489 only has pra_d
  #   2) BBLID 90158, has one pra_d but two records of every other test on adaptive_v, both incomplete
  #   3) BBLID 22668, datasetid_pra 432 only has pra_d
  #   4) BBLID 22454, has one incomplete adaptive_v battery (datasetid_v 322) and one complete adaptive_v battery (datasetid_v 323) with the same pra_d info
  
  
  dat_combined2 <- left_join(dat_combined,adaptive_cnb1, by = c("bblid" = "BBLID"))    # 100079 is excluded and still trying to figure out about 23064
  
  # all rows that have data from bbl_study_all.csv, but not in the adaptivecnb CSVs  
  no_adaptivecnb <- dat_combined2 %>% filter(is.na(pra.1.00.d.cat_default)) # 13 rows missing (at least), 10/3/22
  
}


# print out CSV of combined CNB data
write.csv(dat_combined2,"data/inputs/cnb_merged/cnb_merged_20221005.csv",row.names = F)







#Phenome ----

{ # using Mrugank's merged files, 9/30/22
  # 301 unique BBLIDs
  CATGOA_ext <- read.csv("data/inputs/goa/externalizing_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled"))    # n = 273 unique BBLID, duplicate BBLIDs 106255
  CATGOA_mood <- read.csv("data/inputs/goa/mood_anxiety_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled"))    # n = 276 unique BBLID, duplicate BBLIDs 22196,23336,100056,106255,107712,114007
  CATGOA_per <- read.csv("data/inputs/goa/personality_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled"))    # n = 273 unique BBLID, duplicate BBLIDs 106255
  CATGOA_phob <- read.csv("data/inputs/goa/phobia_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled"))    # n = 274 unique BBLID, duplicate BBLIDs 22196,106255,107712
  CATGOA_psy <- read.csv("data/inputs/goa/psychosis_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled"))    # n = 273 unique BBLID, duplicate BBLIDs 106255
  
  # Digging on Slack and checking other notes to figure out these duplicates
  # 106255: unclear, two administrations, 4 months apart
  # 
    
  # figuring out why there are duplicates
  
  phenome <- bind_rows(CATGOA_ext,CATGOA_mood,CATGOA_per,CATGOA_phob,CATGOA_psy)
  
  phenome_subset <- phenome %>%
    dplyr::select(bblid:Date,Score) %>%
    filter(!is.na(as.numeric(bblid))) %>%
    filter(bblid>1000) %>%
    pivot_wider(names_from = Domain, values_from = Score)
  
  phenome_subset$bblid<- as.numeric(phenome_subset$bblid)
  names(phenome_subset)[10:14] <- c("Externalizing","Mood","Personality","Phobias","Psychosis")
  
  # unlist the scores
  temp_ext <- unlist(phenome_subset$Externalizing)
  temp_moo <- unlist(phenome_subset$Mood)
  temp_per <- unlist(phenome_subset$Personality)
  temp_pho <- unlist(phenome_subset$Phobias)
  temp_psy <- unlist(phenome_subset$Psychosis)
  
  
  
  # trying something else for now
  CATGOA_ext <- read.csv("data/inputs/goa/externalizing_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled")) %>% 
    rename(CAText_score=Score) # n = 273 unique BBLID, duplicate BBLIDs 106255
  CATGOA_mood <- read.csv("data/inputs/goa/mood_anxiety_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled")) %>% 
    rename(CATmood_score=Score)    # n = 276 unique BBLID, duplicate BBLIDs 22196,23336,100056,106255,107712,114007
  CATGOA_per <- read.csv("data/inputs/goa/personality_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled")) %>% 
    rename(CATper_score=Score)    # n = 273 unique BBLID, duplicate BBLIDs 106255
  CATGOA_phob <- read.csv("data/inputs/goa/phobia_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled")) %>% 
    rename(CATphob_score=Score)    # n = 274 unique BBLID, duplicate BBLIDs 22196,106255,107712
  CATGOA_psy <- read.csv("data/inputs/goa/psychosis_adaptive.csv") %>% arrange(bblid) %>% 
    filter(!is.na(Domain),study_status %notin% c("dropout","excluded","not enrolled")) %>% 
    rename(CATpsy_score=Score)    # n = 273 unique BBLID, duplicate BBLIDs 106255
  
  phenome_1 <- left_join(CATGOA_ext %>% dplyr::select(bblid:educ,Date:Time,CAText_score),CATGOA_mood %>% dplyr::select(bblid,CATmood_score),by="bblid")
  phenome_1 <- left_join(phenome_1,CATGOA_per %>% dplyr::select(bblid,CATper_score),by="bblid")
  phenome_1 <- left_join(phenome_1,CATGOA_phob %>% dplyr::select(bblid,CATphob_score),by="bblid")
  phenome_1 <- left_join(phenome_1,CATGOA_psy %>% dplyr::select(bblid,CATpsy_score),by="bblid")
  
  
  # still need to work through the following duplicate BBLIDs: 16527,22593,106255,107712
  
  # all(phenome_subset$bblid %in% demographics_adaptive$bblid) TRUE
  phenome_tojoin <- left_join(demographics_adaptive %>% dplyr::select(bblid,proto_1:proto_2),phenome_1,by="bblid") %>% 
    dplyr::select(bblid:Date,proto_1:proto_2,CAText_score:CATpsy_score)
}

#GOASSESS Full
{
  common_interview<- read.csv("data/inputs/goa/Goassess.csv") %>% filter(bblid>9999)
  x <- common_interview
  
  ext <- data.frame(
    x$bblid,
    x$redcapid,
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
    x$add013_short,
    x$add014_short,
    x$add022_short,
    x$add020_short,
    x$add015_short,
    x$cdd011) 
  
  mood <- data.frame(
    x$bblid,
    x$redcapid,
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
    x$dep004_short,
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
  
  # missing info for now
  per <- data.frame(
    x$bblid,
    x$redcapid
  )
  
  fear <- data.frame(
    x$bblid,
    x$redcapid,
    x$agr006,
    x$agr005,
    x$soc004,
    x$agr008,
    x$soc003,
    x$agr004,
    x$agr001,
    x$soc005,
    x$soc001,
    x$phb004,
    x$agr002,
    x$agr003,
    x$phb007,
    x$soc002,
    x$phb006,
    x$phb001,
    x$agr007,
    x$phb002,
    x$sep509,
    x$phb003,
    x$phb005,
    x$sep508,
    x$phb008,
    x$sep500,
    x$sep511)
  
  psy <- data.frame(
    x$bblid,
    x$redcapid,
    x$sip012,
    x$sip007,
    x$sip010,
    x$sip008,
    x$sip011,
    x$sip013,
    x$sip005,
    x$sip003,
    x$sip004,
    x$sip006,
    x$sip009,
    x$psy001,
    x$psy029,
    x$sip014,
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
  
  # make sum scores
  ext <- cbind(ext[,1:2],mutate_all(ext[,3:ncol(ext)], function(x) as.numeric(as.character(x))))
  ext[ext == 9] <- NA   # 9 is for unknown
  ext <- ext %>% mutate(ext_sum = rowSums(.[,3:ncol(.)],na.rm=T))
  
  mood[mood == 9] <- NA   # 9 is for unknown
  mood <- mood %>% mutate(mood_sum = rowSums(.[,3:ncol(.)],na.rm=T))
  
  fear[fear == 9] <- NA   # 9 is for unknown
  fear <- fear %>% mutate(fear_sum = rowSums(.[,3:ncol(.)],na.rm = T))
  
  psy[psy == 9] <- NA   # 9 is for unknown
  psy <- psy %>% mutate(psy_sum = rowSums(.[,3:ncol(.)],na.rm = T))
  
  # combine all sum scores into one dataframe
  fullGOA <- left_join(ext %>% dplyr::select(x.bblid,x.redcapid,ext_sum),
                       mood %>% dplyr::select(x.bblid,x.redcapid,mood_sum),by=c("x.bblid","x.redcapid"))
  fullGOA <- left_join(fullGOA,fear %>% dplyr::select(x.bblid,x.redcapid,fear_sum),by=c("x.bblid","x.redcapid"))
  fullGOA <- left_join(fullGOA,psy %>% dplyr::select(x.bblid,x.redcapid,psy_sum),by=c("x.bblid","x.redcapid"))  # 305 unique BBLIDs, 10/06/22
  names(fullGOA)[1:2] <- c("bblid","redcapid")
  
  # combine Full GOA scores with CAT GOA scores + demos
  allGOA <- left_join(fullGOA,phenome_tojoin,by="bblid") # 345 rows because of the duplicate BBLIDs, 10/06/22
  # can't fix this in time right now so I'm just going to pick random records for the duplicates
  allGOA <- allGOA %>% arrange(bblid)
  allGOA <- allGOA[-c(173,98:100,250,267:297,301:303,314),]
}


write.csv(allGOA,"data/inputs/goa/GOA_merged_221006.csv", row.names = FALSE)






# extra ----

# full GOA records to check, initial QA with Ally
common_interview <- read.csv("data/inputs/goa/axis_common_interview_220829.csv") %>%
  dplyr::select(bblid,assessment, interview_complete, interview_complete_date, interview_complete_notes)%>%
  filter(assessment == 'PNC GOASSESS(GAF Mod)' )

tocheck <- common_interview %>% filter(is.na(interview_complete) | interview_complete == 0) %>% 
  arrange(interview_complete,bblid)

write.csv(tocheck,"data/outputs/fullGOA_tocheck.csv",row.names = F)


