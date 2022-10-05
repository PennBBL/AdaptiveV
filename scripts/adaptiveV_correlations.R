# using Tyler's code from cnb-cat_test-retest_24march2022.r with most up to date data

# Akira Di Sandro, 5.4.22


`%notin%` <- Negate(`%in%`)
count_LRNR <- function(x){
  resps <- x %>% dplyr::select(matches("_RESP$"))
  lrnr <- c()
  
  for (i in 1:nrow(resps)){
    max_count <- 0
    if(rowSums(is.na(resps[i,]))==ncol(resps)){
      lrnr[i] <- NA
    } else {
      for (j in 1:ncol(resps)) {
        count <- 0
        if(resps[i,j] == 0){
          done <- F
          while (resps[i,j] == 0 & !(done)) {
            count <- count+1
            j <- j+1
            if(j>ncol(resps)){
              j <- ncol(resps)
              done <- T}
          }
        }
        max_count <- max(max_count,count)
      }
      lrnr[i] <- max_count
    }
  }
  return(lrnr)
}
# load packages ----

library(psych)
library(qgraph)
library(dplyr)
library(tidyr)
library(PerFit)
library(ggplot2)
library(ggdist)
library(sdamr)
library(kableExtra)
library(matrixStats)
library(irr)
library(lubridate)
library(readr)


# load data ----

# demos, proto order, groups, full cnb + cat cnb (missing rapid tests aka AIM, CPT, GNG, DIGSYM)
all_cnb <- read.csv("data/inputs/cnb_merged/cnb_merged_20221005.csv",na.strings=c(""," ","NA"))  # 288 rows as of 10/5/22
x <- all_cnb

# full CNB for itemwise DISC tasks
# f_cnb <- read.csv("data/inputs/cnb/cnb_merged_webcnp_surveys_allbblprjcts_longform.csv",na.strings=c(""," ","NA"))  # old name
f_cnb <- read.csv("data/inputs/cnb/cnb_merged_webcnp_surveys_smryscores_allbbl_longform.csv",na.strings=c(""," ","NA"))  # new name same file, just more updated, 10/3/2022
f_cnb <- f_cnb %>% filter(test_sessions_v.battery %in% c("adaptive_v_battery_v0","adaptive_v_battery_v1","adaptive_v_battery_v2"),test_sessions.bblid.clean>9999)
demo_f_cnb <- f_cnb %>% dplyr::select(datasetid_platform:test_sessions.bblid.clean,test_sessions_v.age:test_sessions_v.endtime) %>% 
  rename(BBLID = test_sessions.bblid.clean)  # 281 rows as of 10/3/22

# old itemwise full CNB
{
  old_iw <- read.csv("data/inputs/adaptive_cnb_battery_itemlevel_report_335.csv",na.strings=c(""," ","NA"))
  old_iw <- old_iw %>% filter(test_sessions.datasetid %notin% c(49043, 49044, 47906, 47859,48036))  # 5 duplicates for 4 different bblid's
  
  fullcnb_iw1 <- read.csv("data/inputs/athena_3360_2096_220713.csv",na.strings=c(""," ","NA"))         # no pra,pvrt,cpw,ddisc,rdisc,edisc
  fullcnb_iw2 <- read.csv("data/inputs/athena_254_360_220713.csv",na.strings=c(""," ","NA"))           # for cpw, still no pra,pvrt,ddisc,rdisc,edisc
  
  fcnb_iw1 <- fullcnb_iw1 %>% filter(test_sessions_v.battery %in% c("adaptive_v_battery_v1","adaptive_v_battery_v0","adaptive_v_battery_v2"))  # only 10 rows from adaptive and 15 rows for siteid == "adaptive_v"
  fcnb_iw2 <- fullcnb_iw2 %>% filter(test_sessions_v.battery %in% c("adaptive_v_battery_v1","adaptive_v_battery_v0","adaptive_v_battery_v2"))
  
}

# new IW dataset from 7/13/22
new_iw <- read.csv("data/inputs/cnb/athena_195_360.csv",na.strings=c(""," ","NA"))
new_iw_adaptive <- new_iw %>% mutate(bblid = as.numeric(test_sessions_v.bblid)) %>% 
  filter(test_sessions_v.battery %in% c("adaptive_v_battery_v1","adaptive_v_battery_v2","cat_adaptive_v_battery_link1","cat_adaptive_v_battery_link3"),
         bblid > 9999) %>% arrange(bblid) %>% 
  dplyr::select(test_sessions.datasetid:test_sessions.famid,bblid,test_sessions_v.admin_comments:test_sessions_v.batt_consent,test_sessions_v.deleted_flag:PRA_D.AGE)

new_iw_ad <- new_iw_adaptive %>% filter(test_sessions_v.battery %in% c("adaptive_v_battery_v1","adaptive_v_battery_v2"))     # 252 unique bblid, ADT, CPF, MEDF, PMAT, VSPLOT, SVOLT, SPCPTNL, AIM, GNG, CPW, DIGSYM, ER40, PVRT, PRA (but empty) no MPRACT, PCET, SLNB
# new_iw_cat_1 <- new_iw_adaptive %>% filter(test_sessions_v.battery == "cat_adaptive_v_battery_link1")      # 259 unique bblid, but no actual data
# new_iw_cat_3 <- new_iw_adaptive %>% filter(test_sessions_v.battery == "cat_adaptive_v_battery_link3")      # 261 unique bblid, but no actual data

# get rid of datasetid 49043, 49044, 47905, 47859, 48036, 48505 for duplicate bblid, look into bblid 104265
new_iw_ad <- new_iw_ad %>% filter(test_sessions_v.battery_complete == 1, is.na(test_sessions_v.deleted_flag))
# bblid 23402 has two sessions: take first session (datasetid 52273) and replace EDISC (doesn't exist in this CSV),DIGSYM,GNG150 data with second session (52277) & remove second session
new_iw_ad[which(new_iw_ad$test_sessions.datasetid == 52273),c(4339:5102,5342:7446)] <- new_iw_ad[which(new_iw_ad$test_sessions.datasetid == 52277),c(4339:5102,5342:7446)]
new_iw_ad <- new_iw_ad %>% filter(test_sessions.datasetid != 52277)

# make new demos variable for new_iw
demo_new_iw <- new_iw_ad %>% dplyr::select(test_sessions.datasetid:test_sessions_v.endtime)  # 269 rows with unique bblid as of 10/3/22  

# old stuff
{
  demo_from_iw1 <- fullcnb_iw1 %>% dplyr::select(test_sessions.datasetid:test_sessions_v.endtime)
  demo_from_iw2 <- fullcnb_iw2 %>% dplyr::select(test_sessions.datasetid:test_sessions_v.endtime)
}

# for old_iw
{
  fullcnb_iw <- old_iw
  demo_from_iw <- old_iw %>% dplyr::select(BBLID:test_sessions_v.endtime)
  
  # separate out iw data into each respective test for SMVE measure calculation
  ADT_iw <- fullcnb_iw %>% dplyr::select(matches("^ADT36_A.ADT36A_CR$|^ADT36_A.ADT36A_PC$|^ADT36_A.ADT36A_RTCR$|ADT36_A.ADT36A_QID")) %>% cbind(demo_from_iw,.)
  CPF_iw <- fullcnb_iw %>% dplyr::select(matches("^CPF_B.CPF_CR$|^CPF_B.CPF_RTCR$|^CPF_B.CPF_W_RTCR$|CPF_B.CPF_TRIAL")) %>% cbind(demo_from_iw,.)
  CPW_iw <- fullcnb_iw %>% dplyr::select(matches("^CPW_A.CPW_CR$|^CPW_A.CPW_RTCR$|^CPW_A.CPW_W_RTCR$|CPW_A.CPW_TRIAL")) %>% cbind(demo_from_iw,.)
  DDISC_iw <- f_cnb %>% dplyr::select(matches("DDISC")) %>% cbind(demo_f_cnb,.)
  EDISC_iw <- f_cnb %>% dplyr::select(matches("EDISC")) %>% cbind(demo_f_cnb,.)
  # ER40_iw
  MEDF_iw <- fullcnb_iw %>% dplyr::select(matches("^MEDF36_A.MEDF36A_CR$|^MEDF36_A.MEDF36A_PC$|^MEDF36_A.MEDF36A_RTCR$|MEDF36_A.MEDF36A_QID")) %>% cbind(demo_from_iw,.)
  PMAT_iw <- fullcnb_iw %>% dplyr::select(matches("^PMAT24_A.PMAT24_A_CR$|^PMAT24_A.PMAT24_A_PC$|^PMAT24_A.PMAT24_A_RTCR$|PMAT24_A.PMAT24_A_QID0000")) %>% cbind(demo_from_iw,.)
  PLOT_iw <- fullcnb_iw %>% dplyr::select(matches("VSPLOT15.VSPLOT15_"))
  PLOT_iw <- PLOT_iw %>% dplyr::select(matches("_CR$|_PC$|RTCR$|_DEG_OFF_|_RT_|_ANS_|_EXCESS_|_CORR")) %>% cbind(demo_from_iw,.)
  RDISC_iw <- f_cnb %>% dplyr::select(matches("RDISC")) %>% cbind(demo_f_cnb,.)
  VOLT_iw <- fullcnb_iw %>% dplyr::select(matches("^SVOLT_A.SVOLT_CR$|^SVOLT_A.SVOLT_RTCR$|^SVOLT_A.SVOLT_W_RTCR$|SVOLT_A.SVOLT_TRIAL")) %>% cbind(demo_from_iw,.)
  
  # missing PRA, PVRT, ER40? in old_iw
  
  # rapid tests
  CPT_iw <- fullcnb_iw %>% dplyr::select(matches("^SPCPTNL.SCPL_TP$|^SPCPTNL.SCPL_TPRT$|^SPCPTNL.SCPT_QID")) %>% cbind(demo_from_iw,.)
  AIM_iw <- fullcnb_iw %>% dplyr::select(matches("^AIM.AIMTOT$|^AIM.AIMTOTRT$|^AIM.AIM_QID")) %>% cbind(demo_from_iw,.)
  GNG_iw <- fullcnb_iw %>% dplyr::select(matches("^GNG150.GNG150_CR$|^GNG150.GNG150_RTCR$|^GNG150.GNG150_QID")) %>% cbind(demo_from_iw,.)
  DIGSYM_iw <- fullcnb_iw %>% dplyr::select(matches("^DIGSYM.DSCOR$|^DIGSYM.DSCORRT$|^DIGSYM.DS_TP$|^DIGSYM.DS_TPRT$|^DIGSYM.DSMEMCR$|^DIGSYM.DSMCRRT$|DIGSYM.DS_MEM_QID|DIGSYM.DS_QID")) %>% cbind(demo_from_iw,.)
  DIGSYM_iw <- DIGSYM_iw[,1:559] # from 194 x 2124 --> 559
  DIGSYM_iw <- DIGSYM_iw[,1:469]
}

# separate itemwise data into each specific test
{
  dat <- new_iw_ad
  demos <- demo_new_iw
  
  ADT_iw <- dat %>% dplyr::select(matches("^ADT36_A.ADT36A_CR$|^ADT36_A.ADT36A_PC$|^ADT36_A.ADT36A_RTCR$|ADT36_A.ADT36A_QID")) %>% cbind(demos,.)   # 252 rows, 8/3/22
  CPF_iw <- dat %>% dplyr::select(matches("^CPF_B.CPF_CR$|^CPF_B.CPF_RTCR$|^CPF_B.CPF_W_RTCR$|CPF_B.CPF_TRIAL")) %>% cbind(demos,.)
  CPW_iw <- dat %>% dplyr::select(matches("^CPW_A.CPW_CR$|^CPW_A.CPW_RTCR$|^CPW_A.CPW_W_RTCR$|CPW_A.CPW_TRIAL")) %>% cbind(demos,.)
  DDISC_iw <- f_cnb %>% dplyr::select(matches("DDISC")) %>% cbind(demo_f_cnb,.) %>%  # duplicates exist for 23402 (delete datasetid 52277) and 22012 (delete datasetid 47859 since missing data)
    filter(test_sessions.datasetid %notin% c(52277,47859))   # 279 rows, 10/3/22
  EDISC_iw <- f_cnb %>% dplyr::select(matches("EDISC")) %>% cbind(demo_f_cnb,.) %>% dplyr::select(datasetid_platform:EDISC.valid_code,EDISC.q_101_resp:EDISC.test) %>%  # duplicates exist for 23402 (delete datasetid 52273 since EDISC crashed) and 22012 (delete datasetid 47859 since missing data)
    filter(test_sessions.datasetid %notin% c(52273,47859))
  ER40_iw <- dat %>% dplyr::select(matches("^ER40_D.ER40D_CR$|^ER40_D.ER40D_RTCR$|ER40_D.ER40D_QID")) %>% cbind(demos,.)
  MEDF_iw <- dat %>% dplyr::select(matches("^MEDF36_A.MEDF36A_CR$|^MEDF36_A.MEDF36A_PC$|^MEDF36_A.MEDF36A_RTCR$|MEDF36_A.MEDF36A_QID")) %>% cbind(demos,.)
  PMAT_iw <- dat %>% dplyr::select(matches("^PMAT24_A.PMAT24_A_CR$|^PMAT24_A.PMAT24_A_PC$|^PMAT24_A.PMAT24_A_RTCR$|PMAT24_A.PMAT24_A_QID0000")) %>% cbind(demos,.)
  PLOT_iw <- dat %>% dplyr::select(matches("VSPLOT15.VSPLOT15_"))
  PLOT_iw <- PLOT_iw %>% dplyr::select(matches("_CR$|_PC$|RTCR$|_DEG_OFF_|_RT_|_ANS_|_EXCESS_|_CORR")) %>% cbind(demos,.)
  PRA_iw <- new_iw %>% mutate(bblid = as.numeric(test_sessions_v.bblid)) %>% arrange(bblid) %>%   # 259 rows, 8/3/22
    filter(test_sessions.siteid == "adaptive_v",test_sessions_v.battery == "PRA_D", bblid > 9999) %>% 
    dplyr::select(matches("test_session|^bblid|PRA_D")) %>% dplyr::select(test_sessions.datasetid:test_sessions.famid,bblid,test_sessions_v.admin_comments:test_sessions_v.batt_consent,test_sessions_v.deleted_flag:PRA_D.AGE)
  PVRT_iw <- dat %>% dplyr::select(matches("^SPVRT_A.SPVRTA_CR$|^SPVRT_A.SPVRTA_PC$|^SPVRT_A.SPVRTA_W_CR$|^SPVRT_A.SPVRTA_W_PC$|^SPVRT_A.SPVRTA_RTCR$|SPVRT_A.SPVRTA_QID")) %>% cbind(demos,.)
  RDISC_iw <- f_cnb %>% dplyr::select(matches("RDISC")) %>% cbind(demo_f_cnb,.) %>% dplyr::select(datasetid_platform:KRDISC.test) %>%  # duplicates exist for 23402 (delete datasetid 52277) and 22012 (delete datasetid 47859 since missing data)
    filter(test_sessions.datasetid %notin% c(52277,47859))
  VOLT_iw <- dat %>% dplyr::select(matches("^SVOLT_A.SVOLT_CR$|^SVOLT_A.SVOLT_RTCR$|^SVOLT_A.SVOLT_W_RTCR$|SVOLT_A.SVOLT_TRIAL")) %>% cbind(demos,.)
  
  # rapid tests
  AIM_iw <- dat %>% dplyr::select(matches("^AIM.AIMTOT$|^AIM.AIMTOTRT$|^AIM.AIM_QID")) %>% cbind(demos,.)
  CPT_iw <- dat %>% dplyr::select(matches("^SPCPTNL.SCPL_TP$|^SPCPTNL.SCPL_TPRT$|^SPCPTNL.SCPT_QID")) %>% cbind(demos,.)
  GNG_iw <- dat %>% dplyr::select(matches("^GNG150.GNG150_CR$|^GNG150.GNG150_RTCR$|^GNG150.GNG150_QID")) %>% cbind(demos,.)
  DIGSYM_iw <- dat %>% dplyr::select(matches("^DIGSYM.DSCOR$|^DIGSYM.DSCORRT$|^DIGSYM.DS_TP$|^DIGSYM.DS_TPRT$|^DIGSYM.DSMEMCR$|^DIGSYM.DSMCRRT$|DIGSYM.DS_MEM_QID|DIGSYM.DS_QID")) %>% cbind(demos,.)
  # last real response in DIGSYM_iw is DIGSYM.DS_QID000170_TTR
  
  # temp <- DIGSYM_iw %>% dplyr::select(!matches("_QID$"))
  # temp$not_NA <- rowSums(!is.na(temp[,23:ncol(temp)]))
  # head(sort(temp$not_NA,decreasing = T),15)
  
  DIGSYM_iw <- DIGSYM_iw %>% dplyr::select(test_sessions.datasetid:DIGSYM.DS_QID000170_TTR) # still true, 8/3/22
  
  
  # all memory + ER40 + MEDF columns to make separate full CNB target/foil scores
  {
    
    # ER40 separation
    er40_all <- dat %>% dplyr::select(matches("ER40")) %>% cbind(demos,.)   # emotive vs neutral
    er40_emo <- er40_all %>% dplyr::select(matches("_ANG|_FEAR|_HAP|_SAD")) %>% cbind(demos,.) %>% 
      mutate(ER40_D.ER40D_EMO = rowSums(.[,15:18]))
    # use iw data for ER40_D.ER40D_EMORTCR
    # code below was used to check that the item ordering that I got from Lucky matched preexisting data. it did! yay:)
    er40_ang <- er40_all %>% dplyr::select(ER40_D.ER40D_QID000001_RESP:ER40_D.ER40D_QID000004_CORR,ER40_D.ER40D_QID000021_RESP:ER40_D.ER40D_QID000024_CORR)    # %>% mutate(ang = rowSums(.[,c(5,10,15,20,25,30,35,40)]))
    er40_fear <- er40_all %>% dplyr::select(ER40_D.ER40D_QID000005_RESP:ER40_D.ER40D_QID000008_CORR,ER40_D.ER40D_QID000025_RESP:ER40_D.ER40D_QID000028_CORR)   #  %>% mutate(fear = rowSums(.[,c(5,10,15,20,25,30,35,40)]))
    er40_hap <- er40_all %>% dplyr::select(ER40_D.ER40D_QID000009_RESP:ER40_D.ER40D_QID000012_CORR,ER40_D.ER40D_QID000029_RESP:ER40_D.ER40D_QID000032_CORR)    # %>% mutate(hap = rowSums(.[,c(5,10,15,20,25,30,35,40)]))
    er40_sad <- er40_all %>% dplyr::select(ER40_D.ER40D_QID000017_RESP:ER40_D.ER40D_QID000020_CORR,ER40_D.ER40D_QID000037_RESP:ER40_D.ER40D_QID000040_CORR)    # %>% mutate(sad = rowSums(.[,c(5,10,15,20,25,30,35,40)]))
    
    er40_noe <- er40_all %>% dplyr::select(ER40_D.ER40D_QID000013_RESP:ER40_D.ER40D_QID000016_CORR,ER40_D.ER40D_QID000033_RESP:ER40_D.ER40D_QID000036_CORR)    # %>% mutate(noe = rowSums(.[,c(5,10,15,20,25,30,35,40)]))
    
    er40_resp <- er40_all %>% dplyr::select(matches("_TTR"))
    er40_resp$ER40_D.ER40D_EMORTCR <- rowMedians(as.matrix(er40_resp %>% dplyr::select(ER40_D.ER40D_QID000001_TTR:ER40_D.ER40D_QID000012_TTR,
                                                                                       ER40_D.ER40D_QID000017_TTR:ER40_D.ER40D_QID000032_TTR,
                                                                                       ER40_D.ER40D_QID000037_TTR:ER40_D.ER40D_QID000040_TTR)))
    ER40_EMO_iw <- cbind(demos,er40_ang,er40_fear,er40_hap,er40_sad,er40_emo$ER40_D.ER40D_EMO,er40_resp$ER40_D.ER40D_EMORTCR)
    ER40_NEU_iw <- cbind(demos,er40_noe,er40_all %>% dplyr::select(matches("_NOE")))
    
    
    # MEDDF separation
    medf_all <- dat %>% dplyr::select(matches("MEDF")) %>% cbind(demos,.)   # same vs different
    # medf_emo <- medf_all %>% dplyr::select(matches("_ANG|_FEAR|_HAP|_SAD")) %>% cbind(demos,.) # initially used this to calculate MEDF36_A.MEDF36A_DIF, but realized this doesn't work 
    # use iw data for MEDF36_A.MEDF36A_SAMERTCR
    medf_resp <- medf_all %>% dplyr::select(matches("_TTR"))
    medf_resp$MEDF36_A.MEDF36A_DIFRT <- rowMedians(as.matrix(medf_resp %>% dplyr::select(MEDF36_A.MEDF36A_QID000001_TTR:MEDF36_A.MEDF36A_QID000008_TTR,
                                                                                         MEDF36_A.MEDF36A_QID000010_TTR:MEDF36_A.MEDF36A_QID000012_TTR,
                                                                                         MEDF36_A.MEDF36A_QID000014_TTR:MEDF36_A.MEDF36A_QID000020_TTR,
                                                                                         MEDF36_A.MEDF36A_QID000022_TTR:MEDF36_A.MEDF36A_QID000035_TTR)))
    
    medf_corr <- medf_all %>% dplyr::select(matches("_CORR"))
    medf_corr$MEDF36_A.MEDF36A_DIF <- rowSums(medf_corr %>% dplyr::select(MEDF36_A.MEDF36A_QID000001_CORR:MEDF36_A.MEDF36A_QID000008_CORR,
                                                                          MEDF36_A.MEDF36A_QID000010_CORR:MEDF36_A.MEDF36A_QID000012_CORR,
                                                                          MEDF36_A.MEDF36A_QID000014_CORR:MEDF36_A.MEDF36A_QID000020_CORR,
                                                                          MEDF36_A.MEDF36A_QID000022_CORR:MEDF36_A.MEDF36A_QID000035_CORR))
    
    
    MEDF_DIF_iw <- cbind(demos,medf_all %>% dplyr::select(MEDF36_A.MEDF36A_QID000001_RESP:MEDF36_A.MEDF36A_QID000008_CORR,
                                                          MEDF36_A.MEDF36A_QID000010_RESP:MEDF36_A.MEDF36A_QID000012_CORR,
                                                          MEDF36_A.MEDF36A_QID000014_RESP:MEDF36_A.MEDF36A_QID000020_CORR,
                                                          MEDF36_A.MEDF36A_QID000022_RESP:MEDF36_A.MEDF36A_QID000035_CORR),
                         medf_corr$MEDF36_A.MEDF36A_DIF,medf_resp$MEDF36_A.MEDF36A_DIFRT)
    
    # recalculating RTCR for same items
    medf_same_resp <- medf_all %>% dplyr::select(MEDF36_A.MEDF36A_QID000009_RESP:MEDF36_A.MEDF36A_QID000009_CORR,
                                                 MEDF36_A.MEDF36A_QID000013_RESP:MEDF36_A.MEDF36A_QID000013_CORR,
                                                 MEDF36_A.MEDF36A_QID000021_RESP:MEDF36_A.MEDF36A_QID000021_CORR,
                                                 MEDF36_A.MEDF36A_QID000036_RESP:MEDF36_A.MEDF36A_QID000036_CORR)
    medf_same_resp$MEDF36_A.MEDF36A_SAME_RTCR <- rowMedians(as.matrix(medf_same_resp %>% dplyr::select(matches("_TTR"))))
    
    MEDF_SAME_iw <- cbind(demos, medf_same_resp %>% dplyr::select(MEDF36_A.MEDF36A_QID000009_RESP:MEDF36_A.MEDF36A_QID000036_CORR),
                          medf_all %>% dplyr::select(MEDF36_A.MEDF36A_SAME_CR),medf_same_resp %>% dplyr::select(MEDF36_A.MEDF36A_SAME_RTCR))
    
    # CPF separation
    cpf_all <- dat %>% dplyr::select(matches("CPF")) %>% cbind(demos,.)     # targets vs foils (TP vs TN)
    
    CPF_targets <- cpf_all %>% dplyr::select(CPF_B.CPF_TRIAL000003_RESP:CPF_B.CPF_TRIAL000006_CORR,
                                             CPF_B.CPF_TRIAL000008_RESP:CPF_B.CPF_TRIAL000010_CORR,
                                             CPF_B.CPF_TRIAL000014_RESP:CPF_B.CPF_TRIAL000014_CORR,
                                             CPF_B.CPF_TRIAL000019_RESP:CPF_B.CPF_TRIAL000022_CORR,
                                             CPF_B.CPF_TRIAL000025_RESP:CPF_B.CPF_TRIAL000025_CORR,
                                             CPF_B.CPF_TRIAL000028_RESP:CPF_B.CPF_TRIAL000029_CORR,
                                             CPF_B.CPF_TRIAL000031_RESP:CPF_B.CPF_TRIAL000031_CORR,
                                             CPF_B.CPF_TRIAL000034_RESP:CPF_B.CPF_TRIAL000036_CORR,
                                             CPF_B.CPF_TRIAL000040_RESP:CPF_B.CPF_TRIAL000040_CORR,
                                             CPF_B.CPF_TP,CPF_B.CPF_TPRT) %>% cbind(demos,.)  # checked 8/3/22
    
    CPF_foils <- cpf_all %>% dplyr::select(CPF_B.CPF_TRIAL000001_RESP:CPF_B.CPF_TRIAL000002_CORR,
                                           CPF_B.CPF_TRIAL000007_RESP:CPF_B.CPF_TRIAL000007_CORR,
                                           CPF_B.CPF_TRIAL000011_RESP:CPF_B.CPF_TRIAL000013_CORR,
                                           CPF_B.CPF_TRIAL000015_RESP:CPF_B.CPF_TRIAL000018_CORR,
                                           CPF_B.CPF_TRIAL000023_RESP:CPF_B.CPF_TRIAL000024_CORR,
                                           CPF_B.CPF_TRIAL000026_RESP:CPF_B.CPF_TRIAL000027_CORR,
                                           CPF_B.CPF_TRIAL000030_RESP:CPF_B.CPF_TRIAL000030_CORR,
                                           CPF_B.CPF_TRIAL000032_RESP:CPF_B.CPF_TRIAL000033_CORR,
                                           CPF_B.CPF_TRIAL000037_RESP:CPF_B.CPF_TRIAL000039_CORR,
                                           CPF_B.CPF_TN,CPF_B.CPF_TNRT) %>% cbind(demos,.)   # checked 8/3/22
    
    # temporary to check TP and TPRT, TN and TNR
    { # TPRT are a little bit off but TP matches exactly, so i think we're safe
      temp_resp <- CPF_targets %>% dplyr::select(matches("TTR"))
      temp_resp$CPF_B.CPF_Tscore <- rowMedians(as.matrix(temp_resp)) 
      temp_resp <- cbind(temp_resp,CPF_targets$CPF_B.CPF_TPRT)
      
      temp_corr <- CPF_targets %>% dplyr::select(matches("CORR"))
      temp_corr$CPF_B.CPF_Tscore <- rowSums(temp_corr) 
      temp_corr <- cbind(temp_corr,CPF_targets$CPF_B.CPF_TP)
      
      # TNRT are a little bit off but TN matches exactly, so i think we're safe
      temp_resp <- CPF_foils %>% dplyr::select(matches("TTR"))
      temp_resp$CPF_B.CPF_Tscore <- rowMedians(as.matrix(temp_resp)) 
      temp_resp <- cbind(temp_resp,CPF_foils$CPF_B.CPF_TNRT)
      
      temp_corr <- CPF_foils %>% dplyr::select(matches("CORR"))
      temp_corr$CPF_B.CPF_Tscore <- rowSums(temp_corr) 
      temp_corr <- cbind(temp_corr,CPF_foils$CPF_B.CPF_TN)
                                           }
    
    cpw_all <- dat %>% dplyr::select(matches("CPW")) %>% cbind(demos,.)     # targets vs foils (TP vs TN)
    
    CPW_targets <- cpw_all %>% dplyr::select(CPW_A.CPW_TRIAL000001_RESP:CPW_A.CPW_TRIAL000001_CORR,
                                             CPW_A.CPW_TRIAL000003_RESP:CPW_A.CPW_TRIAL000003_CORR,
                                             CPW_A.CPW_TRIAL000005_RESP:CPW_A.CPW_TRIAL000005_CORR,
                                             CPW_A.CPW_TRIAL000009_RESP:CPW_A.CPW_TRIAL000010_CORR,
                                             CPW_A.CPW_TRIAL000012_RESP:CPW_A.CPW_TRIAL000013_CORR,
                                             CPW_A.CPW_TRIAL000016_RESP:CPW_A.CPW_TRIAL000017_CORR,
                                             CPW_A.CPW_TRIAL000019_RESP:CPW_A.CPW_TRIAL000019_CORR,
                                             CPW_A.CPW_TRIAL000021_RESP:CPW_A.CPW_TRIAL000022_CORR,
                                             CPW_A.CPW_TRIAL000024_RESP:CPW_A.CPW_TRIAL000025_CORR,
                                             CPW_A.CPW_TRIAL000029_RESP:CPW_A.CPW_TRIAL000029_CORR,
                                             CPW_A.CPW_TRIAL000031_RESP:CPW_A.CPW_TRIAL000032_CORR,
                                             CPW_A.CPW_TRIAL000036_RESP:CPW_A.CPW_TRIAL000036_CORR,
                                             CPW_A.CPW_TRIAL000039_RESP:CPW_A.CPW_TRIAL000040_CORR,
                                             CPW_A.CPW_TP,CPW_A.CPW_TPRT) %>% cbind(demos,.)  # checked 8/3/22
    
    CPW_foils <- cpw_all %>% dplyr::select(CPW_A.CPW_TRIAL000002_RESP:CPW_A.CPW_TRIAL000002_CORR,
                                           CPW_A.CPW_TRIAL000004_RESP:CPW_A.CPW_TRIAL000004_CORR,
                                           CPW_A.CPW_TRIAL000006_RESP:CPW_A.CPW_TRIAL000008_CORR,
                                           CPW_A.CPW_TRIAL000011_RESP:CPW_A.CPW_TRIAL000011_CORR,
                                           CPW_A.CPW_TRIAL000014_RESP:CPW_A.CPW_TRIAL000015_CORR,
                                           CPW_A.CPW_TRIAL000018_RESP:CPW_A.CPW_TRIAL000018_CORR,
                                           CPW_A.CPW_TRIAL000020_RESP:CPW_A.CPW_TRIAL000020_CORR,
                                           CPW_A.CPW_TRIAL000023_RESP:CPW_A.CPW_TRIAL000023_CORR,
                                           CPW_A.CPW_TRIAL000026_RESP:CPW_A.CPW_TRIAL000028_CORR,
                                           CPW_A.CPW_TRIAL000030_RESP:CPW_A.CPW_TRIAL000030_CORR,
                                           CPW_A.CPW_TRIAL000033_RESP:CPW_A.CPW_TRIAL000035_CORR,
                                           CPW_A.CPW_TRIAL000037_RESP:CPW_A.CPW_TRIAL000038_CORR,
                                           CPW_A.CPW_TN,CPW_A.CPW_TNRT) %>% cbind(demos,.)   # checked 8/3/22
    
    # temporary to check TP and TPRT, TN and TNR
    { # TPRT are a little bit off but TP matches exactly, so i think we're safe
      temp_resp <- CPW_targets %>% dplyr::select(matches("TTR"))
      temp_resp$CPW_A.CPW_Tscore <- rowMedians(as.matrix(temp_resp)) 
      temp_resp <- cbind(temp_resp,CPW_targets$CPW_A.CPW_TPRT)
      
      temp_corr <- CPW_targets %>% dplyr::select(matches("CORR"))
      temp_corr$CPW_A.CPW_Tscore <- rowSums(temp_corr) 
      temp_corr <- cbind(temp_corr,CPW_targets$CPW_A.CPW_TP)
      
      # TNRT are a little bit off but TN matches exactly, so i think we're safe
      temp_resp <- CPW_foils %>% dplyr::select(matches("TTR"))
      temp_resp$CPW_A.CPW_Tscore <- rowMedians(as.matrix(temp_resp)) 
      temp_resp <- cbind(temp_resp,CPW_foils$CPW_A.CPW_TNRT)
      
      temp_corr <- CPW_foils %>% dplyr::select(matches("CORR"))
      temp_corr$CPW_A.CPW_Tscore <- rowSums(temp_corr) 
      temp_corr <- cbind(temp_corr,CPW_foils$CPW_A.CPW_TN)
                                           }
    
    volt_all <- dat %>% dplyr::select(matches("VOLT")) %>% cbind(demos,.)   # targets vs foils (TP vs TN)
    
    VOLT_targets <- volt_all %>% dplyr::select(SVOLT_A.SVOLT_TRIAL000001_RESP:SVOLT_A.SVOLT_TRIAL000002_CORR,
                                               SVOLT_A.SVOLT_TRIAL000004_RESP:SVOLT_A.SVOLT_TRIAL000005_CORR,
                                               SVOLT_A.SVOLT_TRIAL000008_RESP:SVOLT_A.SVOLT_TRIAL000010_CORR,
                                               SVOLT_A.SVOLT_TRIAL000013_RESP:SVOLT_A.SVOLT_TRIAL000013_CORR,
                                               SVOLT_A.SVOLT_TRIAL000015_RESP:SVOLT_A.SVOLT_TRIAL000015_CORR,
                                               SVOLT_A.SVOLT_TRIAL000019_RESP:SVOLT_A.SVOLT_TRIAL000019_CORR,
                                               SVOLT_A.SVOLT_TP,SVOLT_A.SVOLT_TPRT) %>% cbind(demos,.)  # checked 8/3/22
    
    VOLT_foils <- volt_all %>% dplyr::select(SVOLT_A.SVOLT_TRIAL000003_RESP:SVOLT_A.SVOLT_TRIAL000003_CORR,
                                             SVOLT_A.SVOLT_TRIAL000006_RESP:SVOLT_A.SVOLT_TRIAL000007_CORR,
                                             SVOLT_A.SVOLT_TRIAL000011_RESP:SVOLT_A.SVOLT_TRIAL000012_CORR,
                                             SVOLT_A.SVOLT_TRIAL000014_RESP:SVOLT_A.SVOLT_TRIAL000014_CORR,
                                             SVOLT_A.SVOLT_TRIAL000016_RESP:SVOLT_A.SVOLT_TRIAL000018_CORR,
                                             SVOLT_A.SVOLT_TRIAL000020_RESP:SVOLT_A.SVOLT_TRIAL000020_CORR,
                                             SVOLT_A.SVOLT_TN,SVOLT_A.SVOLT_TNRT) %>% cbind(demos,.)   # checked 8/3/22
    
    # temporary to check TP and TPRT, TN and TNR
    { # TPRT are a little bit off but TP matches exactly, so i think we're safe
      temp_resp <- VOLT_targets %>% dplyr::select(matches("TTR"))
      temp_resp$SVOLT_A.SVOLT_Tscore <- rowMedians(as.matrix(temp_resp)) 
      temp_resp <- cbind(temp_resp,VOLT_targets$SVOLT_A.SVOLT_TPRT)
      
      temp_corr <- VOLT_targets %>% dplyr::select(matches("CORR"))
      temp_corr$SVOLT_A.SVOLT_Tscore <- rowSums(temp_corr) 
      temp_corr <- cbind(temp_corr,VOLT_targets$SVOLT_A.SVOLT_TP)
      
      # TNRT are a little bit off but TN matches exactly, so i think we're safe
      temp_resp <- VOLT_foils %>% dplyr::select(matches("TTR"))
      temp_resp$SVOLT_A.SVOLT_Tscore <- rowMedians(as.matrix(temp_resp)) 
      temp_resp <- cbind(temp_resp,VOLT_foils$SVOLT_A.SVOLT_TNRT)
      
      temp_corr <- VOLT_foils %>% dplyr::select(matches("CORR"))
      temp_corr$SVOLT_A.SVOLT_Tscore <- rowSums(temp_corr) 
      temp_corr <- cbind(temp_corr,VOLT_foils$SVOLT_A.SVOLT_TN)
                                             }
  }
  
  # adding summary scores to DISC tasks
  # disc_tasks <- all_cnb %>% dplyr::select(bblid,ddisc_sum:edisc_mcr) # for some reason, there are some scores missing when using this method
  # 
  # DDISC_iw <- left_join(DDISC_iw,disc_tasks %>% dplyr::select(bblid,ddisc_sum,ddisc_mcr),by=c("BBLID" = "bblid"))
  # EDISC_iw <- left_join(EDISC_iw,disc_tasks %>% dplyr::select(bblid,edisc_sum,edisc_mcr),by=c("BBLID" = "bblid"))
  # RDISC_iw <- left_join(RDISC_iw,disc_tasks %>% dplyr::select(bblid,rdisc_sum,rdisc_mcr),by=c("BBLID" = "bblid"))
  
  DDISC_iw$ddisc_sum <- rowSums((DDISC_iw %>% dplyr::select(matches("KDDISC.q")))-1)
  DDISC_iw$ddisc_mcr <- rowMedians(as.matrix(DDISC_iw %>% dplyr::select(matches("KDDISC.trr"))))
  
  RDISC_iw$rdisc_sum <- rowSums((RDISC_iw %>% dplyr::select(matches("RDISC.q")))-1)
  RDISC_iw$rdisc_mcr <- rowMedians(as.matrix(RDISC_iw %>% dplyr::select(matches("RDISC.trr"))))
  
  EDISC_iw$edisc_sum <- rowSums((EDISC_iw %>% dplyr::select(matches("_resp")))-1)
  EDISC_iw$edisc_mcr <- rowMedians(as.matrix(EDISC_iw %>% dplyr::select(matches("_ttr"))))
}






# QC data ----

# * SMVE ----
# SMVE for non-rapid tests

# creating performance validity scores for full CNB scores to:
#   1) show distributions of performance validity scores
#   2) keep only the best quality data
#   3) compare this with other validity measures (manual + autovalidation rules)

{
  # ADT
  dat <-ADT_iw[,c(grep("_CORR",colnames(ADT_iw)),grep("_TTR",colnames(ADT_iw)))]
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
  ADT_iw <- data.frame(ADT_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(ADT_iw)[(ncol(ADT_iw)-2):ncol(ADT_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0% or 100% so no adjustment needed
  
  
  # AIM
  
  # different method, look at distribution of scores
  
  
  # CPF
  dat <-CPF_iw[,c(grep("_CORR",colnames(CPF_iw)),grep("_TTR",colnames(CPF_iw)))]
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
  CPF_iw <- data.frame(CPF_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(CPF_iw)[(ncol(CPF_iw)-2):ncol(CPF_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  CPF_iw$PFscore1 <- ifelse(CPF_iw$CPF_B.CPF_CR == 40, max(CPF_iw$PFscore1,na.rm = T),CPF_iw$PFscore1)
  CPF_iw$PFscore2 <- ifelse(CPF_iw$CPF_B.CPF_CR == 40, max(CPF_iw$PFscore2,na.rm = T),CPF_iw$PFscore2)
  CPF_iw$SMVE <- (0.42 * CPF_iw$outlier_score_2cut) + (0.02 * CPF_iw$acc3e) + (0.05 * CPF_iw$PFscore1) + (0.50 * CPF_iw$PFscore2)
  
  
  # CPF targets
  dat <-CPF_targets[,c(grep("_CORR",colnames(CPF_targets)),grep("_TTR",colnames(CPF_targets)))]
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
  CPF_targets <- data.frame(CPF_targets,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(CPF_targets)[(ncol(CPF_targets)-2):ncol(CPF_targets)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  CPF_targets$PFscore1 <- ifelse(CPF_targets$CPF_B.CPF_TP == 20, max(CPF_targets$PFscore1,na.rm = T),CPF_targets$PFscore1)
  CPF_targets$PFscore2 <- ifelse(CPF_targets$CPF_B.CPF_TP == 20, max(CPF_targets$PFscore2,na.rm = T),CPF_targets$PFscore2)
  CPF_targets$SMVE <- (0.42 * CPF_targets$outlier_score_2cut) + (0.02 * CPF_targets$acc3e) + (0.05 * CPF_targets$PFscore1) + (0.50 * CPF_targets$PFscore2)
  
  
  # CPF foils
  dat <-CPF_foils[,c(grep("_CORR",colnames(CPF_foils)),grep("_TTR",colnames(CPF_foils)))]
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
  CPF_foils <- data.frame(CPF_foils,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(CPF_foils)[(ncol(CPF_foils)-2):ncol(CPF_foils)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  CPF_foils$PFscore1 <- ifelse(CPF_foils$CPF_B.CPF_TN == 20, max(CPF_foils$PFscore1,na.rm = T),CPF_foils$PFscore1)
  CPF_foils$PFscore2 <- ifelse(CPF_foils$CPF_B.CPF_TN == 20, max(CPF_foils$PFscore2,na.rm = T),CPF_foils$PFscore2)
  CPF_foils$SMVE <- (0.42 * CPF_foils$outlier_score_2cut) + (0.02 * CPF_foils$acc3e) + (0.05 * CPF_foils$PFscore1) + (0.50 * CPF_foils$PFscore2)
  
  
  # CPT
  
  # different method, look at distribution of scores
  
  
  # CPW
  dat <-CPW_iw[,c(grep("_CORR",colnames(CPW_iw)),grep("_TTR",colnames(CPW_iw)))]
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
  CPW_iw <- data.frame(CPW_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(CPW_iw)[(ncol(CPW_iw)-2):ncol(CPW_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  CPW_iw$PFscore1 <- ifelse(CPW_iw$CPW_A.CPW_CR == 40, max(CPW_iw$PFscore1,na.rm = T),CPW_iw$PFscore1)
  CPW_iw$PFscore2 <- ifelse(CPW_iw$CPW_A.CPW_CR == 40, max(CPW_iw$PFscore2,na.rm = T),CPW_iw$PFscore2)
  CPW_iw$SMVE <- (0.42 * CPW_iw$outlier_score_2cut) + (0.02 * CPW_iw$acc3e) + (0.05 * CPW_iw$PFscore1) + (0.50 * CPW_iw$PFscore2)
  
  
  # CPW targets
  dat <-CPW_targets[,c(grep("_CORR",colnames(CPW_targets)),grep("_TTR",colnames(CPW_targets)))]
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
  CPW_targets <- data.frame(CPW_targets,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(CPW_targets)[(ncol(CPW_targets)-2):ncol(CPW_targets)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  CPW_targets$PFscore1 <- ifelse(CPW_targets$CPW_A.CPW_TP == 20, max(CPW_targets$PFscore1,na.rm = T),CPW_targets$PFscore1)
  CPW_targets$PFscore2 <- ifelse(CPW_targets$CPW_A.CPW_TP == 20, max(CPW_targets$PFscore2,na.rm = T),CPW_targets$PFscore2)
  CPW_targets$SMVE <- (0.42 * CPW_targets$outlier_score_2cut) + (0.02 * CPW_targets$acc3e) + (0.05 * CPW_targets$PFscore1) + (0.50 * CPW_targets$PFscore2)
  
  
  # CPW foils
  dat <-CPW_foils[,c(grep("_CORR",colnames(CPW_foils)),grep("_TTR",colnames(CPW_foils)))]
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
  CPW_foils <- data.frame(CPW_foils,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(CPW_foils)[(ncol(CPW_foils)-2):ncol(CPW_foils)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  CPW_foils$PFscore1 <- ifelse(CPW_foils$CPW_A.CPW_TN == 20, max(CPW_foils$PFscore1,na.rm = T),CPW_foils$PFscore1)
  CPW_foils$PFscore2 <- ifelse(CPW_foils$CPW_A.CPW_TN == 20, max(CPW_foils$PFscore2,na.rm = T),CPW_foils$PFscore2)
  CPW_foils$SMVE <- (0.42 * CPW_foils$outlier_score_2cut) + (0.02 * CPW_foils$acc3e) + (0.05 * CPW_foils$PFscore1) + (0.50 * CPW_foils$PFscore2)
  
  
  # DDISC 
  dat <-DDISC_iw[,c(grep(".q_",colnames(DDISC_iw)),grep(".trr_",colnames(DDISC_iw)))]
  dat[,1:34] <- dat[,1:34]-1
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
  DDISC_iw <- data.frame(DDISC_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(DDISC_iw)[(ncol(DDISC_iw)-2):ncol(DDISC_iw)] <- c("PFscore1","PFscore2","SMVE")
  # 0% and 100% scores need adjusting
  
  DDISC_iw$PFscore1 <- ifelse(DDISC_iw$ddisc_sum == 34, max(DDISC_iw$PFscore1,na.rm = T),
                              ifelse(DDISC_iw$ddisc_sum == 0, min(DDISC_iw$PFscore1,na.rm = T),DDISC_iw$PFscore1))
  DDISC_iw$PFscore2 <- ifelse(DDISC_iw$ddisc_sum == 34, max(DDISC_iw$PFscore2,na.rm = T),
                              ifelse(DDISC_iw$ddisc_sum == 0, min(DDISC_iw$PFscore2,na.rm = T),DDISC_iw$PFscore2))
  DDISC_iw$SMVE <- (0.34 * DDISC_iw$outlier_score_2cut) + (0 * DDISC_iw$acc3e) + (0.22 * DDISC_iw$PFscore1) + (0.44 * DDISC_iw$PFscore2)
  # still 1 row with no SMVE because of missing itemwise
  
  
  # DIGSYM
  
  # different method, look at distribution of scores
  
  
  # EDISC 
  dat <-EDISC_iw[,c(grep("_resp",colnames(EDISC_iw)),grep("_ttr",colnames(EDISC_iw)))]
  dat[,1:34] <- dat[,1:34]-1
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
  EDISC_iw <- data.frame(EDISC_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(EDISC_iw)[(ncol(EDISC_iw)-2):ncol(EDISC_iw)] <- c("PFscore1","PFscore2","SMVE")
  # 0% and 100% scores need adjusting
  
  EDISC_iw$PFscore1 <- ifelse(EDISC_iw$edisc_sum == 34, max(EDISC_iw$PFscore1,na.rm = T),
                              ifelse(EDISC_iw$edisc_sum == 0, min(EDISC_iw$PFscore1,na.rm = T),EDISC_iw$PFscore1))
  EDISC_iw$PFscore2 <- ifelse(EDISC_iw$edisc_sum == 34, max(EDISC_iw$PFscore2,na.rm = T),
                              ifelse(EDISC_iw$edisc_sum == 0, min(EDISC_iw$PFscore2,na.rm = T),EDISC_iw$PFscore2))
  EDISC_iw$SMVE <- (0.34 * EDISC_iw$outlier_score_2cut) + (0 * EDISC_iw$acc3e) + (0.22 * EDISC_iw$PFscore1) + (0.44 * EDISC_iw$PFscore2)
  # still 2 rows with no SMVE because of missing itemwise
  
  
  # ER40
  dat <-ER40_iw[,c(grep("_CORR",colnames(ER40_iw)),grep("_TTR",colnames(ER40_iw)))]
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
  ER40_iw <- data.frame(ER40_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(ER40_iw)[(ncol(ER40_iw)-2):ncol(ER40_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  ER40_iw$PFscore1 <- ifelse(ER40_iw$ER40_D.ER40D_CR == 40, max(ER40_iw$PFscore1,na.rm = T),ER40_iw$PFscore1)
  ER40_iw$PFscore2 <- ifelse(ER40_iw$ER40_D.ER40D_CR == 40, max(ER40_iw$PFscore2,na.rm = T),ER40_iw$PFscore2)
  ER40_iw$SMVE <- (0.34 * ER40_iw$outlier_score_2cut) + (0 * ER40_iw$acc3e) + (0.22 * ER40_iw$PFscore1) + (0.44 * ER40_iw$PFscore2)
  
  
  # separated ER40, emotive
  dat <-ER40_EMO_iw[,c(grep("_CORR",colnames(ER40_EMO_iw)),grep("_TTR",colnames(ER40_EMO_iw)))]
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
  ER40_EMO_iw <- data.frame(ER40_EMO_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(ER40_EMO_iw)[(ncol(ER40_EMO_iw)-2):ncol(ER40_EMO_iw)] <- c("PFscore1","PFscore2","SMVE")
  # both 0% and 100% scores need adjusting
  
  ER40_EMO_iw$PFscore1 <- ifelse(ER40_EMO_iw$er40_emo.ER40_D.ER40D_EMO == 0, min(ER40_EMO_iw$PFscore1,na.rm = T),
                                 ifelse(ER40_EMO_iw$er40_emo.ER40_D.ER40D_EMO == 32, max(ER40_EMO_iw$PFscore1,na.rm = T),ER40_EMO_iw$PFscore1))
  ER40_EMO_iw$PFscore2 <- ifelse(ER40_EMO_iw$er40_emo.ER40_D.ER40D_EMO == 0, min(ER40_EMO_iw$PFscore2,na.rm = T),
                                 ifelse(ER40_EMO_iw$er40_emo.ER40_D.ER40D_EMO == 32, max(ER40_EMO_iw$PFscore2,na.rm = T),ER40_EMO_iw$PFscore2))
  ER40_EMO_iw$SMVE <- (0.34 * ER40_EMO_iw$outlier_score_2cut) + (0 * ER40_EMO_iw$acc3e) + (0.22 * ER40_EMO_iw$PFscore1) + (0.44 * ER40_EMO_iw$PFscore2)
  
  
  # separated ER40, neutral
  dat <-ER40_NEU_iw[,c(grep("_CORR",colnames(ER40_NEU_iw)),grep("_TTR",colnames(ER40_NEU_iw)))]
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
  ER40_NEU_iw <- data.frame(ER40_NEU_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(ER40_NEU_iw)[(ncol(ER40_NEU_iw)-2):ncol(ER40_NEU_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  ER40_NEU_iw$PFscore1 <- ifelse(ER40_NEU_iw$ER40_D.ER40D_NOE == 8, max(ER40_NEU_iw$PFscore1,na.rm = T),ER40_NEU_iw$PFscore1)
  ER40_NEU_iw$PFscore2 <- ifelse(ER40_NEU_iw$ER40_D.ER40D_NOE == 8, max(ER40_NEU_iw$PFscore2,na.rm = T),ER40_NEU_iw$PFscore2)
  ER40_NEU_iw$SMVE <- (0.34 * ER40_NEU_iw$outlier_score_2cut) + (0 * ER40_NEU_iw$acc3e) + (0.22 * ER40_NEU_iw$PFscore1) + (0.44 * ER40_NEU_iw$PFscore2)
  
  
  
  # GNG
  
  # different method, look at distribution of scores
  
  
  # MEDF
  dat <-MEDF_iw[,c(grep("_CORR",colnames(MEDF_iw)),grep("_TTR",colnames(MEDF_iw)))]
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
  MEDF_iw <- data.frame(MEDF_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(MEDF_iw)[(ncol(MEDF_iw)-2):ncol(MEDF_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  MEDF_iw$PFscore1 <- ifelse(MEDF_iw$MEDF36_A.MEDF36A_PC == 100, max(MEDF_iw$PFscore1,na.rm = T),MEDF_iw$PFscore1)
  MEDF_iw$PFscore2 <- ifelse(MEDF_iw$MEDF36_A.MEDF36A_PC == 100, max(MEDF_iw$PFscore2,na.rm = T),MEDF_iw$PFscore2)
  MEDF_iw$SMVE <- (0.34 * MEDF_iw$outlier_score_2cut) + (0 * MEDF_iw$acc3e) + (0.22 * MEDF_iw$PFscore1) + (0.44 * MEDF_iw$PFscore2)
  
  # separated MEDF, different
  dat <-MEDF_DIF_iw[,c(grep("_CORR",colnames(MEDF_DIF_iw)),grep("_TTR",colnames(MEDF_DIF_iw)))]
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
  MEDF_DIF_iw <- data.frame(MEDF_DIF_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(MEDF_DIF_iw)[(ncol(MEDF_DIF_iw)-2):ncol(MEDF_DIF_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  MEDF_DIF_iw$PFscore1 <- ifelse(MEDF_DIF_iw$medf_corr.MEDF36_A.MEDF36A_DIF == 32, max(MEDF_DIF_iw$PFscore1,na.rm = T),MEDF_DIF_iw$PFscore1)
  MEDF_DIF_iw$PFscore2 <- ifelse(MEDF_DIF_iw$medf_corr.MEDF36_A.MEDF36A_DIF == 32, max(MEDF_DIF_iw$PFscore2,na.rm = T),MEDF_DIF_iw$PFscore2)
  MEDF_DIF_iw$SMVE <- (0.34 * MEDF_DIF_iw$outlier_score_2cut) + (0 * MEDF_DIF_iw$acc3e) + (0.22 * MEDF_DIF_iw$PFscore1) + (0.44 * MEDF_DIF_iw$PFscore2)
  
  # separated MEDF, same
  dat <-MEDF_SAME_iw[,c(grep("_CORR",colnames(MEDF_SAME_iw)),grep("_TTR",colnames(MEDF_SAME_iw)))]
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
  MEDF_SAME_iw <- data.frame(MEDF_SAME_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(MEDF_SAME_iw)[(ncol(MEDF_SAME_iw)-2):ncol(MEDF_SAME_iw)] <- c("PFscore1","PFscore2","SMVE")
  # both 0% and 100% scores need adjusting
  
  MEDF_SAME_iw$PFscore1 <- ifelse(MEDF_SAME_iw$MEDF36_A.MEDF36A_SAME_CR == 0, min(MEDF_SAME_iw$PFscore1,na.rm = T),
                                  ifelse(MEDF_SAME_iw$MEDF36_A.MEDF36A_SAME_CR == 4, max(MEDF_SAME_iw$PFscore1,na.rm = T),MEDF_SAME_iw$PFscore1))
  MEDF_SAME_iw$PFscore2 <- ifelse(MEDF_SAME_iw$MEDF36_A.MEDF36A_SAME_CR == 0, min(MEDF_SAME_iw$PFscore2,na.rm = T),
                                  ifelse(MEDF_SAME_iw$MEDF36_A.MEDF36A_SAME_CR == 4, max(MEDF_SAME_iw$PFscore2,na.rm = T),MEDF_SAME_iw$PFscore2))
  MEDF_SAME_iw$SMVE <- (0.34 * MEDF_SAME_iw$outlier_score_2cut) + (0 * MEDF_SAME_iw$acc3e) + (0.22 * MEDF_SAME_iw$PFscore1) + (0.44 * MEDF_SAME_iw$PFscore2)
  
  
  # PMAT
  dat <-PMAT_iw[,c(grep("_CORR",colnames(PMAT_iw)),grep("_TTR",colnames(PMAT_iw)))]
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
  PMAT_iw <- data.frame(PMAT_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(PMAT_iw)[(ncol(PMAT_iw)-2):ncol(PMAT_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  PMAT_iw$PFscore1 <- ifelse(PMAT_iw$PMAT24_A.PMAT24_A_PC == 100, max(PMAT_iw$PFscore1,na.rm = T),PMAT_iw$PFscore1)
  PMAT_iw$PFscore2 <- ifelse(PMAT_iw$PMAT24_A.PMAT24_A_PC == 100, max(PMAT_iw$PFscore2,na.rm = T),PMAT_iw$PFscore2)
  PMAT_iw$SMVE <- (0.34 * PMAT_iw$outlier_score_2cut) + (0 * PMAT_iw$acc3e) + (0.22 * PMAT_iw$PFscore1) + (0.44 * PMAT_iw$PFscore2)
  
  
  # PLOT
  dat <-PLOT_iw[,c(grep("_CORR",colnames(PLOT_iw)),grep("_RT",colnames(PLOT_iw)))]
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
  PLOT_iw <- data.frame(PLOT_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(PLOT_iw)[(ncol(PLOT_iw)-2):ncol(PLOT_iw)] <- c("PFscore1","PFscore2","SMVE")
  # 0% and 100% scores need adjusting
  
  PLOT_iw$PFscore1 <- ifelse(PLOT_iw$VSPLOT15.VSPLOT15_PC == 100, max(PLOT_iw$PFscore1,na.rm = T),
                             ifelse(PLOT_iw$VSPLOT15.VSPLOT15_PC == 0, min(PLOT_iw$PFscore1,na.rm = T),PLOT_iw$PFscore1))
  PLOT_iw$PFscore2 <- ifelse(PLOT_iw$VSPLOT15.VSPLOT15_PC == 100, max(PLOT_iw$PFscore2,na.rm = T),
                             ifelse(PLOT_iw$VSPLOT15.VSPLOT15_PC == 0, min(PLOT_iw$PFscore2,na.rm = T),PLOT_iw$PFscore2))
  PLOT_iw$SMVE <- (0.34 * PLOT_iw$outlier_score_2cut) + (0 * PLOT_iw$acc3e) + (0.22 * PLOT_iw$PFscore1) + (0.44 * PLOT_iw$PFscore2)
  
  
  # PRA - no RT data, no SMVE
  
  
  # PVRT
  dat <-PVRT_iw[,c(grep("_CORR",colnames(PVRT_iw)),grep("_TTR",colnames(PVRT_iw)))]
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
  PVRT_iw <- data.frame(PVRT_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(PVRT_iw)[(ncol(PVRT_iw)-2):ncol(PVRT_iw)] <- c("PFscore1","PFscore2","SMVE")
  # 0% and 100% scores need adjusting
  
  PVRT_iw$PFscore1 <- ifelse(PVRT_iw$SPVRT_A.SPVRTA_PC == 100, max(PVRT_iw$PFscore1,na.rm = T),
                             ifelse(PVRT_iw$SPVRT_A.SPVRTA_PC == 0, min(PVRT_iw$PFscore1,na.rm = T),PVRT_iw$PFscore1))
  PVRT_iw$PFscore2 <- ifelse(PVRT_iw$SPVRT_A.SPVRTA_PC == 100, max(PVRT_iw$PFscore2,na.rm = T),
                             ifelse(PVRT_iw$SPVRT_A.SPVRTA_PC == 0, min(PVRT_iw$PFscore2,na.rm = T),PVRT_iw$PFscore2))
  PVRT_iw$SMVE <- (0.34 * PVRT_iw$outlier_score_2cut) + (0 * PVRT_iw$acc3e) + (0.22 * PVRT_iw$PFscore1) + (0.44 * PVRT_iw$PFscore2)
  
  
  # RDISC 
  dat <- RDISC_iw
  dat <-dat[,c(grep(".q_",colnames(dat)),grep(".trr_",colnames(dat)))]
  dat[,1:41] <- dat[,1:41]-1
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
  RDISC_iw <- data.frame(RDISC_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(RDISC_iw)[(ncol(RDISC_iw)-2):ncol(RDISC_iw)] <- c("PFscore1","PFscore2","SMVE")
  # 0% and 100% scores need adjusting
  
  RDISC_iw$PFscore1 <- ifelse(RDISC_iw$rdisc_sum == 41, max(RDISC_iw$PFscore1,na.rm = T),
                              ifelse(RDISC_iw$rdisc_sum == 0, min(RDISC_iw$PFscore1,na.rm = T),RDISC_iw$PFscore1))
  RDISC_iw$PFscore2 <- ifelse(RDISC_iw$rdisc_sum == 41, max(RDISC_iw$PFscore2,na.rm = T),
                              ifelse(RDISC_iw$rdisc_sum == 0, min(RDISC_iw$PFscore2,na.rm = T),RDISC_iw$PFscore2))
  RDISC_iw$SMVE <- (0.34 * RDISC_iw$outlier_score_2cut) + (0 * RDISC_iw$acc3e) + (0.22 * RDISC_iw$PFscore1) + (0.44 * RDISC_iw$PFscore2)
  # still 1 row with no SMVE because of missing itemwise
  
  
  # VOLT
  dat <-VOLT_iw[,c(grep("_CORR",colnames(VOLT_iw)),grep("_TTR",colnames(VOLT_iw)))]
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
  VOLT_iw <- data.frame(VOLT_iw,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(VOLT_iw)[(ncol(VOLT_iw)-2):ncol(VOLT_iw)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  VOLT_iw$PFscore1 <- ifelse(VOLT_iw$SVOLT_A.SVOLT_CR == 20, max(VOLT_iw$PFscore1,na.rm = T),VOLT_iw$PFscore1)
  VOLT_iw$PFscore2 <- ifelse(VOLT_iw$SVOLT_A.SVOLT_CR == 20, max(VOLT_iw$PFscore2,na.rm = T),VOLT_iw$PFscore2)
  VOLT_iw$SMVE <- (0.42 * VOLT_iw$outlier_score_2cut) + (0.02 * VOLT_iw$acc3e) + (0.05 * VOLT_iw$PFscore1) + (0.50 * VOLT_iw$PFscore2)
  
  
  # VOLT targets
  dat <-VOLT_targets[,c(grep("_CORR",colnames(VOLT_targets)),grep("_TTR",colnames(VOLT_targets)))]
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
  VOLT_targets <- data.frame(VOLT_targets,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(VOLT_targets)[(ncol(VOLT_targets)-2):ncol(VOLT_targets)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  VOLT_targets$PFscore1 <- ifelse(VOLT_targets$SVOLT_A.SVOLT_TP == 10, max(VOLT_targets$PFscore1,na.rm = T),VOLT_targets$PFscore1)
  VOLT_targets$PFscore2 <- ifelse(VOLT_targets$SVOLT_A.SVOLT_TP == 10, max(VOLT_targets$PFscore2,na.rm = T),VOLT_targets$PFscore2)
  VOLT_targets$SMVE <- (0.42 * VOLT_targets$outlier_score_2cut) + (0.02 * VOLT_targets$acc3e) + (0.05 * VOLT_targets$PFscore1) + (0.50 * VOLT_targets$PFscore2)
  
  
  # VOLT foils
  dat <-VOLT_foils[,c(grep("_CORR",colnames(VOLT_foils)),grep("_TTR",colnames(VOLT_foils)))]
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
  VOLT_foils <- data.frame(VOLT_foils,outlier_score_2cut,acc3e,pfit1,pfit2,sc)
  names(VOLT_foils)[(ncol(VOLT_foils)-2):ncol(VOLT_foils)] <- c("PFscore1","PFscore2","SMVE")
  # no 0%, but 100% scores need adjusting
  
  VOLT_foils$PFscore1 <- ifelse(VOLT_foils$SVOLT_A.SVOLT_TN == 10, max(VOLT_foils$PFscore1,na.rm = T),VOLT_foils$PFscore1)
  VOLT_foils$PFscore2 <- ifelse(VOLT_foils$SVOLT_A.SVOLT_TN == 10, max(VOLT_foils$PFscore2,na.rm = T),VOLT_foils$PFscore2)
  VOLT_foils$SMVE <- (0.42 * VOLT_foils$outlier_score_2cut) + (0.02 * VOLT_foils$acc3e) + (0.05 * VOLT_foils$PFscore1) + (0.50 * VOLT_foils$PFscore2)
}



# * * plot out distribution of SMVE subscores as well as score ----

{
  # variable to keep track of whose data not to use
  no_good <- data.frame(BBLID = demo_new_iw$bblid,ADT = NA,AIM = NA,CPF = NA,
                        CPT = NA,CPW = NA,DDISC = NA,DIGSYM = NA,EDISC = NA,
                        ER40 = NA,GNG = NA,MEDF = NA,PLOT = NA,PMAT = NA,
                        PRA = NA,PVRT = NA,RDISC = NA,VOLT = NA)
  
  # ADT distribution of performance validity
  dat <- ADT_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 225 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) + 
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form ADT", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") + 
      # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
      # scale_y_continuous(breaks = seq(0,175000,25000)) +
      coord_flip() 
    
    # pdf("data/outputs/ADT_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum <- data.frame(cutoff = qu,n = n_left)
  no_good[,2] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # AIM distribution of performance
  dat <- AIM_iw
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = AIM.AIMTOT)) + 
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
      theme_minimal() + labs(title = "AIM.TOT Distribution for full-form AIM",
                             x = "", y = "TOTAL CORRECT") + 
      coord_flip() 
    
    # pdf("data/outputs/AIM_TOT_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
    
    
    my_plot1 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = AIM.AIMTOTRT)) + 
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
      theme_minimal() + labs(title = "AIM.TOTRT Distribution for full-form AIM",
                             x = "", y = "TOTRT") + 
      coord_flip() 
    
    # pdf("data/outputs/AIM_TOTRT_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot1
    # dev.off()
  }
  
  # get rid of bottom outlier, AIMTOT < 40
  bad_aim <- dat %>% filter(AIM.AIMTOT < 40) %>% pull(bblid)
  no_good[,3] <- ifelse(no_good$BBLID %in% bad_aim,1,0)
  
  
  # CPF distribution of performance validity
  dat <- CPF_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 185 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) + 
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form CPF", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") + 
      coord_flip() 
    
    # pdf("data/outputs/CPF_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[3,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,4] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # CPT distribution of performance
  dat <- CPT_iw
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SPCPTNL.SCPL_TP)) + 
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
      theme_minimal() + labs(title = "SPCPTNL.SCPL_TP Distribution for full-form CPT",
                             x = "", y = "TOTAL CORRECT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/CPT_TP_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
    
    
    my_plot1 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SPCPTNL.SCPL_TPRT)) + 
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
      theme_minimal() + labs(title = "SPCPTNL.SCPL_TPRT Distribution for full-form CPT",
                             x = "", y = "TOTRT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/CPT_TPRT_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot1
    # dev.off()
  }
  
  # for LRNR, need to reorder everything by trial first
  CPT_new_order <- data.frame(bblid = dat$bblid, dat$SPCPTNL.SCPL_TP,dat$SPCPTNL.SCPL_TPRT,dat[,grepl("TRIAL",colnames(dat))|grepl("RESP",colnames(dat))])
  
  for (i in 1:nrow(CPT_new_order)) {
    row <- CPT_new_order[i,4:363]
    
    row_t <- row %>% dplyr::select(matches("TRIAL"))
    row_r <- row %>% dplyr::select(matches("RESP"))
    
    CPT_new_order[i,grepl("RESP",colnames(CPT_new_order))] <- row_r[,order(row_t)]
  }
  
  CPT_iw$SPCPTNL.SCPT_LRNR <- count_LRNR(CPT_new_order)
  dat <- CPT_iw
  
  # one extreme outlier, code below removes that (LRNR = 25)
  {
    my_plot2 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SPCPTNL.SCPT_LRNR)) + 
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
      theme_minimal() + labs(title = "SPCPTNL.SCPT_LRNR Distribution for full-form CPT",
                             x = "", y = "LRNR") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/CPT_LRNR_dist_noext_220720.pdf",height = 7,width = 10)
    # my_plot2
    # dev.off()
  }
  
  bad_cpt <- dat %>% filter(SPCPTNL.SCPT_LRNR == 25) %>% pull(bblid)
  no_good[,5] <- ifelse(no_good$BBLID %in% bad_cpt,1,0)
  
  
  # CPW distribution of performance validity
  dat <- CPW_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 185 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) + 
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form CPW", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") + 
      coord_flip() 
    
    # pdf("data/outputs/CPW_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[5,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,6] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # DDISC distribution of performance validity
  dat <- DDISC_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 230 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions.siteid, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form DDISC", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/DDISC_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[6,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,7] <- ifelse(no_good$BBLID %in% data_left$BBLID,0,1)
  
  
  # DIGSYM distribution of performance 
  dat <- DIGSYM_iw
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = DIGSYM.DSCOR)) + 
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
      theme_minimal() + labs(title = "DIGSYM.DSCOR Distribution for full-form DIGSYM",
                             x = "", y = "TOTAL CORRECT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/DS_COR_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
    
    
    my_plot1 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = DIGSYM.DSCORRT)) + 
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
      theme_minimal() + labs(title = "DIGSYM.DSCORRT Distribution for full-form DIGSYM",
                             x = "", y = "COR RT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/DS_CORRT_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot1
    # dev.off()
    
    my_plot2 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = DIGSYM.DS_TP)) + 
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
      theme_minimal() + labs(title = "DIGSYM.DS_TP Distribution for full-form DIGSYM",
                             x = "", y = "TOTAL CORRECT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/DS_TP_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot2
    # dev.off()
    
    
    my_plot3 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = DIGSYM.DS_TPRT)) + 
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
      theme_minimal() + labs(title = "DIGSYM.DS_TPRT Distribution for full-form DIGSYM",
                             x = "", y = "TPRT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/DS_TPRT_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot3
    # dev.off()
    
    my_plot4 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = DIGSYM.DSMEMCR)) + 
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
      theme_minimal() + labs(title = "DIGSYM.DSMEMCR Distribution for full-form DIGSYM",
                             x = "", y = "TOTAL CORRECT (memory)") + 
      scale_y_continuous(breaks=seq(0, 10, 1)) +
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/DS_MEMCR_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot4
    # dev.off()
    
    
    my_plot5 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = DIGSYM.DSMCRRT)) + 
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
      theme_minimal() + labs(title = "DIGSYM.DSMCRRT Distribution for full-form DIGSYM",
                             x = "", y = "COR RT (memory)") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/DS_MCRRT_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot5
    # dev.off()
  }
  
  bad_ds <- dat %>% filter(DIGSYM.DSCOR > 90) %>% pull(bblid)
  no_good[,8] <- ifelse(no_good$BBLID %in% bad_ds,1,0)
  
  
  # EDISC distribution of performance validity
  dat <- EDISC_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 229 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions.siteid, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form EDISC", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/EDISC_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[8,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,9] <- ifelse(no_good$BBLID %in% data_left$BBLID,0,1)
  
  
  # ER40 distribution of performance validity
  dat <- ER40_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 225 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form ER40", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/ER40_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[9,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,10] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # GNG distribution of performance 
  dat <- GNG_iw
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = GNG150.GNG150_CR)) + 
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
      theme_minimal() + labs(title = "GNG150.GNG150_CR Distribution for full-form GNG",
                             x = "", y = "TOTAL CORRECT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/GNG_CR_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
    
    
    my_plot1 <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = GNG150.GNG150_RTCR)) + 
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
      theme_minimal() + labs(title = "GNG150.GNG150_RTCR Distribution for full-form DIGSYM",
                             x = "", y = "COR RT") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/GNG_RTCR_dist_6Jul22.pdf",height = 7,width = 10)
    # my_plot1
    # dev.off()
  }
  
  # looking at LRNR
  GNG_iw$GNG150.GNG150_LRNR <- count_LRNR(dat)
  dat <- GNG_iw
  
  # one extreme outlier, code below removes that (LRNR = 150)
  {
    my_plot2 <- ggplot(dat %>% filter(GNG150.GNG150_LRNR<150), aes(x = test_sessions_v.battery_complete, y = GNG150.GNG150_LRNR)) + 
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
      theme_minimal() + labs(title = "GNG150.GNG150_LRNR Distribution for full-form GNG",
                             x = "", y = "LRNR") + 
      coord_flip() 
    
    # pdf("data/outputs/rapid_dist/GNG_LRNR_dist_noext_8Jul22.pdf",height = 7,width = 10)
    # my_plot2
    # dev.off()
  }
  
  bad_gng <- dat %>% filter(GNG150.GNG150_LRNR == 150) %>% pull(bblid)
  no_good[,11] <- ifelse(no_good$BBLID %in% bad_gng,1,0)
  
  
  # MEDF distribution of performance validity
  dat <- MEDF_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 225 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form MEDF", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/MEDF_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[11,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,12] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # PLOT distribution of performance validity
  dat <- PLOT_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 225 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form PLOT", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/PLOT_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[12,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,13] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # PMAT distribution of performance validity
  dat <- PMAT_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 225 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form PMAT", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/PMAT_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[13,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,14] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # PRA distribution of performance validity
  # no clear QC method yet
  
  
  # PVRT distribution of performance validity
  dat <- PVRT_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 225 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form PVRT", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/PVRT_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[15,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,16] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
  
  
  # RDISC distribution of performance validity
  dat <- RDISC_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 265 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions.siteid, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form RDISC", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/RDISC_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[16,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,17] <- ifelse(no_good$BBLID %in% data_left$BBLID,0,1)
  
  
  # VOLT distribution of performance validity
  dat <- VOLT_iw
  n_tot <- dim(dat)[1]
  qu <- quantile(dat$SMVE,0.05,na.rm=TRUE)
  
  data_left <- dat[which(dat$SMVE > qu),]  # 254 left
  n_left <- dim(data_left)[1]
  
  {
    my_plot <- ggplot(dat, aes(x = test_sessions_v.battery_complete, y = SMVE)) +
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
      geom_hline(yintercept = qu) +
      coord_cartesian(xlim = c(1.2, NA)) +
      theme_minimal() + labs(title = "SMVE Distribution for full-form VOLT", caption = paste("n =",dim(dat)[1], "and", n_left, "rows left after removing bottom 5%"),
                             x = "", y = "SMVE") +
      coord_flip()
    
    # pdf("data/outputs/VOLT_SMVE_dist_26Jun22.pdf",height = 7,width = 10)
    # my_plot
    # dev.off()
  }
  
  smve_sum[17,] <- data.frame(cutoff = qu,n = n_left)
  no_good[,18] <- ifelse(no_good$BBLID %in% data_left$bblid,0,1)
}


rownames(smve_sum) <- c("ADT","AIM","CPF","CPT","CPW","DDISC","DIGSYM","EDISC","ER40",
                        "GNG","MEDF","PLOT","PMAT","PRA","PVRT","RDISC","SVOLT")
smve_sum[,1] <- sapply(smve_sum[,1],round,3)

# smve_sum %>%
#   kbl(caption = "SMVE Cutoff value + Rows left per Test", align = rep("c", 8),
#       col.names = c("Cutoff", "Rows Left")) %>%
#   kable_classic(full_width = F, html_font = "Cambria") %>%
#   save_kable(file = "data/outputs/SMVE_table_221005.pdf", self_contained = T)


# sum of all test SMVE_drops
no_good$sum <- rowSums(no_good[,2:18],na.rm = T) # max: 6 tests flagged, 10/3/2022

# NA for rows with SMVE below 5%
{
  lower_names <- tolower(names(no_good[,2:18]))
  
  xx <- all_cnb
  x <- left_join(xx,no_good[,1:18],by=c("bblid"="BBLID"))
  # 1:18 of x are demos, 19:91 are full CNB, 92:[140] are CAT CNB, 141:157 are no_good
  
  for (i in 1:length(lower_names)) {
    tname <- lower_names[i]
    for (j in 1:nrow(x)) {
      if (!is.na(x[j,i+140]) & x[j,i+140] == 1) {   # need to update this number every time more columns are added to all_cnb, use bracketed number from above comment
        x[j,grepl(tname,colnames(x))] <- NA
      }
    }
  }
}


# * Multivariate Outlier Removal ----


# ** Extreme Outliers (seen from scatters) ----
# basically checking the people that would have been flagged in multivariate outlier removal (removing anyone whose difference between Full/CAT scores are > 3 SD)

# outlier_pra <- x_all %>% filter(pra_cr_Oreg < -4) %>% dplyr::select(bblid:proto_4,pra_cr,pra.1.00.d.cat_default,pra_cr_Oreg,pra.1.00.d.cat_default_Oreg)
# i checked the BBLIDs (23230,105176,122277) that have 0 for pra_cr and they all have C1 on their CNB notes, ignoring these peoples' PRA data for now
# x[which(x$pra_cr == 0),grepl("^pra",colnames(x))] <- NA

# also checking the PRA outlier that has a near 0 on Full but very low score on CAT 
# outlier_pra2 <- x_all %>% filter(pra_cr_Oreg > -0.15, pra.1.00.d.cat_default_Oreg < -3.5) %>% dplyr::select(bblid:proto_4,pra_cr,pra.1.00.d.cat_default,pra_cr_Oreg,pra.1.00.d.cat_default_Oreg)
# BBLID 95116: C1 on Full and CAT CNB
# x[which(x$bblid == 95116),grepl("^pra",colnames(x))] <- NA

# outlier_edisc <- x_all %>% filter(edisc_sum_Oreg > 1.85,edisc.1.00.cat_default_Oreg < -1) %>% dplyr::select(bblid:proto_4,edisc_sum,edisc.1.00.cat_default,edisc_sum_Oreg,edisc.1.00.cat_default_Oreg)
# BBLID (22016,23463) both have C1 on Full and CAT CNB notes, therefore ignoring for now

# outlier_gng <- x_all %>% filter(gng_cr_Oreg < -3.9,GNG60.GNG60_CR_Oreg > 0) %>% dplyr::select(bblid:proto_4,gng_cr,GNG60.GNG60_CR,gng_cr_Oreg,GNG60.GNG60_CR_Oreg)
# BBLID 91335: glitches in the Full CNB, C3 good enough to remove? this will be captured by manual validation codes

# outlier_plot <- x_all %>% filter(plot_pc_Oreg < -1.9,plot.1.00.cat_default_Oreg > 1) %>% dplyr::select(bblid:proto_4,plot_pc,plot.1.00.cat_default,plot_pc_Oreg,plot.1.00.cat_default_Oreg)
# BBLID 112598: C1 for both Full and CAT CNB



# Select relevant columns for scatters ----


temp <- x$gng_cr
temp[temp<100] <- NA # 2 turned into NA with this, 10/3/22
x$gng_cr <- temp

temp <- x$GNG60.GNG60_CR
temp[temp<50] <- NA # 2 turned into NA with this, 10/3/22
x$GNG60.GNG60_CR <- temp

cpt_acc <- x$cpt_ptp - x$cpt_pfp

# looking at distributions of cpt_acc
{
  temp_cpt_dat <- data.frame(cpt_acc,temp = 1)
  my_plot <- ggplot(temp_cpt_dat, aes(x = temp, y = cpt_acc)) + 
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
    theme_minimal() + labs(title = "cpt_acc Distribution for full-form CPT",
                           x = "", y = "CPT Accuracy") + 
    coord_flip() 
  
  # pdf("data/outputs/rapid_dist/CPT_acc_dist_7Jul22.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
}

acc <- x %>% dplyr::select(matches("_cr$|tot$|cor$|pc$|memcr$|sum$",ignore.case = F)) %>% cbind(cpt_acc) # all acc scores from full CNB

er40_cat <- ((x$er40.1.00.cat_neutral + x$er40.1.00.cat_emotive + x$er40.1.00.cat_emotive + x$er40.1.00.cat_emotive)/4)*sqrt(4)
medf_cat <- ((x$medf.1.00.cat_same + x$medf.1.00.cat_different + x$medf.1.00.cat_different + x$medf.1.00.cat_different)/4)*sqrt(4)
adt_cat <- ((x$adt.1.00.cat_same + x$adt.1.00.cat_different + x$adt.1.00.cat_different + x$adt.1.00.cat_different)/4)*sqrt(4)

# new_er40_cat <- ((x$er40.1.00.cat_neutral)*(2/5) + (x$er40.1.00.cat_emotive)*(3/5))*sqrt(5)
# new_medf_cat <- ((x$medf.1.00.cat_same)*(5/17) + (x$medf.1.00.cat_different)*(12/17))*sqrt(17)
# new_adt_cat <- ((x$adt.1.00.cat_same)*(5/17) + (x$adt.1.00.cat_different)*(12/17))*sqrt(17)

cpf_cat <- ((x$cpf.1.00.v1.cat_target + x$cpf.1.00.v1.cat_foil)/2)*sqrt(2)
cpw_cat <- ((x$cpw.1.00.v1.cat_target + x$cpw.1.00.v1.cat_foil)/2)*sqrt(2)
volt_cat <- ((x$volt.1.00.v1.cat_targets + x$volt.1.00.v1.cat_foils)/2)*sqrt(2)
cpt_cat <- (x$CPT108.CATCPTT_TP/36) - (x$CPT108.CATCPTT_FP/72)

# using the "real"/fixed cpf and er40 (still n=62 only)
er40_2_cat <- ((x$er40.2.00.cat_neutral + x$er40.2.00.cat_emotive + x$er40.2.00.cat_emotive + x$er40.2.00.cat_emotive)/4)*sqrt(4)
cpf_2_cat <- ifelse(!is.na(x$cpf2.1.00.v2.cat_target),
                  ((x$cpf2.1.00.v2.cat_target + x$cpf2.1.00.v2.cat_foil)/2)*sqrt(2),
                  ((x$cpf2.1.00.v1.cat_target + x$cpf2.1.00.v1.cat_foil)/2)*sqrt(2))

cpf_2_t <- ifelse(!is.na(x$cpf2.1.00.v2.cat_target),x$cpf2.1.00.v2.cat_target,x$cpf2.1.00.v1.cat_target)
cpf_2_f <- ifelse(!is.na(x$cpf2.1.00.v2.cat_foil),x$cpf2.1.00.v2.cat_foil,x$cpf2.1.00.v1.cat_foil)

# cpf2.1.00.v1.cat_foil is the cat_v2 (full fixed link data, )

cat_acc <- x %>% dplyr::select(matches("default|TOT$|_CR$|DSCOR$|MEMCR",ignore.case = F)) # all CAT CNB task scores that aren't split


x99 <- data.frame(x %>% dplyr::select(bblid,age_enroll,study_group,sex,proto_3,proto_4),acc,
                  cat_acc,er40_cat,medf_cat,adt_cat,cpf_cat,cpw_cat,volt_cat,er40_2_cat,cpf_2_cat,cpt_cat, 
                  x %>% dplyr::select(er40.1.00.cat_emotive,er40.1.00.cat_neutral,medf.1.00.cat_same,medf.1.00.cat_different,
                                      er40.2.00.cat_emotive,er40.2.00.cat_neutral,cpw.1.00.v1.cat_target,cpw.1.00.v1.cat_foil,
                                      volt.1.00.v1.cat_targets,volt.1.00.v1.cat_foils),
                  cpf_2_t,cpf_2_f)

# temporarily (no QC-ing whatsoever of CAT vs full data to match un-QC-ed extralong full vs full)
x_xl <- x99

# split tasks: need to add separated memory tasks 
{
  temp_er40_emo <- ER40_EMO_iw %>% dplyr::select(matches("^bblid$|er40_emo.ER40_D.ER40D_EMO|er40_resp.ER40_D.ER40D_EMORTCR|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|ER40")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_er40_emo$SMVE,0.05,na.rm=TRUE)
  temp_er40_emo$er40_emo_splitSMVE <- ifelse(temp_er40_emo$SMVE > qu,temp_er40_emo$er40_emo.ER40_D.ER40D_EMO,NA)
  temp_er40_emo$er40_emo_SMVE <- ifelse(temp_er40_emo$ER40 != 1,temp_er40_emo$er40_emo.ER40_D.ER40D_EMO,NA)
  
  temp_er40_neu <- ER40_NEU_iw %>% dplyr::select(matches("^bblid$|ER40_D.ER40D_NOE|ER40_D.ER40D_NOERT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|ER40")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_er40_neu$SMVE,0.05,na.rm=TRUE)
  temp_er40_neu$er40_noe_splitSMVE <- ifelse(temp_er40_neu$SMVE > qu,temp_er40_neu$ER40_D.ER40D_NOE,NA)
  temp_er40_neu$er40_noe_SMVE <- ifelse(temp_er40_neu$ER40 != 1,temp_er40_neu$ER40_D.ER40D_NOE,NA)
  
  temp_medf_dif <- MEDF_DIF_iw %>% dplyr::select(matches("^bblid$|medf_corr.MEDF36_A.MEDF36A_DIF|medf_resp.MEDF36_A.MEDF36A_DIFRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|MEDF")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_medf_dif$SMVE,0.05,na.rm=TRUE)
  temp_medf_dif$medf_dif_splitSMVE <- ifelse(temp_medf_dif$SMVE > qu,temp_medf_dif$medf_corr.MEDF36_A.MEDF36A_DIF,NA)
  temp_medf_dif$medf_dif_SMVE <- ifelse(temp_medf_dif$MEDF != 1,temp_medf_dif$medf_corr.MEDF36_A.MEDF36A_DIF,NA)
  
  temp_medf_same <- MEDF_SAME_iw %>% dplyr::select(matches("^bblid$|MEDF36_A.MEDF36A_SAME_CR|MEDF36_A.MEDF36A_SAME_RTCR|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|MEDF")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_medf_same$SMVE,0.05,na.rm=TRUE)
  temp_medf_same$medf_same_splitSMVE <- ifelse(temp_medf_same$SMVE > qu,temp_medf_same$MEDF36_A.MEDF36A_SAME_CR,NA)
  temp_medf_same$medf_same_SMVE <- ifelse(temp_medf_same$MEDF != 1,temp_medf_same$MEDF36_A.MEDF36A_SAME_CR,NA)
  
  split_er40_medf <- left_join(all_cnb %>% dplyr::select("bblid"),temp_er40_emo %>% dplyr::select(matches("bblid|er40_emo_splitSMVE|er40_emo_SMVE")),by="bblid")
  split_er40_medf <- left_join(split_er40_medf,temp_er40_neu %>% dplyr::select(matches("bblid|er40_noe_splitSMVE|er40_noe_SMVE")),by="bblid")
  split_er40_medf <- left_join(split_er40_medf,temp_medf_dif %>% dplyr::select(matches("bblid|medf_dif_splitSMVE|medf_dif_SMVE")),by="bblid")
  split_er40_medf <- left_join(split_er40_medf,temp_medf_same %>% dplyr::select(matches("bblid|medf_same_splitSMVE|medf_same_SMVE")),by="bblid")
  
  x99_split <- left_join(x99,split_er40_medf,by="bblid")
  
  
  temp_cpf_t <- CPF_targets %>% dplyr::select(matches("^bblid$|CPF_B.CPF_TP|CPF_B.CPF_TPRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|CPF")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_cpf_t$SMVE,0.05,na.rm=TRUE)
  temp_cpf_t$cpf_t_splitSMVE <- ifelse(temp_cpf_t$SMVE > qu,temp_cpf_t$CPF_B.CPF_TP,NA)
  temp_cpf_t$cpf_t_SMVE <- ifelse(temp_cpf_t$CPF != 1,temp_cpf_t$CPF_B.CPF_TP,NA)
  
  temp_cpf_f <- CPF_foils %>% dplyr::select(matches("^bblid$|CPF_B.CPF_TN|CPF_B.CPF_TNRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|CPF")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_cpf_f$SMVE,0.05,na.rm=TRUE)
  temp_cpf_f$cpf_f_splitSMVE <- ifelse(temp_cpf_f$SMVE > qu,temp_cpf_f$CPF_B.CPF_TN,NA)
  temp_cpf_f$cpf_f_SMVE <- ifelse(temp_cpf_f$CPF != 1,temp_cpf_f$CPF_B.CPF_TN,NA)
  
  temp_cpw_t <- CPW_targets %>% dplyr::select(matches("^bblid$|CPW_A.CPW_TP|CPW_A.CPW_TPRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|CPW")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_cpw_t$SMVE,0.05,na.rm=TRUE)
  temp_cpw_t$cpw_t_splitSMVE <- ifelse(temp_cpw_t$SMVE > qu,temp_cpw_t$CPW_A.CPW_TP,NA)
  temp_cpw_t$cpw_t_SMVE <- ifelse(temp_cpw_t$CPW != 1,temp_cpw_t$CPW_A.CPW_TP,NA)
  
  temp_cpw_f <- CPW_foils %>% dplyr::select(matches("^bblid$|CPW_A.CPW_TN|CPW_A.CPW_TNRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|CPW")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_cpw_f$SMVE,0.05,na.rm=TRUE)
  temp_cpw_f$cpw_f_splitSMVE <- ifelse(temp_cpw_f$SMVE > qu,temp_cpw_f$CPW_A.CPW_TN,NA)
  temp_cpw_f$cpw_f_SMVE <- ifelse(temp_cpw_f$CPW != 1,temp_cpw_f$CPW_A.CPW_TN,NA)
  
  temp_volt_t <- VOLT_targets %>% dplyr::select(matches("^bblid$|SVOLT_A.SVOLT_TP|SVOLT_A.SVOLT_TPRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|VOLT")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_volt_t$SMVE,0.05,na.rm=TRUE)
  temp_volt_t$volt_t_splitSMVE <- ifelse(temp_volt_t$SMVE > qu,temp_volt_t$SVOLT_A.SVOLT_TP,NA)
  temp_volt_t$volt_t_SMVE <- ifelse(temp_volt_t$VOLT != 1,temp_volt_t$SVOLT_A.SVOLT_TP,NA)
  
  temp_volt_f <- VOLT_foils %>% dplyr::select(matches("^bblid$|SVOLT_A.SVOLT_TN|SVOLT_A.SVOLT_TNRT|SMVE")) %>% 
    filter(bblid %in% all_cnb$bblid) %>% left_join(no_good %>% dplyr::select(matches("BBLID|VOLT")),by=c("bblid" = "BBLID"))
  qu <- quantile(temp_volt_f$SMVE,0.05,na.rm=TRUE)
  temp_volt_f$volt_f_splitSMVE <- ifelse(temp_volt_f$SMVE > qu,temp_volt_f$SVOLT_A.SVOLT_TN,NA)
  temp_volt_f$volt_f_SMVE <- ifelse(temp_volt_f$VOLT != 1,temp_volt_f$SVOLT_A.SVOLT_TN,NA)
  
  split_mem <- left_join(all_cnb %>% dplyr::select("bblid"),temp_cpf_t %>% dplyr::select(matches("bblid|cpf_t_splitSMVE|cpf_t_SMVE")),by="bblid")
  split_mem <- left_join(split_mem,temp_cpf_f %>% dplyr::select(matches("bblid|cpf_f_splitSMVE|cpf_f_SMVE")),by="bblid")
  split_mem <- left_join(split_mem,temp_cpw_t %>% dplyr::select(matches("bblid|cpw_t_splitSMVE|cpw_t_SMVE")),by="bblid")
  split_mem <- left_join(split_mem,temp_cpw_f %>% dplyr::select(matches("bblid|cpw_f_splitSMVE|cpw_f_SMVE")),by="bblid")
  split_mem <- left_join(split_mem,temp_volt_t %>% dplyr::select(matches("bblid|volt_t_splitSMVE|volt_t_SMVE")),by="bblid")
  split_mem <- left_join(split_mem,temp_volt_f %>% dplyr::select(matches("bblid|volt_f_splitSMVE|volt_f_SMVE")),by="bblid")
  
  x99_split <- left_join(x99_split,split_mem,by="bblid")
}

# order regress
{
  sc <- matrix(NA,nrow(x99_split),ncol(x99_split)-6)
  
  for (i in 1:(ncol(x99_split)-6)) {
    mod <- lm(x99_split[,(i+6)]~proto_3,data=x99_split,na.action=na.exclude)
    sc[,i] <- scale(residuals(mod,na.action=na.exclude))
  }
  
  colnames(sc) <- paste0(colnames(x99_split[,7:ncol(x99_split)]),"_Oreg")
}

x_all <- data.frame(x99_split,sc)
x_TD <- data.frame(x99_split,sc) %>% filter(study_group %in% c("Healthy Controls"))
x_PS <- data.frame(x99_split,sc) %>% filter(study_group %in% c("Psychosis"))
x_MD <- data.frame(x99_split,sc) %>% filter(study_group %in% c("Mood-Anx-BP"))

# write.csv(x_all,"CNB-CAT_test-retest_with_order-regressed_221003.csv",na="",row.names=FALSE)


# scatters spread by test
{
  pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_ALL_220715_test.pdf",height=9,width=12)
  # par(mfrow=c(2,2))   # still need to figure out how to have all four scatters in one plot
  pairs.panels(x %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(x_TD %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(x_PS %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
  pairs.panels(x_MD %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
  dev.off()
}

# * all individual tests printed out by condition ---- 

# overall
pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_ALL_bytest_221005.pdf",height=9,width=12)
pairs.panels(x_all %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("aim_tot_Oreg|S_AIM.AIMTOT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
# pairs.panels(x_all %>% dplyr::select(matches("cpf_cr_Oreg|cpf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("cpf_cr_Oreg|cpf_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_all %>% dplyr::select(matches("cpf_t_SMVE_Oreg|cpf_2_t_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_all %>% dplyr::select(matches("cpf_f_SMVE_Oreg|cpf_2_f_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
pairs.panels(x_all %>% dplyr::select(matches("cpt_acc_Oreg|cpt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_all %>% dplyr::select(matches("cpw_cr_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_all %>% dplyr::select(matches("cpw_t_SMVE_Oreg|cpw.1.00.v1.cat_target_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_all %>% dplyr::select(matches("cpw_f_SMVE_Oreg|cpw.1.00.v1.cat_foil_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("ddisc_sum_Oreg|ddisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("dscor_Oreg|S_DIGSYM.DSCOR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_all %>% dplyr::select(matches("dsmemcr_Oreg|S_DIGSYM.DSMEMCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_all %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% filter(bblid %notin% c(22016,23463)) %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_all %>% dplyr::select(matches("er40_cr_Oreg|er40_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("er40_cr_Oreg|er40_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                         # new, corrected er40 (er40v2)
# pairs.panels(x_all %>% dplyr::select(matches("er40_emo_SMVE_Oreg|er40.2.00.cat_emotive_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, emotive items, by general SMVE
# pairs.panels(x_all %>% dplyr::select(matches("er40_noe_SMVE_Oreg|er40.2.00.cat_neutral_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, neutral items, by general SMVE
pairs.panels(x_all %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_all %>% filter(bblid !=91335) %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("medf_pc_Oreg|medf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_all %>% dplyr::select(matches("medf_dif_SMVE_Oreg|medf.1.00.cat_different_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)      # medf split, dif items, by general SMVE
# pairs.panels(x_all %>% dplyr::select(matches("medf_same_SMVE_Oreg|medf.1.00.cat_same_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)          # medf split, same items, by general SMVE
pairs.panels(x_all %>% dplyr::select(matches("plot_pc_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("pmat_pc_Oreg|pmat.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra
pairs.panels(x_all %>% filter(pra_cr!=0, bblid!=95116) %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra with pra_cr==0 removed
pairs.panels(x_all %>% dplyr::select(matches("pvrt_cr_Oreg|pvrt.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("rdisc_sum_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("volt_cr_Oreg|volt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_all %>% dplyr::select(matches("volt_t_SMVE_Oreg|volt.1.00.v1.cat_targets_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_all %>% dplyr::select(matches("volt_f_SMVE_Oreg|volt.1.00.v1.cat_foils_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# TD
pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_TD_bytest_221005.pdf",height=9,width=12)
pairs.panels(x_TD %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("aim_tot_Oreg|S_AIM.AIMTOT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
# pairs.panels(x_TD %>% dplyr::select(matches("cpf_cr_Oreg|cpf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("cpf_cr_Oreg|cpf_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_TD %>% dplyr::select(matches("cpf_t_SMVE_Oreg|cpf_2_t_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_TD %>% dplyr::select(matches("cpf_f_SMVE_Oreg|cpf_2_f_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
pairs.panels(x_TD %>% dplyr::select(matches("cpt_acc_Oreg|cpt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_TD %>% dplyr::select(matches("cpw_cr_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_TD %>% dplyr::select(matches("cpw_t_SMVE_Oreg|cpw.1.00.v1.cat_target_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_TD %>% dplyr::select(matches("cpw_f_SMVE_Oreg|cpw.1.00.v1.cat_foil_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("ddisc_sum_Oreg|ddisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("dscor_Oreg|S_DIGSYM.DSCOR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_TD %>% dplyr::select(matches("dsmemcr_Oreg|S_DIGSYM.DSMEMCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_TD %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% filter(bblid %notin% c(22016,23463)) %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_TD %>% dplyr::select(matches("er40_cr_Oreg|er40_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("er40_cr_Oreg|er40_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                         # new, corrected er40 (er40v2)
# pairs.panels(x_TD %>% dplyr::select(matches("er40_emo_SMVE_Oreg|er40.2.00.cat_emotive_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, emotive items, by general SMVE
# pairs.panels(x_TD %>% dplyr::select(matches("er40_noe_SMVE_Oreg|er40.2.00.cat_neutral_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, neutral items, by general SMVE
pairs.panels(x_TD %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_TD %>% filter(bblid !=91335) %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("medf_pc_Oreg|medf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_TD %>% dplyr::select(matches("medf_dif_SMVE_Oreg|medf.1.00.cat_different_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)      # medf split, dif items, by general SMVE
# pairs.panels(x_TD %>% dplyr::select(matches("medf_same_SMVE_Oreg|medf.1.00.cat_same_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)          # medf split, same items, by general SMVE
pairs.panels(x_TD %>% dplyr::select(matches("plot_pc_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("pmat_pc_Oreg|pmat.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra
pairs.panels(x_TD %>% filter(pra_cr!=0, bblid!=95116) %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra with pra_cr==0 removed
pairs.panels(x_TD %>% dplyr::select(matches("pvrt_cr_Oreg|pvrt.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("rdisc_sum_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("volt_cr_Oreg|volt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_TD %>% dplyr::select(matches("volt_t_SMVE_Oreg|volt.1.00.v1.cat_targets_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_TD %>% dplyr::select(matches("volt_f_SMVE_Oreg|volt.1.00.v1.cat_foils_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# PS
pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_PS_bytest_221005.pdf",height=9,width=12)
pairs.panels(x_PS %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("aim_tot_Oreg|S_AIM.AIMTOT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
# pairs.panels(x_PS %>% dplyr::select(matches("cpf_cr_Oreg|cpf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("cpf_cr_Oreg|cpf_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_PS %>% dplyr::select(matches("cpf_t_SMVE_Oreg|cpf_2_t_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_PS %>% dplyr::select(matches("cpf_f_SMVE_Oreg|cpf_2_f_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
pairs.panels(x_PS %>% dplyr::select(matches("cpt_acc_Oreg|cpt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_PS %>% dplyr::select(matches("cpw_cr_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_PS %>% dplyr::select(matches("cpw_t_SMVE_Oreg|cpw.1.00.v1.cat_target_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_PS %>% dplyr::select(matches("cpw_f_SMVE_Oreg|cpw.1.00.v1.cat_foil_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("ddisc_sum_Oreg|ddisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("dscor_Oreg|S_DIGSYM.DSCOR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_PS %>% dplyr::select(matches("dsmemcr_Oreg|S_DIGSYM.DSMEMCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_PS %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% filter(bblid %notin% c(22016,23463)) %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_PS %>% dplyr::select(matches("er40_cr_Oreg|er40_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("er40_cr_Oreg|er40_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                         # new, corrected er40 (er40v2)
# pairs.panels(x_PS %>% dplyr::select(matches("er40_emo_SMVE_Oreg|er40.2.00.cat_emotive_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, emotive items, by general SMVE
# pairs.panels(x_PS %>% dplyr::select(matches("er40_noe_SMVE_Oreg|er40.2.00.cat_neutral_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, neutral items, by general SMVE
pairs.panels(x_PS %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_PS %>% filter(bblid !=91335) %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("medf_pc_Oreg|medf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_PS %>% dplyr::select(matches("medf_dif_SMVE_Oreg|medf.1.00.cat_different_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)      # medf split, dif items, by general SMVE
# pairs.panels(x_PS %>% dplyr::select(matches("medf_same_SMVE_Oreg|medf.1.00.cat_same_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)          # medf split, same items, by general SMVE
pairs.panels(x_PS %>% dplyr::select(matches("plot_pc_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("pmat_pc_Oreg|pmat.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra
pairs.panels(x_PS %>% filter(pra_cr!=0, bblid!=95116) %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra with pra_cr==0 removed
pairs.panels(x_PS %>% dplyr::select(matches("pvrt_cr_Oreg|pvrt.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("rdisc_sum_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("volt_cr_Oreg|volt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_PS %>% dplyr::select(matches("volt_t_SMVE_Oreg|volt.1.00.v1.cat_targets_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_PS %>% dplyr::select(matches("volt_f_SMVE_Oreg|volt.1.00.v1.cat_foils_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# MD
pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_MD_bytest_221005.pdf",height=9,width=12)
pairs.panels(x_MD %>% dplyr::select(matches("adt_pc_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("aim_tot_Oreg|S_AIM.AIMTOT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
# pairs.panels(x_MD %>% dplyr::select(matches("cpf_cr_Oreg|cpf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("cpf_cr_Oreg|cpf_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_MD %>% dplyr::select(matches("cpf_t_SMVE_Oreg|cpf_2_t_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
# pairs.panels(x_MD %>% dplyr::select(matches("cpf_f_SMVE_Oreg|cpf_2_f_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                           # new, corrected cpf (cpfv2)
pairs.panels(x_MD %>% dplyr::select(matches("cpt_acc_Oreg|cpt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_MD %>% dplyr::select(matches("cpw_cr_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_MD %>% dplyr::select(matches("cpw_t_SMVE_Oreg|cpw.1.00.v1.cat_target_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_MD %>% dplyr::select(matches("cpw_f_SMVE_Oreg|cpw.1.00.v1.cat_foil_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("ddisc_sum_Oreg|ddisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("dscor_Oreg|S_DIGSYM.DSCOR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_MD %>% dplyr::select(matches("dsmemcr_Oreg|S_DIGSYM.DSMEMCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_MD %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% filter(bblid %notin% c(22016,23463)) %>% dplyr::select(matches("edisc_sum_Oreg|edisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_MD %>% dplyr::select(matches("er40_cr_Oreg|er40_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("er40_cr_Oreg|er40_2_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                         # new, corrected er40 (er40v2)
# pairs.panels(x_MD %>% dplyr::select(matches("er40_emo_SMVE_Oreg|er40.2.00.cat_emotive_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, emotive items, by general SMVE
# pairs.panels(x_MD %>% dplyr::select(matches("er40_noe_SMVE_Oreg|er40.2.00.cat_neutral_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)        # er40 split, neutral items, by general SMVE
pairs.panels(x_MD %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_MD %>% filter(bblid !=91335) %>% dplyr::select(matches("gng_cr_Oreg|GNG60.GNG60_CR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("medf_pc_Oreg|medf_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_MD %>% dplyr::select(matches("medf_dif_SMVE_Oreg|medf.1.00.cat_different_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)      # medf split, dif items, by general SMVE
# pairs.panels(x_MD %>% dplyr::select(matches("medf_same_SMVE_Oreg|medf.1.00.cat_same_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)          # medf split, same items, by general SMVE
pairs.panels(x_MD %>% dplyr::select(matches("plot_pc_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("pmat_pc_Oreg|pmat.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra
pairs.panels(x_MD %>% filter(pra_cr!=0, bblid!=95116) %>% dplyr::select(matches("pra_cr_Oreg|pra.1.00.d.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # pra with pra_cr==0 removed
pairs.panels(x_MD %>% dplyr::select(matches("pvrt_cr_Oreg|pvrt.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("rdisc_sum_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("volt_cr_Oreg|volt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_MD %>% dplyr::select(matches("volt_t_SMVE_Oreg|volt.1.00.v1.cat_targets_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_MD %>% dplyr::select(matches("volt_f_SMVE_Oreg|volt.1.00.v1.cat_foils_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()




# for AdaptiveV vs Extralong (XL) comparison
# overall
pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_ALL_bytest_forXL_221005.pdf",height=9,width=12)
pairs.panels(x_xl %>% dplyr::select(matches("adt_pc|adt_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("aim_tot|S_AIM.AIMTOT")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
# pairs.panels(x_xl %>% dplyr::select(matches("cpf_cr|cpf_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("cpf_cr|cpf_2_cat")),lm=TRUE,scale=TRUE,ci=TRUE)  # new, corrected cpf (cpfv2)
pairs.panels(x_xl %>% dplyr::select(matches("cpt_acc|cpt_cat")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_xl %>% dplyr::select(matches("cpw_cr|cpw_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("ddisc_sum|ddisc.1.00.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("dscor|S_DIGSYM.DSCOR")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_xl %>% dplyr::select(matches("dsmemcr|S_DIGSYM.DSMEMCR")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_xl %>% dplyr::select(matches("edisc_sum|edisc.1.00.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(x_xl %>% dplyr::select(matches("er40_cr|er40_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("er40_cr|er40_2_cat")),lm=TRUE,scale=TRUE,ci=TRUE)                         # new, corrected er40 (er40v2)
pairs.panels(x_xl %>% dplyr::select(matches("gng_cr|GNG60.GNG60_CR")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_xl %>% dplyr::select(matches("medf_pc|medf_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("plot_pc|plot.1.00.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("pmat_pc|pmat.1.00.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("pra_cr|pra.1.00.d.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE) # pra
pairs.panels(x_xl %>% dplyr::select(matches("pvrt_cr|pvrt.1.00.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("rdisc_sum|rdisc.1.00.cat_default")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_xl %>% dplyr::select(matches("volt_cr|volt_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# table to accompany above plots
adapt_XL <- data.frame(matrix(NA,ncol = 1,nrow = 19))
rownames(adapt_XL) <- c("ADT","AIM","CPF","CPF v2","CPT","CPW","DDISC","DIGSYM","EDISC","ER40","ER40 v2",
                        "GNG","MEDF","PLOT","PMAT","PRA","PVRT","RDISC","SVOLT")

adapt_XL[1,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("adt_pc")))),sum(!is.na(x_xl %>% dplyr::select(matches("adt_cat")))))
adapt_XL[2,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("aim_tot")))),sum(!is.na(x_xl %>% dplyr::select(matches("S_AIM.AIMTOT")))))
adapt_XL[3,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("cpf_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("cpf_cat")))))
adapt_XL[4,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("cpf_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("cpf_2_cat")))))
adapt_XL[5,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("cpt_acc")))),sum(!is.na(x_xl %>% dplyr::select(matches("cpt_cat")))))
adapt_XL[6,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("cpw_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("cpw_cat")))))
adapt_XL[7,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("ddisc_sum")))),sum(!is.na(x_xl %>% dplyr::select(matches("ddisc.1.00.cat_default")))))
adapt_XL[8,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("dscor")))),sum(!is.na(x_xl %>% dplyr::select(matches("S_DIGSYM.DSCOR")))))
adapt_XL[9,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("edisc_sum")))),sum(!is.na(x_xl %>% dplyr::select(matches("edisc.1.00.cat_default")))))
adapt_XL[10,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("er40_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("er40_cat")))))
adapt_XL[11,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("er40_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("er40_2_cat")))))
adapt_XL[12,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("gng_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("GNG60.GNG60_CR")))))
adapt_XL[13,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("medf_pc")))),sum(!is.na(x_xl %>% dplyr::select(matches("medf_cat")))))
adapt_XL[14,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("plot_pc")))),sum(!is.na(x_xl %>% dplyr::select(matches("plot.1.00.cat_default")))))
adapt_XL[15,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("pmat_pc")))),sum(!is.na(x_xl %>% dplyr::select(matches("pmat.1.00.cat_default")))))
adapt_XL[16,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("pra_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("pra.1.00.d.cat_default")))))
adapt_XL[17,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("pvrt_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("pvrt.1.00.cat_default")))))
adapt_XL[18,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("rdisc_sum")))),sum(!is.na(x_xl %>% dplyr::select(matches("rdisc.1.00.cat_default")))))
adapt_XL[19,1] <- min(sum(!is.na(x_xl %>% dplyr::select(matches("volt_cr")))),sum(!is.na(x_xl %>% dplyr::select(matches("volt_cat")))))

# adapt_XL %>% 
#   kbl(caption = "Number of Rows for each Test", align = rep("c", 8),
#       col.names = "N") %>%
#   kable_classic(full_width = F, html_font = "Cambria") %>%
#   column_spec(1, width = "12em") %>% 
#   save_kable(file = "data/outputs/AdaptiveV_table_220810.pdf", self_contained = T)







# print all tests separated by category of test



pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_ALL_220707.pdf",height=9,width=12)
pairs.panels(x %>% dplyr::select(matches("pvrt_pc_Oreg|pmat_pc_Oreg|plot_pc_Oreg|pvrt.1.00.cat_default_Oreg|pmat.1.00.cat_default_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # pvrt, pmat, plot
pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # er40, medf, adt
# pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|dsmemcr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg|dsmemcr_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, dsmemcr
pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, no dsmemcr for now
# pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # no dscor, gng, aim, cpt for now
pairs.panels(x %>% dplyr::select(matches("ddisc_sum_Oreg|edisc_sum_Oreg|rdisc_sum_Oreg|ddisc.1.00.cat_default_Oreg|edisc.1.00.cat_default_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # ddisc, edisc, rdisc
dev.off()


# healthy controls only
x <- data.frame(x99,sc) %>% filter(study_group == "Healthy Controls")

pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_CONTROLS_220707.pdf",height=9,width=12)
pairs.panels(x %>% dplyr::select(matches("pvrt_pc_Oreg|pmat_pc_Oreg|plot_pc_Oreg|pvrt.1.00.cat_default_Oreg|pmat.1.00.cat_default_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # pvrt, pmat, plot
pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # er40, medf, adt
# pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|dsmemcr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg|dsmemcr_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, dsmemcr
pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, no dsmemcr for now
# pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # no dscor, gng, aim, cpt for now
pairs.panels(x %>% dplyr::select(matches("ddisc_sum_Oreg|edisc_sum_Oreg|rdisc_sum_Oreg|ddisc.1.00.cat_default_Oreg|edisc.1.00.cat_default_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # ddisc, edisc, rdisc
dev.off()



# PS only
x <- data.frame(x99,sc) %>% filter(study_group == "Psychosis")

pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_PS_220707.pdf",height=9,width=12)
pairs.panels(x %>% dplyr::select(matches("pvrt_pc_Oreg|pmat_pc_Oreg|plot_pc_Oreg|pvrt.1.00.cat_default_Oreg|pmat.1.00.cat_default_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # pvrt, pmat, plot
pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # er40, medf, adt
# pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|dsmemcr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg|dsmemcr_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, dsmemcr
pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, no dsmemcr for now
# pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # no dscor, gng, aim, cpt for now
pairs.panels(x %>% dplyr::select(matches("ddisc_sum_Oreg|edisc_sum_Oreg|rdisc_sum_Oreg|ddisc.1.00.cat_default_Oreg|edisc.1.00.cat_default_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # ddisc, edisc, rdisc
dev.off()



# Mood only
x <- data.frame(x99,sc) %>% filter(study_group == "Mood-Anx-BP")

pdf("data/outputs/scatters/CNB-CAT_test-retest_scatter_matrices_MOOD_220707.pdf",height=9,width=12)
pairs.panels(x %>% dplyr::select(matches("pvrt_pc_Oreg|pmat_pc_Oreg|plot_pc_Oreg|pvrt.1.00.cat_default_Oreg|pmat.1.00.cat_default_Oreg|plot.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # pvrt, pmat, plot
pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # er40, medf, adt
# pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|dsmemcr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg|dsmemcr_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, dsmemcr
pairs.panels(x %>% dplyr::select(matches("volt_cr_Oreg|cpf_cr_Oreg|cpw_cr_Oreg|volt_cat_Oreg|cpf_cat_Oreg|cpw_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # volt, cpf, cpw, no dsmemcr for now
# pairs.panels(x %>% dplyr::select(matches("er40_cr_Oreg|medf_pc_Oreg|adt_pc_Oreg|er40_cat_Oreg|medf_cat_Oreg|adt_cat_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # no dscor, gng, aim, cpt for now
pairs.panels(x %>% dplyr::select(matches("ddisc_sum_Oreg|edisc_sum_Oreg|rdisc_sum_Oreg|ddisc.1.00.cat_default_Oreg|edisc.1.00.cat_default_Oreg|rdisc.1.00.cat_default_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)   # ddisc, edisc, rdisc
dev.off()





# ICCs ----

for_icc <- read.csv("data/inputs/crowd_sourced_norms/for_ICCs.csv",na.strings=c(""," ","NA"))


# * scaling full CNB measures ----
# order regress after scaling

# ADT (adt_pc, 36 total) 
temp <- x %>% dplyr::select(adt_pc)
temp$adt_cr <- temp*36*.01
x$adt_cr_scaled <- (temp$adt_cr - for_icc[which(for_icc$X=="ADT"),"mean"]) / for_icc[which(for_icc$X=="ADT"),"sd"]
colnames(x)[ncol(x)] <- "adt_cr_scaled"


# AIM (aim_tot, 36 total) 
x$aim_tot_scaled <- (x$aim_tot - for_icc[which(for_icc$X=="AIM"),"mean"]) / for_icc[which(for_icc$X=="AIM"),"sd"]
colnames(x)[ncol(x)] <- "aim_tot_scaled"


# CPF (cpf_cr, 40 total) 
x$cpf_cr_scaled <- (x$cpf_cr - for_icc[which(for_icc$X=="CPF"),"mean"]) / for_icc[which(for_icc$X=="CPF"),"sd"]
colnames(x)[ncol(x)] <- "cpf_cr_scaled"


# CPT (cpt_acc, 40 total)  what do I use to scale?
x$cpt_acc_scaled <- NA # for now
# x$cpt_acc_scaled <- (x$cpf_cr - for_icc[which(for_icc$X=="CPT"),"mean"]) / for_icc[which(for_icc$X=="CPT"),"sd"]
# colnames(x)[ncol(x)] <- "cpt_acc_scaled"


# CPW (cpw_cr, 40 total)  
x$cpw_cr_scaled <- (x$cpw_cr - for_icc[which(for_icc$X=="CPW"),"mean"]) / for_icc[which(for_icc$X=="CPW"),"sd"]
colnames(x)[ncol(x)] <- "cpw_cr_scaled"


# DDISC (ddisc_sum, 34 total) 
x$ddisc_sum_scaled <- (x$ddisc_sum - for_icc[which(for_icc$X=="DDISC"),"mean"]) / for_icc[which(for_icc$X=="DDISC"),"sd"]
colnames(x)[ncol(x)] <- "ddisc_sum_scaled"


# DIGSYM (dscor)  
x$dscor_scaled <- (x$dscor - for_icc[which(for_icc$X=="DIGSYM"),"mean"]) / for_icc[which(for_icc$X=="DIGSYM"),"sd"]
colnames(x)[ncol(x)] <- "dscor_scaled"


# EDISC (ddisc_sum, 34 total)  
x$edisc_sum_scaled <- (x$edisc_sum - for_icc[which(for_icc$X=="EDISC"),"mean"]) / for_icc[which(for_icc$X=="EDISC"),"sd"]
colnames(x)[ncol(x)] <- "edisc_sum_scaled"


# ER40 (er40_cr, 40 total)  
x$er40_cr_scaled <- (x$er40_cr - for_icc[which(for_icc$X=="ER40"),"mean"]) / for_icc[which(for_icc$X=="ER40"),"sd"]
colnames(x)[ncol(x)] <- "er40_cr_scaled"


# GNG (gng_cr, 150 total)  
x$gng_cr_scaled <- (x$gng_cr - for_icc[which(for_icc$X=="GNG"),"mean"]) / for_icc[which(for_icc$X=="GNG"),"sd"]
colnames(x)[ncol(x)] <- "gng_cr_scaled"


# MEDF (medf_pc, 36 total) 
temp <- x %>% dplyr::select(medf_pc)
temp$medf_cr <- temp*36*.01
x$medf_pc_scaled <- (temp$medf_cr - for_icc[which(for_icc$X=="MEDF"),"mean"]) / for_icc[which(for_icc$X=="MEDF"),"sd"]
colnames(x)[ncol(x)] <- "medf_pc_scaled"


# PLOT (plot_pc, 15 total)  
temp <- x %>% dplyr::select(plot_pc)
temp$plot_cr <- round(temp*15*.01,2)
x$plot_pc_scaled <- (temp$plot_cr - for_icc[which(for_icc$X=="PLOT"),"mean"]) / for_icc[which(for_icc$X=="PLOT"),"sd"]
colnames(x)[ncol(x)] <- "plot_pc_scaled"


# PMAT (pmat_pc, 24 total)  
temp <- x %>% dplyr::select(pmat_pc)
temp$pmat_cr <- temp*24*.01
x$pmat_pc_scaled <- (temp$pmat_cr - for_icc[which(for_icc$X=="PMAT"),"mean"]) / for_icc[which(for_icc$X=="PMAT"),"sd"]
colnames(x)[ncol(x)] <- "pmat_pc_scaled"


# PRA (pra_acc, _ total)  missing
x$pra_acc_scaled <- NA # for now
# x$pra_acc_scaled <- (x$pra_acc - for_icc[which(for_icc$X=="PRA"),"mean"]) / for_icc[which(for_icc$X=="PRA"),"sd"]
# colnames(x)[ncol(x)] <- "pra_acc_scaled"


# PVRT (ddisc_sum, 8 total)  
x$pvrt_cr_scaled <- (x$pvrt_cr - for_icc[which(for_icc$X=="PVRT"),"mean"]) / for_icc[which(for_icc$X=="PVRT"),"sd"]
colnames(x)[ncol(x)] <- "pvrt_cr_scaled"


# RDISC (rdisc_sum, 40 total) 
x$rdisc_sum_scaled <- (x$rdisc_sum - for_icc[which(for_icc$X=="RDISC"),"mean"]) / for_icc[which(for_icc$X=="RDISC"),"sd"]
colnames(x)[ncol(x)] <- "rdisc_sum_scaled"


# VOLT (volt_cr, 20 total)  
x$volt_cr_scaled <- (x$volt_cr - for_icc[which(for_icc$X=="SVOLT"),"mean"]) / for_icc[which(for_icc$X=="SVOLT"),"sd"]
colnames(x)[ncol(x)] <- "volt_cr_scaled"


x_for_ICC <- x %>% mutate(adt_cr_scaled = adt_cr_scaled$adt_pc,
                          medf_pc_scaled = medf_pc_scaled$medf_pc,
                          plot_pc_scaled = plot_pc_scaled$plot_pc,
                          pmat_pc_scaled = pmat_pc_scaled$pmat_pc) %>% 
  dplyr::select(bblid:medf_same_SMVE_Oreg,adt_cr_scaled,aim_tot_scaled:gng_cr_scaled,medf_pc_scaled:pmat_pc_scaled,pra_acc_scaled:volt_cr_scaled)



# * order regressing scaled full CNB scores ----

# use mean of the scaled variable, add this back in before regressing order out
new_scaled <- data.frame(matrix(rep(NA),nrow = nrow(x_for_ICC),ncol = 17))

for (i in 1:17) {
  new_scaled[,i] <- x_for_ICC[,i+114] + mean(x_for_ICC[,i+114],na.rm  = T)   # from er40_cat_Oreg
}

names(new_scaled) <- paste0(names(x_for_ICC)[115:131],"2")
x_for_ICC <- cbind(x_for_ICC,new_scaled)

sc1 <- matrix(NA,nrow(x_for_ICC),17)

for (i in 1:17) {
  mod <- lm(x_for_ICC[,(i+99)]~proto_3,data=x_for_ICC,na.action=na.exclude)
  sc1[,i] <- scale(residuals(mod,na.action=na.exclude))
}
# i will eventually want to use the code above, but for now, I have to leave out the CAT forms that don't have real scores 

for (i in c(1:3,5:13,15:17)) {
  mod <- lm(x_for_ICC[,(i+131)]~proto_3,data=x_for_ICC,na.action=na.exclude)
  sc1[,i] <- scale(residuals(mod,na.action=na.exclude))
}

colnames(sc1) <- paste0(colnames(x_for_ICC[,115:131]),"_Oreg")

x <- data.frame(x_for_ICC,sc1)
x_TD <- data.frame(x_for_ICC,sc1) %>% filter(study_group == "Healthy Controls")
x_PS <- data.frame(x_for_ICC,sc1) %>% filter(study_group == "Psychosis")
x_MD <- data.frame(x_for_ICC,sc1) %>% filter(study_group == "Mood-Anx-BP")


# * calculating the ICCs ----

# ADT
dat_icc <- x %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")
icc_adt <- icc(dat_icc,type="agreement",model="twoway")$value

# TD only
dat_icc <- x_TD %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# AIM
# dat_icc <- x %>% dplyr::select(aim_tot_scaled_Oreg,aim_cat) # AIM
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")

# TD only
# dat_icc <- x_TD %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")

# PS only
# dat_icc <- x_PS %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")

# MD only
# dat_icc <- x_MD %>% dplyr::select(adt_cr_scaled_Oreg,adt_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")


# CPF
dat_icc <- x %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# CPF v2 -- still need to fix  this
# dat_icc <- x %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# # TD only
# dat_icc <- x_TD %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# # PS only
# dat_icc <- x_PS %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# # MD only
# dat_icc <- x_MD %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")


# CPT
# dat_icc <- x %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat) # CPT
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")

# TD only
# dat_icc <- x_TD %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# PS only
# dat_icc <- x_PS %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# MD only
# dat_icc <- x_MD %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")


# CPW
dat_icc <- x %>% dplyr::select(cpw_cr_scaled_Oreg,cpw_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(cpw_cr_scaled_Oreg,cpw_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(cpw_cr_scaled_Oreg,cpw_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(cpw_cr_scaled_Oreg,cpw_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# DDISC
dat_icc <- x %>% dplyr::select(ddisc_sum_scaled_Oreg,ddisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(ddisc_sum_scaled_Oreg,ddisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(ddisc_sum_scaled_Oreg,ddisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(ddisc_sum_scaled_Oreg,ddisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# DIGSYM
# dat_icc <- x %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat) # DIGSYM
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")


# EDISC
dat_icc <- x %>% dplyr::select(edisc_sum_scaled_Oreg,edisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(edisc_sum_scaled_Oreg,edisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(edisc_sum_scaled_Oreg,edisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(edisc_sum_scaled_Oreg,edisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# ER40
dat_icc <- x %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# GNG
# dat_icc <- x %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat) # GNG
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")


# MEDF
dat_icc <- x %>% dplyr::select(medf_pc_scaled_Oreg,medf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(medf_pc_scaled_Oreg,medf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(medf_pc_scaled_Oreg,medf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(medf_pc_scaled_Oreg,medf_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# PLOT
dat_icc <- x %>% dplyr::select(plot_pc_scaled_Oreg,plot.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(plot_pc_scaled_Oreg,plot.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(plot_pc_scaled_Oreg,plot.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(plot_pc_scaled_Oreg,plot.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# PMAT
dat_icc <- x %>% dplyr::select(pmat_pc_scaled_Oreg,pmat.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(pmat_pc_scaled_Oreg,pmat.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(pmat_pc_scaled_Oreg,pmat.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(pmat_pc_scaled_Oreg,pmat.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# PRA
# dat_icc <- x %>% dplyr::select(cpf_cr_scaled_Oreg,cpf_cat) # PRA
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")

# TD only
# dat_icc <- x_TD %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# # PS only
# dat_icc <- x_PS %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")
# 
# # MD only
# dat_icc <- x_MD %>% dplyr::select(er40_cr_scaled_Oreg,er40_cat_Oreg)
# dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
# icc(dat_icc,type="agreement",model="twoway")


# PVRT
dat_icc <- x %>% dplyr::select(pvrt_cr_scaled_Oreg,pvrt.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(pvrt_cr_scaled_Oreg,pvrt.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(pvrt_cr_scaled_Oreg,pvrt.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(pvrt_cr_scaled_Oreg,pvrt.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# RDISC
dat_icc <- x %>% dplyr::select(rdisc_sum_scaled_Oreg,rdisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(rdisc_sum_scaled_Oreg,rdisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(rdisc_sum_scaled_Oreg,rdisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(rdisc_sum_scaled_Oreg,rdisc.1.00.cat_default_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")


# VOLT
dat_icc <- x %>% dplyr::select(volt_cr_scaled_Oreg,volt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# TD only
dat_icc <- x_TD %>% dplyr::select(volt_cr_scaled_Oreg,volt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# PS only
dat_icc <- x_PS %>% dplyr::select(volt_cr_scaled_Oreg,volt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")

# MD only
dat_icc <- x_MD %>% dplyr::select(volt_cr_scaled_Oreg,volt_cat_Oreg)
dat_icc <- dat_icc[which(rowSums(is.na(dat_icc))==0),]   # keeping only rows that have scores for both full and CAT
icc(dat_icc,type="agreement",model="twoway")








# for XL comparison (specifically looking at memory tasks) ----

all_cnb2 <- all_cnb %>% mutate(dotest = as.Date(dotest,format = "%Y-%m-%d"),date_pra = as.Date(date_pra,format = "%Y-%m-%d"))

all_cnb2$dotest_num <- as.numeric(ymd(all_cnb2$dotest))
all_cnb2$date_pra_num <- as.numeric(ymd(all_cnb2$date_pra))
all_cnb2$date_diff <- abs(all_cnb2$dotest_num - all_cnb2$date_pra_num)

# quick distribution plot of date_diff
{
  dat <- all_cnb2 %>% mutate(same = 1)
  my_plot <- ggplot(dat, aes(x = same, y = date_diff)) + 
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
    geom_hline(yintercept = 7) +
    coord_cartesian(xlim = c(1.2, NA)) +
    theme_minimal() + labs(title = "Distribution of differences in CAT-full CNB test dates", caption = paste("n =",dim(dat %>% filter(date_diff <= 7))[1], "CAT/full tests <= 7 days apart;", "n =",dim(dat %>% filter(date_diff > 7))[1], "more than a week apart"),
                           x = "", y = "Difference (days)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/CAT-full_diffdays_dist_220822.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
}

# split into <= 7d and > 7d difference, memory tasks only
{
  mem_week <- all_cnb2 %>% filter(date_diff <= 7) %>% dplyr::select(bblid:dotest,volt_genus:cpf_w_rtcr,cpw_genus:cpw_w_rtcr,cpf.1.00.v1.cat_target:cpf.1.00.v1.cat_foil,cpw.1.00.v1.cat_target:cpw.1.00.v1.cat_foil,
                                                                    volt.1.00.v1.cat_targets:volt.1.00.v1.cat_foils,cpf2.1.00.v2.cat_target:cpf2.1.00.v2.cat_foil,cpf2.1.00.v1.cat_target:cpf2.1.00.v1.cat_foil) # n=148 as of 8/22/22
  mem_week$cpf_v2_target <- ifelse(!is.na(mem_week$cpf2.1.00.v1.cat_target),mem_week$cpf2.1.00.v1.cat_target,mem_week$cpf2.1.00.v2.cat_target)
  mem_week$cpf_v2_foil <- ifelse(!is.na(mem_week$cpf2.1.00.v1.cat_foil),mem_week$cpf2.1.00.v1.cat_foil,mem_week$cpf2.1.00.v2.cat_foil)
  
  mem_week <- mem_week %>% mutate(cpf_cat = ((cpf.1.00.v1.cat_target + cpf.1.00.v1.cat_foil)/2)*sqrt(2),
                                  cpf2_cat = ((cpf_v2_target + cpf_v2_foil)/2)*sqrt(2),
                                  cpw_cat = ((cpw.1.00.v1.cat_target + cpw.1.00.v1.cat_foil)/2)*sqrt(2),
                                  volt_cat = ((volt.1.00.v1.cat_targets + volt.1.00.v1.cat_foils)/2)*sqrt(2))
  
  mem_more <- all_cnb2 %>% filter(date_diff > 7) %>% dplyr::select(bblid:dotest,volt_genus:cpf_w_rtcr,cpw_genus:cpw_w_rtcr,cpf.1.00.v1.cat_target:cpf.1.00.v1.cat_foil,cpw.1.00.v1.cat_target:cpw.1.00.v1.cat_foil,
                                                                   volt.1.00.v1.cat_targets:volt.1.00.v1.cat_foils,cpf2.1.00.v2.cat_target:cpf2.1.00.v2.cat_foil,cpf2.1.00.v1.cat_target:cpf2.1.00.v1.cat_foil) # n=101 as of 8/22/22
  mem_more$cpf_v2_target <- ifelse(!is.na(mem_more$cpf2.1.00.v1.cat_target),mem_more$cpf2.1.00.v1.cat_target,mem_more$cpf2.1.00.v2.cat_target)
  mem_more$cpf_v2_foil <- ifelse(!is.na(mem_more$cpf2.1.00.v1.cat_foil),mem_more$cpf2.1.00.v1.cat_foil,mem_more$cpf2.1.00.v2.cat_foil)
  
  mem_more <- mem_more %>% mutate(cpf_cat = ((cpf.1.00.v1.cat_target + cpf.1.00.v1.cat_foil)/2)*sqrt(2),
                                  cpf2_cat = ((cpf_v2_target + cpf_v2_foil)/2)*sqrt(2),
                                  cpw_cat = ((cpw.1.00.v1.cat_target + cpw.1.00.v1.cat_foil)/2)*sqrt(2),
                                  volt_cat = ((volt.1.00.v1.cat_targets + volt.1.00.v1.cat_foils)/2)*sqrt(2))
  
  
  # scatters for testing within a week 
  # pdf("data/outputs/scatters/CAT-full_diffdays_scatters_week_220822.pdf",height=9,width=12)
  # pairs.panels(mem_week %>% dplyr::select(matches("cpf_cr|cpf_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(mem_week %>% dplyr::select(matches("cpf_cr|cpf2_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(mem_week %>% dplyr::select(matches("cpw_cr|cpw_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(mem_week %>% dplyr::select(matches("volt_cr|volt_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
  
  # scatters for testing more than a week apart
  # pdf("data/outputs/scatters/CAT-full_diffdays_scatters_overweek_220822.pdf",height=9,width=12)
  # pairs.panels(mem_more %>% dplyr::select(matches("cpf_cr|cpf_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(mem_more %>% dplyr::select(matches("cpf_cr|cpf2_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(mem_more %>% dplyr::select(matches("cpw_cr|cpw_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # pairs.panels(mem_more %>% dplyr::select(matches("volt_cr|volt_cat")),lm=TRUE,scale=TRUE,ci=TRUE)
  # dev.off()
}


# distribution of raw RT for AdaptiveV memory tasks (another possible method of QA)
adaptive_v2 <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v2_20220721_144825.csv")
adaptive_v <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v_20220428_125332.csv")
adaptive_cpfv2_er40v2 <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v_cpfv2_er40v2_20220721_144812.csv")

adaptive_v2 <- rename(adaptive_v2, Dataset.ID = 1) %>% 
  mutate(version = "adaptive_cnb_V2")
adaptive_v <- rename(adaptive_v, Dataset.ID = 1) %>% 
  mutate(version = "adaptive_cnb_V")
adaptive_cpfv2_er40v2 <- rename(adaptive_cpfv2_er40v2, Dataset.ID = 1) %>% 
  mutate(version = "adaptive_cnb_V2_cpf_er40")

adaptive_cpfv2_er40v2 <- adaptive_cpfv2_er40v2 %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v2CE=Dataset.ID,date=`Response Logged Date`,score=`Running Score`) %>% 
  filter(BBLID > 9999) 
adaptive_v <- adaptive_v %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v=Dataset.ID,date=`Response Logged Date`,score=`Running Score`) %>% 
  filter(BBLID > 9999) 
adaptive_v2 <- adaptive_v2 %>% mutate(BBLID=as.numeric(`BBL ID`)) %>% rename(testcode=`Test Code`,decisiontreename=`Decision Tree Collection Name`,datasetid_v2=Dataset.ID,date=`Response Logged Date`,score=`Running Score`) %>% 
  filter(BBLID > 9999) 

{
  # extracting only memory tests from adaptive_v2
  mem_v2 <- adaptive_v2 %>% filter(testcode %in% c("cpf2-1.00-v1-cat","cpw-1.00-v1-cat","volt-1.00-v1-cat"))
  cpf2_v2 <- adaptive_v2 %>% filter(testcode == "cpf2-1.00-v1-cat")
  cpw_v2 <- adaptive_v2 %>% filter(testcode == "cpw-1.00-v1-cat")
  volt_v2 <- adaptive_v2 %>% filter(testcode == "volt-1.00-v1-cat")
  
  # extracting only memory tests from adaptive_v
  mem_v <- adaptive_v %>% filter(testcode %in% c("cpf-1.00-v1-cat","cpw-1.00-v1-cat","volt-1.00-v1-cat"),datasetid_v != 322,BBLID %notin% c(90158,22591))
  cpf_v <- adaptive_v %>% filter(testcode == "cpf-1.00-v1-cat",datasetid_v != 322,BBLID %notin% c(90158,22591))
  cpw_v <- adaptive_v %>% filter(testcode == "cpw-1.00-v1-cat",datasetid_v != 322,BBLID %notin% c(90158,22591))
  volt_v <- adaptive_v %>% filter(testcode == "volt-1.00-v1-cat",datasetid_v != 322,BBLID %notin% c(90158,22591))
  # datasetid 322 is not complete, use 323 only
  # CAT CNB for 90158 was too glitchy (first time aka datasetid 542 got stuck at CPW, second time aka datasetid 543 got stuck at ER40)
  
  # extracting only memory tests from adaptive_cpfv2_er40v2
  cpf2_v <- adaptive_cpfv2_er40v2 %>% filter(testcode == "cpf2-1.00-v2-cat",datasetid_v2CE != 597)
  # datasetid 597 is  an incomplete record
  
  # combine CPW and VOLT since they are the same versions
  cpw_combo <- rbind(cpw_v %>% rename(datasetid = datasetid_v), cpw_v2 %>% rename(datasetid = datasetid_v2))
  volt_combo <- rbind(volt_v %>% rename(datasetid = datasetid_v), volt_v2 %>% rename(datasetid = datasetid_v2))
  
  
  # median RT per pt
  cpf_v <- cpf_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpf_v,.,by="BBLID")
  cpf2_v2 <- cpf2_v2 %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpf2_v2,.,by="BBLID")
  cpf2_v <- cpf2_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpf2_v,.,by="BBLID")
  cpw_combo <- cpw_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpw_combo,.,by="BBLID")
  volt_combo <- volt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(volt_combo,.,by="BBLID")
  
  # save mRT in its own vectors
  cpf_v_mRT <- cpf_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpf2_v2_mRT <- cpf2_v2 %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpf2_v_mRT <- cpf2_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpw_combo_mRT <- cpw_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  volt_combo_mRT <- volt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
}

# plots of RT distribution for CAT CNB memory tasks
{ # CPF v1
  dat <- cpf_v %>% mutate(same=1)
  my_plot <- ggplot(dat, aes(x = same, y = mRT)) + 
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
    theme_minimal() + labs(title = "Distribution of CPFv1 mRT in CAT CNB", caption = paste("n =",length(unique(dat$BBLID))),
                           x = "", y = "Median response time (ms)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/cat_cnb_rt_dist/CAT-CNB_mRT_cpfv1_dist_220825.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
  
  
  
  # CPF v2.1
  dat <- cpf2_v2 %>% mutate(same=1)
  my_plot <- ggplot(dat, aes(x = same, y = mRT)) + 
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
    theme_minimal() + labs(title = "Distribution of CPFv2.1 RT in CAT CNB", caption = paste("n =",length(unique(dat$BBLID))),
                           x = "", y = "Median response time (ms)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/cat_cnb_rt_dist/CAT-CNB_mRT_cpfv2_1_dist_220825.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
  
  
  
  # CPF v2.2
  dat <- cpf2_v %>% mutate(same=1)
  my_plot <- ggplot(dat, aes(x = same, y = mRT)) + 
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
    theme_minimal() + labs(title = "Distribution of CPFv2.2 mRT in CAT CNB", caption = paste("n =",length(unique(dat$BBLID))),
                           x = "", y = "Median response time (ms)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/cat_cnb_rt_dist/CAT-CNB_mRT_cpfv2_2_dist_220825.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
  
  
  
  # CPW
  dat <- cpw_combo %>% mutate(same=1)
  my_plot <- ggplot(dat, aes(x = same, y = mRT)) + 
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
    theme_minimal() + labs(title = "Distribution of CPW mRT in CAT CNB", caption = paste("n =",length(unique(dat$BBLID))),
                           x = "", y = "Median response time (ms)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/cat_cnb_rt_dist/CAT-CNB_mRT_cpw_dist_220825.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
  
  
  
  # VOLT
  dat <- volt_combo %>% mutate(same=1)
  my_plot <- ggplot(dat, aes(x = same, y = mRT)) + 
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
    theme_minimal() + labs(title = "Distribution of VOLT mRT in CAT CNB", caption = paste("n =",length(unique(dat$BBLID))),
                           x = "", y = "Median response time (ms)") + 
    # scale_color_manual(values = wes_palette("Darjeeling1",n=1)) + scale_fill_manual(values = wes_palette("Darjeeling1",n=1)) +
    # scale_y_continuous(breaks = seq(0,175000,25000)) +
    coord_flip() 
  
  # pdf("data/outputs/cat_cnb_rt_dist/CAT-CNB_mRT_volt_dist_220825.pdf",height = 7,width = 10)
  # my_plot
  # dev.off()
}

# table to describe basic stats for RT
{
  #  median, mean, +/-2SD, max, min
  rt_tab <- data.frame(matrix(NA,nrow = 5,ncol = 7))
  rownames(rt_tab) <- c("CPW","VOLT","CPF v1","CPF v2.1","CPF v2.2")
  names(rt_tab) <- c("median","mean","+/- 1SD","max","min","1%","5%")
  
  rt_tab[1,] <- c(median(cpw_combo$mRT,na.rm = T), round(mean(cpw_combo$mRT,na.rm = T),2),
                  paste(as.character(round(mean(cpw_combo$mRT,na.rm = T) - sd(cpw_combo$mRT,na.rm = T),2)),as.character(round(mean(cpw_combo$mRT,na.rm = T) + sd(cpw_combo$mRT,na.rm = T),2)),sep=" , "),
                  max(cpw_combo$mRT,na.rm = T),min(cpw_combo$mRT,na.rm = T),round(quantile(cpw_combo$mRT,0.01,na.rm=TRUE),2),round(quantile(cpw_combo$mRT,0.05,na.rm=TRUE),2))
  rt_tab[2,] <- c(median(volt_combo$mRT,na.rm = T), round(mean(volt_combo$mRT,na.rm = T),2),
                  paste(as.character(round(mean(volt_combo$mRT,na.rm = T) - sd(volt_combo$mRT,na.rm = T),2)),as.character(round(mean(volt_combo$mRT,na.rm = T) + sd(volt_combo$mRT,na.rm = T),2)),sep=" , "),
                  max(volt_combo$mRT,na.rm = T),min(volt_combo$mRT,na.rm = T),round(quantile(volt_combo$mRT,0.01,na.rm=TRUE),2),round(quantile(volt_combo$mRT,0.05,na.rm=TRUE),2))
  rt_tab[3,] <- c(median(cpf_v$mRT,na.rm = T), round(mean(cpf_v$mRT,na.rm = T),2),
                  paste(as.character(round(mean(cpf_v$mRT,na.rm = T) - sd(cpf_v$mRT,na.rm = T),2)),as.character(round(mean(cpf_v$mRT,na.rm = T) + sd(cpf_v$mRT,na.rm = T),2)),sep=" , "),
                  max(cpf_v$mRT,na.rm = T),min(cpf_v$mRT,na.rm = T),round(quantile(cpf_v$mRT,0.01,na.rm=TRUE),2),round(quantile(cpf_v$mRT,0.05,na.rm=TRUE),2))
  rt_tab[4,] <- c(median(cpf2_v2$mRT,na.rm = T), round(mean(cpf2_v2$mRT,na.rm = T),2),
                  paste(as.character(round(mean(cpf2_v2$mRT,na.rm = T) - sd(cpf2_v2$mRT,na.rm = T),2)),as.character(round(mean(cpf2_v2$mRT,na.rm = T) + sd(cpf2_v2$mRT,na.rm = T),2)),sep=" , "),
                  max(cpf2_v2$mRT,na.rm = T),min(cpf2_v2$mRT,na.rm = T),round(quantile(cpf2_v2$mRT,0.01,na.rm=TRUE),2),round(quantile(cpf2_v2$mRT,0.05,na.rm=TRUE),2))
  rt_tab[5,] <- c(median(cpf2_v$mRT,na.rm = T), round(mean(cpf2_v$mRT,na.rm = T),2),
                  paste(as.character(round(mean(cpf2_v$mRT,na.rm = T) - sd(cpf2_v$mRT,na.rm = T),2)),as.character(round(mean(cpf2_v$mRT,na.rm = T) + sd(cpf2_v$mRT,na.rm = T),2)),sep=" , "),
                  max(cpf2_v$mRT,na.rm = T),min(cpf2_v$mRT,na.rm = T),round(quantile(cpf2_v$mRT,0.01,na.rm=TRUE),2),round(quantile(cpf2_v$mRT,0.05,na.rm=TRUE),2))
  
  
  rt_tab %>%
    kbl(caption = "CAT CNB Memory task Response Times (ms)", align = rep("c", 8)) %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    save_kable(file = "data/outputs/cat_cnb_rt_dist/mRT_table_220825.pdf", self_contained = T)
}


# make summary/total scores
{
  cpf_v$uniqtest <- paste(cpf_v$testcode,cpf_v$decisiontreename,sep="_")
  cpf_v_sum <- cpf_v %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpf_v_sum <- left_join(cpf_v_sum,cpf_v_mRT,by=c("BBLID"))
  
  cpf2_v2$uniqtest <- paste(cpf2_v2$testcode,cpf2_v2$decisiontreename,sep="_")
  cpf2_v2_sum <- cpf2_v2 %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v2,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpf2_v2_sum <- left_join(cpf2_v2_sum,cpf2_v2_mRT,by=c("BBLID"))
  
  cpf2_v$uniqtest <- paste(cpf2_v$testcode,cpf2_v$decisiontreename,sep="_")
  cpf2_v_sum <- cpf2_v %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v2CE,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpf2_v_sum <- left_join(cpf2_v_sum,cpf2_v_mRT,by=c("BBLID"))
  
  cpw_combo$uniqtest <- paste(cpw_combo$testcode,cpw_combo$decisiontreename,sep="_")
  cpw_combo_sum <- cpw_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpw_combo_sum <- left_join(cpw_combo_sum,cpw_combo_mRT,by=c("BBLID"))
  
  volt_combo$uniqtest <- paste(volt_combo$testcode,volt_combo$decisiontreename,sep="_")
  volt_combo_sum <- volt_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  volt_combo_sum <- left_join(volt_combo_sum,volt_combo_mRT,by=c("BBLID"))
  
}




# RT QA method ----
# using the same RT QA method as was done for the memory tasks

# separate into tasks
{
  # make uniqtest column before separating
  adaptive_v2$uniqtest <- paste(adaptive_v2$testcode,adaptive_v2$decisiontreename,sep="_")
  adaptive_v$uniqtest <- paste(adaptive_v$testcode,adaptive_v$decisiontreename,sep="_")
  adaptive_cpfv2_er40v2$uniqtest <- paste(adaptive_cpfv2_er40v2$testcode,adaptive_cpfv2_er40v2$decisiontreename,sep="_")
  
  # extracting only memory tests from adaptive_v2
  adt_v2 <- adaptive_v2 %>% filter(testcode == "adt-1.00-cat")
  cpf2_v2 <- adaptive_v2 %>% filter(testcode == "cpf2-1.00-v1-cat")
  cpw_v2 <- adaptive_v2 %>% filter(testcode == "cpw-1.00-v1-cat")
  ddisc_v2 <- adaptive_v2 %>% filter(testcode == "ddisc-1.00-cat")
  edisc_v2 <- adaptive_v2 %>% filter(testcode == "edisc-1.00-cat")
  er40_v2 <- adaptive_v2 %>% filter(testcode == "er40-2.00-cat")
  medf_v2 <- adaptive_v2 %>% filter(testcode == "medf-1.00-cat")
  plot_v2 <- adaptive_v2 %>% filter(testcode == "plot-1.00-cat")
  pmat_v2 <- adaptive_v2 %>% filter(testcode == "pmat-1.00-cat")
  pvrt_v2 <- adaptive_v2 %>% filter(testcode == "pvrt-1.00-cat")
  rdisc_v2 <- adaptive_v2 %>% filter(testcode == "rdisc-1.00-cat")
  volt_v2 <- adaptive_v2 %>% filter(testcode == "volt-1.00-v1-cat")
  
  # extracting only memory tests from adaptive_v
  # 22591 has less rows than others, 
  adt_v <- adaptive_v %>% filter(testcode == "adt-1.00-cat")
  cpf_v <- adaptive_v %>% filter(testcode == "cpf-1.00-v1-cat")
  cpw_v <- adaptive_v %>% filter(testcode == "cpw-1.00-v1-cat")
  ddisc_v <- adaptive_v %>% filter(testcode == "ddisc-1.00-cat")
  edisc_v <- adaptive_v %>% filter(testcode == "edisc-1.00-cat")
  er40_v <- adaptive_v %>% filter(testcode == "er40-1.00-cat")
  medf_v <- adaptive_v %>% filter(testcode == "medf-1.00-cat")
  plot_v <- adaptive_v %>% filter(testcode == "plot-1.00-cat")
  pmat_v <- adaptive_v %>% filter(testcode == "pmat-1.00-cat")
  pvrt_v <- adaptive_v %>% filter(testcode == "pvrt-1.00-cat")
  rdisc_v <- adaptive_v %>% filter(testcode == "rdisc-1.00-cat")
  volt_v <- adaptive_v %>% filter(testcode == "volt-1.00-v1-cat")
  # datasetid 322 is not complete, use 323 only
  # CAT CNB for 90158 was too glitchy (first time aka datasetid 542 got stuck at CPW, second time aka datasetid 543 got stuck at ER40)
  
  # extracting only memory tests from adaptive_cpfv2_er40v2
  cpf2_v <- adaptive_cpfv2_er40v2 %>% filter(testcode == "cpf2-1.00-v2-cat",datasetid_v2CE != 597)
  er402_v <- adaptive_cpfv2_er40v2 %>% filter(testcode == "er40-2.00-cat",datasetid_v2CE != 597)
  # datasetid 597 is  an incomplete record
  
  # combine tests with same versions
  adt_combo <- rbind(adt_v %>% rename(datasetid = datasetid_v), adt_v2 %>% rename(datasetid = datasetid_v2))
  cpw_combo <- rbind(cpw_v %>% rename(datasetid = datasetid_v), cpw_v2 %>% rename(datasetid = datasetid_v2))
  ddisc_combo <- rbind(ddisc_v %>% rename(datasetid = datasetid_v), ddisc_v2 %>% rename(datasetid = datasetid_v2))
  edisc_combo <- rbind(edisc_v %>% rename(datasetid = datasetid_v), edisc_v2 %>% rename(datasetid = datasetid_v2))
  er402_combo <- rbind(er402_v %>% rename(datasetid = datasetid_v2CE), er40_v2 %>% rename(datasetid = datasetid_v2))
  medf_combo <- rbind(medf_v %>% rename(datasetid = datasetid_v), medf_v2 %>% rename(datasetid = datasetid_v2))
  plot_combo <- rbind(plot_v %>% rename(datasetid = datasetid_v), plot_v2 %>% rename(datasetid = datasetid_v2))
  pmat_combo <- rbind(pmat_v %>% rename(datasetid = datasetid_v), pmat_v2 %>% rename(datasetid = datasetid_v2))
  pvrt_combo <- rbind(pvrt_v %>% rename(datasetid = datasetid_v), pvrt_v2 %>% rename(datasetid = datasetid_v2))
  rdisc_combo <- rbind(rdisc_v %>% rename(datasetid = datasetid_v), rdisc_v2 %>% rename(datasetid = datasetid_v2))
  volt_combo <- rbind(volt_v %>% rename(datasetid = datasetid_v), volt_v2 %>% rename(datasetid = datasetid_v2))
  
  
  # median RT per pt
  adt_combo <- adt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(adt_combo,.,by="BBLID")
  cpf_v <- cpf_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpf_v,.,by="BBLID")
  cpf2_v2 <- cpf2_v2 %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpf2_v2,.,by="BBLID")
  cpf2_v <- cpf2_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpf2_v,.,by="BBLID")
  cpw_combo <- cpw_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(cpw_combo,.,by="BBLID")
  ddisc_combo <- ddisc_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(ddisc_combo,.,by="BBLID")
  edisc_combo <- edisc_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(edisc_combo,.,by="BBLID")
  er402_combo <- er402_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(er402_combo,.,by="BBLID")
  er40_v <- er40_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(er40_v,.,by="BBLID")
  medf_combo <- medf_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(medf_combo,.,by="BBLID")
  plot_combo <- plot_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(plot_combo,.,by="BBLID")
  pmat_combo <- pmat_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(pmat_combo,.,by="BBLID")
  pvrt_combo <- pvrt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(pvrt_combo,.,by="BBLID")
  rdisc_combo <- rdisc_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(rdisc_combo,.,by="BBLID")
  volt_combo <- volt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`)) %>% left_join(volt_combo,.,by="BBLID")
  
  # save mRT in its own vectors
  adt_combo_mRT <- adt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpf_v_mRT <- cpf_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpf2_v2_mRT <- cpf2_v2 %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpf2_v_mRT <- cpf2_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  cpw_combo_mRT <- cpw_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  ddisc_combo_mRT <- ddisc_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  edisc_combo_mRT <- edisc_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  er402_combo_mRT <- er402_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  er40_v_mRT <- er40_v %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  medf_combo_mRT <- medf_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  plot_combo_mRT <- plot_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  pmat_combo_mRT <- pmat_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  pvrt_combo_mRT <- pvrt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  rdisc_combo_mRT <- rdisc_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
  volt_combo_mRT <- volt_combo %>% group_by(BBLID) %>% summarise(mRT=median(`Response Time (ms)`))
}

# plots of RT distribution for CAT CNB memory tasks
{ # put this in a loop
  texts <- c("adt_combo","cpf_v","cpf2_v2","cpf2_v","cpw_combo","ddisc_combo","edisc_combo","er40_v",
             "er402_combo","medf_combo","plot_combo","pmat_combo","pvrt_combo","rdisc_combo","volt_combo")
  tests <- mget(paste0(texts,"_mRT"))
  captions <- c("ADT","CPFv1","CPFv2.1","CPFv2.2","CPW","DDISC","EDISC","ER40v1","ER40v2","MEDF","PLOT","PMAT","PVRT","RDISC","VOLT")
  
  mRT_plots <- list()
  
  for (i in 1:length(texts)) {
    test_dat <- tests[[i]]
    dat <- test_dat %>% mutate(same=1)
    mRT_plots[[i]] <- ggplot(dat, aes(x = same, y = mRT)) + 
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
      theme_minimal() + labs(title = paste0("Distribution of ",captions[i]," mRT in CAT CNB"), caption = paste("n =",length(unique(dat$BBLID))),
                             x = "", y = "Median response time (ms)") + coord_flip()
  }
  
  # pdf("data/outputs/cat_cnb_rt_dist/CAT-CNB_mRT_dist_220829.pdf",height = 7,width = 10)
  # for (i in 1:length(texts)){
  #   print(mRT_plots[[i]])
  # }
  # dev.off()
}


# table to describe basic stats for RT
{
  #  median, mean, +/-2SD, max, min
  rt_tab <- data.frame(matrix(NA,nrow = length(texts),ncol = 7))
  rownames(rt_tab) <- captions
  names(rt_tab) <- c("median","mean","+/- 1SD","max","min","1%","5%")
  
  for (i in 1:length(texts)) {
    test_dat <- tests[[i]]
    rt_tab[i,] <- c(median(test_dat$mRT,na.rm = T), round(mean(test_dat$mRT,na.rm = T),2),
                    paste(as.character(round(mean(test_dat$mRT,na.rm = T) - sd(test_dat$mRT,na.rm = T),2)),as.character(round(mean(test_dat$mRT,na.rm = T) + sd(test_dat$mRT,na.rm = T),2)),sep=" , "),
                    max(test_dat$mRT,na.rm = T),min(test_dat$mRT,na.rm = T),round(quantile(test_dat$mRT,0.01,na.rm=TRUE),2),round(quantile(test_dat$mRT,0.05,na.rm=TRUE),2))
  }
  rt_tab %>%
    kbl(caption = "CAT CNB tasks, Response Times (ms)", align = rep("c", 8)) %>%
    kable_classic(full_width = F, html_font = "Cambria") 
    # %>% save_kable(file = "data/outputs/cat_cnb_rt_dist/mRT_all_table_220829.pdf", self_contained = T)
}


# make summary/total scores
{
  adt_combo_sum <- adt_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  adt_combo_sum <- left_join(adt_combo_sum,adt_combo_mRT,by=c("BBLID"))
  
  cpf_v_sum <- cpf_v %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpf_v_sum <- left_join(cpf_v_sum,cpf_v_mRT,by=c("BBLID"))
  
  cpf2_v2_sum <- cpf2_v2 %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v2,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpf2_v2_sum <- left_join(cpf2_v2_sum,cpf2_v2_mRT,by=c("BBLID"))
  
  cpf2_v_sum <- cpf2_v %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v2CE,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpf2_v_sum <- left_join(cpf2_v_sum,cpf2_v_mRT,by=c("BBLID"))
  
  cpw_combo_sum <- cpw_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  cpw_combo_sum <- left_join(cpw_combo_sum,cpw_combo_mRT,by=c("BBLID"))
  
  ddisc_combo_sum <- ddisc_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  ddisc_combo_sum <- left_join(ddisc_combo_sum,ddisc_combo_mRT,by=c("BBLID"))
  
  edisc_combo_sum <- edisc_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  edisc_combo_sum <- left_join(edisc_combo_sum,edisc_combo_mRT,by=c("BBLID"))
  
  er40_v_sum <- er40_v %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid_v,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  er40_v_sum <- left_join(er40_v_sum,er40_v_mRT,by=c("BBLID"))
  
  er402_combo_sum <- er402_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  er402_combo_sum <- left_join(er402_combo_sum,er402_combo_mRT,by=c("BBLID"))
  
  medf_combo_sum <- medf_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  medf_combo_sum <- left_join(medf_combo_sum,medf_combo_mRT,by=c("BBLID"))
  
  plot_combo_sum <- plot_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  plot_combo_sum <- left_join(plot_combo_sum,plot_combo_mRT,by=c("BBLID"))
  
  pmat_combo_sum <- pmat_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  pmat_combo_sum <- left_join(pmat_combo_sum,pmat_combo_mRT,by=c("BBLID"))
  
  pvrt_combo_sum <- pvrt_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  pvrt_combo_sum <- left_join(pvrt_combo_sum,pvrt_combo_mRT,by=c("BBLID"))
  
  rdisc_combo_sum <- rdisc_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  rdisc_combo_sum <- left_join(rdisc_combo_sum,rdisc_combo_mRT,by=c("BBLID"))
  
  volt_combo_sum <- volt_combo %>% filter(`Final Score` == "Yes") %>% dplyr::select(BBLID,datasetid,date,score,uniqtest) %>% 
    pivot_wider(names_from = "uniqtest",values_from = "score") %>% arrange(BBLID)
  volt_combo_sum <- left_join(volt_combo_sum,volt_combo_mRT,by=c("BBLID"))
}


# make table of accuracy scores for people under 1% of mRT
{ # table should have test name, n, n removed, mRT 1% cutoff, max / min score 1, max / min score 2
  acc_tab <- data.frame(matrix(NA,nrow = length(texts),ncol = 7))
  rownames(acc_tab) <- captions
  names(acc_tab) <- c("N","Rows removed","1% cutoff","max T/S/E","min T/S/E","max F/D/N","min F/D/N")
  
  tests_sum <- mget(paste0(texts,"_sum"))
  
  for (i in 1:length(texts)) {
    test_dat <- tests_sum[[i]]
    test <- test_dat %>% filter(mRT <= as.numeric(rt_tab[i,6]))
    
    # first add the max/min T/S/E or default
    acc_tab[i,1:5] <- c(dim(test_dat)[1],dim(test)[1],rt_tab[i,6],round(max(as.numeric(unlist(test[,4]))),3),round(min(as.numeric(unlist(test[,4]))),3))
    
    if (ncol(test_dat) == 6) {
      # add F/D/N for those who have it
      acc_tab[i,6:7] <- c(round(max(as.numeric(unlist(test[,5]))),3),round(min(as.numeric(unlist(test[,5]))),3))
    } else {
      # have NA for the last two columns of the table
      acc_tab[i,6:7] <- rep(NA,2)
    }
  }
  acc_tab %>%
    kbl(caption = "CAT CNB tasks, Response Times (ms)", align = rep("c", 8)) %>%
    kable_classic(full_width = F, html_font = "Cambria") # %>%
    # save_kable(file = "data/outputs/cat_cnb_rt_dist/mRT_all_table_220830.pdf", self_contained = T)
}






