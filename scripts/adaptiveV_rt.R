# creating CAT vs Full scatters for RT

# Akira Di Sandro, 8.31.22


`%notin%` <- Negate(`%in%`)
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
library(wesanderson)


# load data ----
{
  adaptive_v2 <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v2_20220829_105100.csv")
  adaptive_v <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v_20220428_125332.csv")
  adaptive_cpfv2_er40v2 <- read_csv("data/inputs/cnb/CNB_CAT_session_adaptive_v_cpfv2_er40v2_20220829_104834.csv")
  # no PRA because we don't need PRA RT
  # adaptive_prad <- read_csv("data/inputs/cnb/CNB_CAT_session_pra-d_20220829_104906.csv")
}

# clean up the adaptive CSVs
{
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
}

# the process
# first, get SD and mean RT from crowdsourced data for each item on each test
# no rapid tests (CPT, GNG) because they're rapid and no AIM/DIGSYM because no itemwise RT available
{
  ADT_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_ADT.csv")
  CPF_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_CPF.csv")
  CPW_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_CPW.csv")
  DDISC_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_DDISC.csv")
  EDISC_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_EDISC.csv")
  ER40_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_ER40.csv")
  MEDF_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_MEDF.csv")
  PLOT_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_PLOT.csv")
  PMAT_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_PMAT.csv")
  PVRT_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_PVRT.csv")
  RDISC_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_RDISC.csv")
  VOLT_rt <- read.csv("data/inputs/crowd_sourced_norms/RT_norms_SVOLT.csv")
}

fullCNB <- read.csv("data/inputs/cnb_merged/cnb_merged_20220906.csv")


# Scaling CAT CNB RT ----
# itemwise RT from adaptive CNB with be scaled using SDs and mean from above (question number is number from itembank)
{
  CAT_ADT <- rbind(adaptive_v %>% filter(testcode == "adt-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "adt-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_CPF <- rbind(# adaptive_v %>% filter(testcode == "cpf-1.00-v1-cat") %>% rename(datasetid = datasetid_v),   # Tyler said not to keep the "wrong" versions of CPF and ER40, leaves n = ~165, 9/6/22
                   adaptive_v2 %>% filter(testcode == "cpf2-1.00-v1-cat") %>% rename(datasetid = datasetid_v2),
                   adaptive_cpfv2_er40v2 %>% filter(testcode == "cpf2-1.00-v2-cat") %>% rename(datasetid = datasetid_v2CE))
  CAT_CPW <- rbind(adaptive_v %>% filter(testcode == "cpw-1.00-v1-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "cpw-1.00-v1-cat") %>% rename(datasetid = datasetid_v2))
  CAT_DDISC <- rbind(adaptive_v %>% filter(testcode == "ddisc-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "ddisc-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_EDISC <- rbind(adaptive_v %>% filter(testcode == "edisc-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "edisc-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_ER40 <- rbind( # adaptive_v %>% filter(testcode == "er40-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "er40-2.00-cat") %>% rename(datasetid = datasetid_v2),
                   adaptive_cpfv2_er40v2 %>% filter(testcode == "er40-2.00-cat") %>% rename(datasetid = datasetid_v2CE))
  CAT_MEDF <- rbind(adaptive_v %>% filter(testcode == "medf-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "medf-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_PLOT <- rbind(adaptive_v %>% filter(testcode == "plot-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "plot-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_PMAT <- rbind(adaptive_v %>% filter(testcode == "pmat-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "pmat-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_PVRT <- rbind(adaptive_v %>% filter(testcode == "pvrt-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "pvrt-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_RDISC <- rbind(adaptive_v %>% filter(testcode == "rdisc-1.00-cat") %>% rename(datasetid = datasetid_v),
                   adaptive_v2 %>% filter(testcode == "rdisc-1.00-cat") %>% rename(datasetid = datasetid_v2))
  CAT_VOLT <- rbind(adaptive_v %>% filter(testcode == "volt-1.00-v1-cat") %>% rename(datasetid = datasetid_v),
                    adaptive_v2 %>% filter(testcode == "volt-1.00-v1-cat") %>% rename(datasetid = datasetid_v2))
}

{ # ADT
  # first check that all items in CAT_ADT have norms in ADT_rt
  CAT_dat <- CAT_ADT
  dat_rt <- ADT_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_ADT <- CAT_dat
}

{ # CPF
  # first check that all items in CAT_CPF have norms in CPF_rt
  CAT_dat <- CAT_CPF
  dat_rt <- CPF_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # F
  
  # need to change item names to match them
  CPF_rt[1:40,"X"] <- gsub("CPF_v1_Item","cpf_a",CPF_rt[1:40,"X"])
  CPF_rt[41:nrow(CPF_rt),"X"] <- gsub("-","_",CPF_rt[41:nrow(CPF_rt),"X"])
  
  dat_rt <- CPF_rt %>% filter(X %in% unique(CAT_dat$Stimulus))
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_CPF <- CAT_dat
}

{ # CPW
  # first check that all items in CAT_CPW have norms in CPW_rt
  CAT_dat <- CAT_CPW
  dat_rt <- CPW_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_CPW <- CAT_dat
}

{ # DDISC
  # first check that all items in CAT_DDISC have norms in DDISC_rt
  CAT_dat <- CAT_DDISC
  dat_rt <- DDISC_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_DDISC <- CAT_dat
}

{ # EDISC
  # first check that all items in CAT_EDISC have norms in EDISC_rt
  CAT_dat <- CAT_EDISC
  dat_rt <- EDISC_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # F
  
  # need to change item names to match them
  EDISC_rt$X <- paste0("E",EDISC_rt$X)
  EDISC_rt[grep("A_",EDISC_rt$X),"X"] <- gsub("A_","A_Item",EDISC_rt[grep("A_",EDISC_rt$X),"X"])
  EDISC_rt[grep("B_",EDISC_rt$X),"X"] <- gsub("B_","B_Item",EDISC_rt[grep("B_",EDISC_rt$X),"X"])
  EDISC_rt[grep("_C_",EDISC_rt$X),"X"] <- gsub("_C_","_C_Item",EDISC_rt[grep("_C_",EDISC_rt$X),"X"])
  
  dat_rt <- EDISC_rt %>% filter(X %in% unique(CAT_dat$Stimulus))
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_EDISC <- CAT_dat
}

{ # ER40
  # first check that all items in CAT_ER40 have norms in ER40_rt
  CAT_dat <- CAT_ER40
  dat_rt <- ER40_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_ER40 <- CAT_dat
}

{ # MEDF
  # first check that all items in CAT_MEDF have norms in MEDF_rt
  CAT_dat <- CAT_MEDF
  dat_rt <- MEDF_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_MEDF <- CAT_dat
}

{ # PLOT
  # first check that all items in CAT_PLOT have norms in PLOT_rt
  CAT_dat <- CAT_PLOT
  dat_rt <- PLOT_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # F
  
  # need to change item names to match them
  PLOT_rt[grep("NewPLOT",PLOT_rt$X),"X"] <- gsub("NewPLOT_v3","new_PLOT",PLOT_rt[grep("NewPLOT",PLOT_rt$X),"X"])
  
  dat_rt <- PLOT_rt %>% filter(X %in% unique(CAT_dat$Stimulus))
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_PLOT <- CAT_dat
}

{ # PMAT
  # first check that all items in CAT_PMAT have norms in PMAT_rt
  CAT_dat <- CAT_PMAT
  dat_rt <- PMAT_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # F
  
  # need to change item names to match them
  PMAT_rt[grep("pmat",PMAT_rt$X),"X"] <- gsub("\\_18.*","",PMAT_rt[grep("pmat",PMAT_rt$X),"X"])
  
  dat_rt <- PMAT_rt %>% filter(X %in% unique(CAT_dat$Stimulus))
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_PMAT <- CAT_dat
}

{ # PVRT
  # first check that all items in CAT_PVRT have norms in PVRT_rt
  CAT_dat <- CAT_PVRT
  dat_rt <- PVRT_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_PVRT <- CAT_dat
}

{ # RDISC
  # first check that all items in CAT_RDISC have norms in RDISC_rt
  CAT_dat <- CAT_RDISC
  dat_rt <- RDISC_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_RDISC <- CAT_dat
}

{ # VOLT
  # first check that all items in CAT_VOLT have norms in VOLT_rt
  CAT_dat <- CAT_VOLT
  dat_rt <- VOLT_rt %>% filter(X %in% unique(CAT_dat$Stimulus)) # only keep items that are used in the CAT version
  all(unique(CAT_dat$Stimulus) %in% unique(dat_rt$X)) # T
  
  # start this as a loop and then maybe make it into a function?
  CAT_dat$scaled <- NA
  
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
  CAT_VOLT <- CAT_dat
}


# Calculate MRT for each pt ----
{ # ADT
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_ADT$BBLID)))
  # calculate MRT
  MRT_ADT <- CAT_ADT %>% group_by(BBLID) %>% summarise(ADT_MRT = median(scaled))
  # merge MRT 
  CAT_ADT <- left_join(CAT_ADT,MRT_ADT, by="BBLID")
  
  # CPF
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_CPF$BBLID)))
  # duplicate data for 15507 (get rid of datasetid 893), 22723 keep original data from when they completed the whole battery (get rid of datasetid 556), 22293 (get rid of datasetid 597 for incomplete data)
  CAT_CPF <- CAT_CPF %>% filter(datasetid %notin% c(893,556,597))
  # calculate MRT
  MRT_CPF <- CAT_CPF %>% group_by(BBLID) %>% summarise(CPF_MRT = median(scaled))
  # merge MRT 
  CAT_CPF <- left_join(CAT_CPF,MRT_CPF, by="BBLID")
  
  # CPW
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_CPW$BBLID)))
  # duplicate data for 22454 (get rid of datasetid 322), 90158 test got stuck at the "smile" screen (get rid of datasetid 542)
  CAT_CPW <- CAT_CPW %>% filter(datasetid %notin% c(322,542))
  # calculate MRT
  MRT_CPW <- CAT_CPW %>% group_by(BBLID) %>% summarise(CPW_MRT = median(scaled))
  # merge MRT 
  CAT_CPW <- left_join(CAT_CPW,MRT_CPW, by="BBLID")
  
  # DDISC
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_DDISC$BBLID)))
  # calculate MRT
  MRT_DDISC <- CAT_DDISC %>% group_by(BBLID) %>% summarise(DDISC_MRT = median(scaled))
  # merge MRT 
  CAT_DDISC <- left_join(CAT_DDISC,MRT_DDISC, by="BBLID")
  
  # EDISC
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_EDISC$BBLID)))
  # calculate MRT
  MRT_EDISC <- CAT_EDISC %>% group_by(BBLID) %>% summarise(EDISC_MRT = median(scaled))
  # merge MRT 
  CAT_EDISC <- left_join(CAT_EDISC,MRT_EDISC, by="BBLID")
  
  # ER40
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_ER40$BBLID)))
  # duplicate data for 15507 (get rid of datasetid 893), 22723 keep original data from when they completed the whole battery (get rid of datasetid 556)
  CAT_ER40 <- CAT_ER40 %>% filter(datasetid %notin% c(893,556))
  # calculate MRT
  MRT_ER40 <- CAT_ER40 %>% group_by(BBLID) %>% summarise(ER40_MRT = median(scaled))
  # merge MRT 
  CAT_ER40 <- left_join(CAT_ER40,MRT_ER40, by="BBLID")
  
  # MEDF
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_MEDF$BBLID)))
  # duplicate data for 22454 (get rid of datasetid 322), 90158 (get rid of datasetid 543)
  CAT_MEDF <- CAT_MEDF %>% filter(datasetid %notin% c(322,543))
  # calculate MRT
  MRT_MEDF <- CAT_MEDF %>% group_by(BBLID) %>% summarise(MEDF_MRT = median(scaled))
  # merge MRT 
  CAT_MEDF <- left_join(CAT_MEDF,MRT_MEDF, by="BBLID")
  
  # PLOT
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_PLOT$BBLID)))
  # calculate MRT
  MRT_PLOT <- CAT_PLOT %>% group_by(BBLID) %>% summarise(PLOT_MRT = median(scaled))
  # merge MRT 
  CAT_PLOT <- left_join(CAT_PLOT,MRT_PLOT, by="BBLID")
  
  # PMAT
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_PMAT$BBLID)))
  # calculate MRT
  MRT_PMAT <- CAT_PMAT %>% group_by(BBLID) %>% summarise(PMAT_MRT = median(scaled))
  # merge MRT 
  CAT_PMAT <- left_join(CAT_PMAT,MRT_PMAT, by="BBLID")
  
  # PVRT
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_PVRT$BBLID)))
  # duplicate data for 22454 (get rid of datasetid 322), 90158 (get rid of datasetid 543)
  CAT_PVRT <- CAT_PVRT %>% filter(datasetid %notin% c(322,543))
  # calculate MRT
  MRT_PVRT <- CAT_PVRT %>% group_by(BBLID) %>% summarise(PVRT_MRT = median(scaled))
  # merge MRT 
  CAT_PVRT <- left_join(CAT_PVRT,MRT_PVRT, by="BBLID")
  
  # RDISC
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_RDISC$BBLID)))
  # calculate MRT
  MRT_RDISC <- CAT_RDISC %>% group_by(BBLID) %>% summarise(RDISC_MRT = median(scaled))
  # merge MRT 
  CAT_RDISC <- left_join(CAT_RDISC,MRT_RDISC, by="BBLID")
  
  # VOLT
  # check that each BBLID only has one complete set of data
  # View(data.frame(table(CAT_VOLT$BBLID)))
  # calculate MRT
  MRT_VOLT <- CAT_VOLT %>% group_by(BBLID) %>% summarise(VOLT_MRT = median(scaled))
  # merge MRT 
  CAT_VOLT <- left_join(CAT_VOLT,MRT_VOLT, by="BBLID")
}

# * try no-scaled median RT ----
{
  noscale <- CAT_ADT %>% group_by(BBLID) %>% summarise(noscale = median(`Response Time (ms)`))
  CAT_ADT <- left_join(CAT_ADT,noscale,by="BBLID")
}


# Merge full and CAT CNB ----
{ # merge with full data
  CNB_merged <- left_join(fullCNB,MRT_ADT,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_CPF,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_CPW,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_DDISC,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_EDISC,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_ER40,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_MEDF,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_PLOT,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_PMAT,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_PVRT,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_RDISC,by=c("bblid"="BBLID"))
  CNB_merged <- left_join(CNB_merged,MRT_VOLT,by=c("bblid"="BBLID"))
  
  # add noscale
  CNB_merged <- left_join(CNB_merged,noscale,by=c("bblid"="BBLID"))
}


# QA ----
# SMVE (not using this method for now)


# RT QA method


# multivariate outlier removal


# winsorize (trim=0.025) ?



# Plot scatters ----

demos <- CNB_merged %>% dplyr::select(bblid:dotest)
merged_RTs <- CNB_merged %>% dplyr::select(matches("tprt|mcr|rtcr|MRT|corrt|mcrrt|totrt|CRRT|noscale")) %>% cbind(demos,.)

# order regress
{
  sc <- matrix(NA,nrow(merged_RTs),ncol(merged_RTs)-ncol(demos))
  
  for (i in 1:(ncol(merged_RTs)-ncol(demos))) {
    mod <- lm(merged_RTs[,(i+ncol(demos))]~proto_3,data=merged_RTs,na.action=na.exclude)
    sc[,i] <- scale(residuals(mod,na.action=na.exclude))
  }
  
  colnames(sc) <- paste0(colnames(merged_RTs[,(ncol(demos)+1):ncol(merged_RTs)]),"_Oreg")
}


# * checking effect size ----
{
  my_gg <- ggplot(merged_RTs, aes(x = proto_3, y = er40_rtcr,color=proto_3)) + 
    ggdist::stat_halfeye(
      adjust = .5, 
      width = .6,
      justification = -.2, 
      .width = 0, 
      point_colour = NA,
      alpha = 0.8
    ) + 
    geom_boxplot(
      width = .12, 
      outlier.color = NA, ## `outlier.shape = NA` works as well
      alpha = 0.5
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
      )
    )  + 
    geom_hline(yintercept = 7) +
    coord_cartesian(xlim = c(1.2, NA)) +
    theme_minimal() + labs(title = "Distribution of RT in full ER40",x = "", y = "RT (ms)") + coord_flip() 
}



x_all <- data.frame(merged_RTs,sc)
x_TD <- data.frame(merged_RTs,sc) %>% filter(study_group %in% c("Healthy Controls"))
x_PS <- data.frame(merged_RTs,sc) %>% filter(study_group %in% c("Psychosis"))
x_MD <- data.frame(merged_RTs,sc) %>% filter(study_group %in% c("Mood-Anx-BP"))

# * all individual tests printed out by condition ---- 

# overall
pdf("data/outputs/scatters/CNB-CAT_RT_test-retest_scatter_matrices_ALL_bytest_220908.pdf",height=9,width=12)
pairs.panels(x_all %>% dplyr::select(matches("adt_rtcr_Oreg|ADT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("aim_totrt_Oreg|S_AIM.AIMTOTRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
pairs.panels(x_all %>% dplyr::select(matches("cpf_w_rtcr_Oreg|CPF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                          # new, corrected cpf (cpfv2)
pairs.panels(x_all %>% dplyr::select(matches("cpt_tprt_Oreg|CPT108.CATCPTT_TPRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_all %>% dplyr::select(matches("cpw_w_rtcr_Oreg|CPW_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("ddisc_mcr_Oreg|DDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("dscorrt_Oreg|S_DIGSYM.DSCORRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_all %>% dplyr::select(matches("dsmcrrt_Oreg|S_DIGSYM.DSMCRRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_all %>% dplyr::select(matches("edisc_mcr_Oreg|EDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("er40_rtcr_Oreg|ER40_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("gng_rtcr_Oreg|GNG60.GNG60_RTCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_all %>% dplyr::select(matches("medf_rtcr_Oreg|MEDF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("plot_rtcr_Oreg|PLOT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("pmat_rtcr_Oreg|PMAT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("pvrt_rtcr_Oreg|PVRT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("rdisc_mcr_Oreg|RDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_all %>% dplyr::select(matches("volt_w_rtcr_Oreg|VOLT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# TD
pdf("data/outputs/scatters/CNB-CAT_RT_test-retest_scatter_matrices_TD_bytest_220908.pdf",height=9,width=12)
pairs.panels(x_TD %>% dplyr::select(matches("adt_rtcr_Oreg|ADT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("aim_totrt_Oreg|S_AIM.AIMTOTRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
pairs.panels(x_TD %>% dplyr::select(matches("cpf_w_rtcr_Oreg|CPF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                          # new, corrected cpf (cpfv2)
pairs.panels(x_TD %>% dplyr::select(matches("cpt_tprt_Oreg|CPT108.CATCPTT_TPRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_TD %>% dplyr::select(matches("cpw_w_rtcr_Oreg|CPW_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("ddisc_mcr_Oreg|DDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("dscorrt_Oreg|S_DIGSYM.DSCORRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_TD %>% dplyr::select(matches("dsmcrrt_Oreg|S_DIGSYM.DSMCRRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_TD %>% dplyr::select(matches("edisc_mcr_Oreg|EDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("er40_rtcr_Oreg|ER40_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("gng_rtcr_Oreg|GNG60.GNG60_RTCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_TD %>% dplyr::select(matches("medf_rtcr_Oreg|MEDF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("plot_rtcr_Oreg|PLOT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("pmat_rtcr_Oreg|PMAT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("pvrt_rtcr_Oreg|PVRT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("rdisc_mcr_Oreg|RDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_TD %>% dplyr::select(matches("volt_w_rtcr_Oreg|VOLT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# PS
pdf("data/outputs/scatters/CNB-CAT_RT_test-retest_scatter_matrices_PS_bytest_220908.pdf",height=9,width=12)
pairs.panels(x_PS %>% dplyr::select(matches("adt_rtcr_Oreg|ADT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("aim_totrt_Oreg|S_AIM.AIMTOTRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
pairs.panels(x_PS %>% dplyr::select(matches("cpf_w_rtcr_Oreg|CPF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                          # new, corrected cpf (cpfv2)
pairs.panels(x_PS %>% dplyr::select(matches("cpt_tprt_Oreg|CPT108.CATCPTT_TPRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_PS %>% dplyr::select(matches("cpw_w_rtcr_Oreg|CPW_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("ddisc_mcr_Oreg|DDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("dscorrt_Oreg|S_DIGSYM.DSCORRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_PS %>% dplyr::select(matches("dsmcrrt_Oreg|S_DIGSYM.DSMCRRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_PS %>% dplyr::select(matches("edisc_mcr_Oreg|EDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("er40_rtcr_Oreg|ER40_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("gng_rtcr_Oreg|GNG60.GNG60_RTCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_PS %>% dplyr::select(matches("medf_rtcr_Oreg|MEDF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("plot_rtcr_Oreg|PLOT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("pmat_rtcr_Oreg|PMAT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("pvrt_rtcr_Oreg|PVRT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("rdisc_mcr_Oreg|RDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_PS %>% dplyr::select(matches("volt_w_rtcr_Oreg|VOLT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()


# MD
pdf("data/outputs/scatters/CNB-CAT_RT_test-retest_scatter_matrices_MD_bytest_220908.pdf",height=9,width=12)
pairs.panels(x_MD %>% dplyr::select(matches("adt_rtcr_Oreg|ADT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("aim_totrt_Oreg|S_AIM.AIMTOTRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # aim
pairs.panels(x_MD %>% dplyr::select(matches("cpf_w_rtcr_Oreg|CPF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)                          # new, corrected cpf (cpfv2)
pairs.panels(x_MD %>% dplyr::select(matches("cpt_tprt_Oreg|CPT108.CATCPTT_TPRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # cpt
pairs.panels(x_MD %>% dplyr::select(matches("cpw_w_rtcr_Oreg|CPW_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("ddisc_mcr_Oreg|DDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("dscorrt_Oreg|S_DIGSYM.DSCORRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym
pairs.panels(x_MD %>% dplyr::select(matches("dsmcrrt_Oreg|S_DIGSYM.DSMCRRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # digsym mem
pairs.panels(x_MD %>% dplyr::select(matches("edisc_mcr_Oreg|EDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("er40_rtcr_Oreg|ER40_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("gng_rtcr_Oreg|GNG60.GNG60_RTCR_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE) # gng
pairs.panels(x_MD %>% dplyr::select(matches("medf_rtcr_Oreg|MEDF_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("plot_rtcr_Oreg|PLOT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("pmat_rtcr_Oreg|PMAT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("pvrt_rtcr_Oreg|PVRT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("rdisc_mcr_Oreg|RDISC_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(x_MD %>% dplyr::select(matches("volt_w_rtcr_Oreg|VOLT_MRT_Oreg")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()













