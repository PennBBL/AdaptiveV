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
  for (i in 1:nrow(CAT_dat)) {
    # get stimulus name (stim) and RT
    stim <- unlist(CAT_dat[i,"Stimulus"])
    RT <- as.numeric(CAT_dat[i,"Response Time (ms)"])
    
    # get relevant row from dat_RT
    rt_norms <- dat_rt %>% filter(X == stim) %>% dplyr::select(mean:sd)
    
    CAT_dat[i,"scaled"] <- (RT - rt_norms$mean) / rt_norms$sd
  }
  
  #  try nested loop method to see if it goes faster
  for (i in 1:nrow(dat_rt)) {
    # item name from norms (stim) + norms in rt_norms
    stim <- dat_rt[i,1]
    rt_norms <- dat_rt[i,2:3]
    
    # rows from CAT_dat where stim was administered
    to_scale <- CAT_dat %>% filter(Stimulus == stim) %>% dplyr::select(BBLID,Stimulus,`Response Time (ms)`)
    CAT_dat[which(CAT_dat$Stimulus == stim),"scaled"] <- (to_scale$`Response Time (ms)` - rt_norms$mean) / rt_norms$sd
  }
}

# QA
# SMVE (not using this method for now)


# RT QA method


# multivariate outlier removal









