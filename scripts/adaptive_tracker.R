#libraries
library(data.table)
library(dplyr)
library(lubridate)
library(readr)
library(gtools)
library(tidyr)

# function for datediff in months
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01"));lt$year*12 + lt$mon }
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }
`%notin%` <- Negate(`%in%`)

# Study Enroll
{
  
  studyAll<-read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/bbl_study_all.csv", comment.char="#")
 
  #pull out grmpy enrolled subjects
  studyEnroll<-studyAll%>%
    filter(PROTOCOL %in% c("834552 - Adaptive-V"))
  studyEnroll$BBLID = as.numeric(studyEnroll$BBLID)
  
  studyEnroll <- studyEnroll%>%
    filter( BBLID>1000, BBLID != 18026)
  studyEnroll$formatted_doenroll<- as.Date(as.character(studyEnroll$DOENROLL), format = "%d-%b-%y")
  
  names(studyEnroll)<-tolower(names(studyEnroll))
  studyEnroll2<-studyEnroll%>%
    select(bblid, study_status, study_group, timepoints, formatted_doenroll)
}

#Demographics
{
  
  demo<-read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/subjectvisitsall_v.csv")

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
  
  nodemo<-demo2[demo2$bblid %notin%demo3$bblid,]
  
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
  
  demographics_adaptive <- left_join(studyEnroll,demofinal, by = "bblid")
  demographics_adaptive$age_enroll<- (mapply(mondf,as.Date(as.character(demographics_adaptive$dobirth), format = "%d-%b-%y"),demographics_adaptive$formatted_doenroll)/12)
  demographics_adaptive<-demographics_adaptive%>%
    select(bblid, study_status, study_group, sex, race, ethnic, age_enroll)

}
#CNB

{

cnb<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/cnb_merged_webcnp_surveys_allbblprjcts_longform.csv")

cnb1<- cnb%>%
  select(matches("datasetid_platform|test_sessions.datasetid|siteid|bblid.clean|age|battery|dotest|
gng|cpf|medf|pvrt|er40|cpt|cpw|aim|adt|plot|volt|pmat|digsym|dscor|dsmemcr|dscorrt|dsmcrrt|KDDISC.valid_code
|KRDISC.valid_code|EDISC.valid_code"))

colnames(cnb1)[colnames(cnb1) %in% c("test_sessions.datasetid", "test_sessions.siteid",
                                     "test_sessions.bblid.clean", "test_sessions_v.age",            
                                     "test_sessions_v.battery",           
                                     "test_sessions_v.dotest" )] <-
                                     c("datasetid", "siteid",       
                                      "bblid.clean", "age","battery","dotest")

cnb2<-cnb1[cnb1$bblid.clean %in% studyEnroll$bblid,]%>%
  filter(siteid == "adaptive_v")

cnb2_subset<-cnb2%>%
  select(bblid.clean,dotest,battery)%>%
  distinct(bblid.clean,dotest,battery)%>%
  mutate(status = 1)%>%
  pivot_wider(names_from = battery, values_from = status)





missing<-studyEnroll2[studyEnroll2$bblid %notin% cnb2$bblid.clean,]
  

tracker<- left_join(demographics_adaptive,cnb2_subset, by = c("bblid" = "bblid.clean"))

}

#Adaptive CNB
{
adaptive_v2<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/CNB_CAT_session_adaptive_v2_20220428_125604.csv")
adaptive_v<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/CNB_CAT_session_adaptive_v_20220428_125332.csv")
adaptive_cpfv2_er40v2<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/CNB_CAT_session_adaptive_v_cpfv2_er40v2_20220428_125641.csv")
adaptive_prad<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/CNB_CAT_session_pra-d_20220428_125704.csv")

adaptive_v2<-rename(adaptive_v2, Dataset.ID = 1)%>%
  mutate(version = "adaptive_cnb_V2")
adaptive_v<-rename(adaptive_v, Dataset.ID = 1)%>%
  mutate(version = "adaptive_cnb_V")
adaptive_cpfv2_er40v2<-rename(adaptive_cpfv2_er40v2, Dataset.ID = 1)%>%
  mutate(version = "adaptive_cnb_V2_cpf_er40")
adaptive_prad<-rename(adaptive_prad, Dataset.ID = 1)%>%
  mutate(version = "adaptive_cnb_pra_d")

adaptive_cnb<- bind_rows(adaptive_v,adaptive_v2,adaptive_cpfv2_er40v2,adaptive_prad)%>%
  rename(BBLID = BBL.ID)%>%
  filter(!is.na(as.numeric(BBLID)))%>%
  filter(BBLID>1000)

adaptive_cnb_subset<-adaptive_cnb%>%
  select(BBLID, version)%>%
  distinct(BBLID,version)%>%
  mutate(status = 1)%>%
  pivot_wider(names_from = version, values_from = status)

adaptive_cnb_subset$BBLID<- as.numeric(adaptive_cnb_subset$BBLID)
  


tracker2<-left_join(tracker,adaptive_cnb_subset, by = c("bblid" = "BBLID"))

}

#Phenome
{
mood_anxiety<-read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/GACAT_completed_Mood_Anxiety_v1_20220502_1605.csv")%>%
  rename(bblid= 1)
externalizing<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/GACAT_completed_Externalizing_v1_20220428_1723.csv")%>%
  rename(bblid= 1)
phobia<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/GACAT_completed_Phobias_v1_20220428_1713.csv")%>%
  rename(bblid= 1)

phenome<-bind_rows(mood_anxiety,externalizing,phobia)

phenome_subset<-phenome%>%
  select(bblid,Domain,Score)%>%
  filter(!is.na(as.numeric(bblid)))%>%
  filter(bblid>1000)%>%
  pivot_wider(names_from = Domain, values_from = Score)

phenome_subset$bblid<- as.numeric(phenome_subset$bblid)

tracker3<- left_join(tracker2,phenome_subset, by = "bblid")
}

#Oracle

{
  adaptive<-read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/adaptive.csv") 
  adaptive<-adaptive%>%
    select(-c("ORD","REAL_LINK","CONDITION_IDENTIFIER","LIST","GRP"))%>%
    unite(Order, c("ORDER_NAME", "LINK_NAME"))%>%
    mutate(BBLID = as.numeric(BBLID))%>%
    filter(BBLID>1000)
  
  adaptive$Order <- paste("Oracle", adaptive$Order, sep="_")
    
  adaptive_oracle_pivot<- adaptive%>%
    pivot_wider(names_from = Order, values_from = DONE)
  
  names(adaptive_oracle_pivot)<-tolower(names(adaptive_oracle_pivot))
  
  tracker4<- left_join(tracker3,adaptive_oracle_pivot, by = "bblid")
}



#GOASSESS Full
{
  common_interview<- read.csv("C:/Users/Mrugank/Desktop/oracle dump/oracle_dump_update/axis_common_interview.csv")%>%
    select(bblid,
           assessment)%>%
    filter(assessment == 'PNC GOASSESS(GAF Mod)' )%>%
    filter(bblid>1000)%>%
    mutate(status = 1)%>%
    distinct(bblid,assessment, status)%>%
    pivot_wider(names_from = assessment, values_from = status)
  
  

  tracker5<- left_join(tracker4, common_interview, by = "bblid")
  
  missing_goa<-common_interview[common_interview$bblid %notin% studyEnroll2$bblid,]
}

tracker6 = data.frame(lapply(tracker5, as.character), stringsAsFactors=FALSE)
write.csv(tracker6,"C:/Users/Mrugank/Desktop/Adaptive_V/adaptive_tracker1.csv", row.names = FALSE)

