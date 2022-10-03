#libraries
library(data.table)
library(dplyr)
library(lubridate)
library(readr)
library(gtools)
library(tidyr)

#This is an aedit test

# function for datediff in months
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01"));lt$year*12 + lt$mon }
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }

# Study Enroll
{
  
  studyAll<-read.csv("/data/secure/bit_dwh/oracle/bbl_study_all.csv", comment.char="#")
  
  #pull out grmpy enrolled subjects
  studyEnroll<-studyAll%>%
    filter(PROTOCOL %in% c("834552 - Adaptive-V"))
  studyEnroll$BBLID = as.numeric(studyEnroll$BBLID)
  
  studyEnroll <- studyEnroll%>%
    filter( BBLID>1000, BBLID != 18026)
  
  
  studyEnroll$formatted_doenroll<- as.Date(as.character(studyEnroll$DOENROLL), format = "%d-%b-%y")
  studyEnroll<- studyEnroll%>%
    group_by(BBLID)%>%
    arrange(formatted_doenroll, .by_group = TRUE)
  
  
  studyEnroll$system_Datediff<- difftime(Sys.Date(),studyEnroll$formatted_doenroll, units = "days")
  
  # The following code gives timepoints, but is currently redundant as we would generally use the timepoint provided in study enroll
  studyEnroll<-within(studyEnroll, mydiff <- ave(system_Datediff, BBLID, FUN=function(x) x - x[1]))
  
  studyEnroll<-studyEnroll %>%
    mutate(timep = ifelse(abs(mydiff)>365,2,ifelse(abs(mydiff)==0,1,3)))
  
}

# Demographics

{

  demo<-read.csv("/data/secure/bit_dwh/oracle/subjectvisitsall_v.csv")
  
  demo2<-demo[demo$BBLID %in% studyEnroll$BBLID,]   # makes a subset of all bblids in stud enroll
  demo2<-demo2[which(demo2$PROTOCOL=='834552 - Adaptive-V'),] # further restricts them by protocol
  demo2<-demo2%>%
    #mutate(formatted_dovisit = as.Date(as.character(DOVISIT)))%>%
    group_by(BBLID)%>%
    mutate(rank = rank(DOVISIT))
  
  # Join study enroll and demo2 by bnlid
  demographics_adaptive <- left_join(studyEnroll,demo2, by = "BBLID")
  names(demographics_adaptive)<-tolower(names(demographics_adaptive))
  
  demographics_adaptive$dobirth<- as.Date(as.character(demographics_adaptive$dobirth), format = "%d-%b-%y") # different files have different date formats, this seems to work for most of them
  
  
  # age of enroll calculated using the fucntions declared above
  demographics_adaptive$age_enroll<- (mapply(mondf,demographics_adaptive$dobirth,demographics_adaptive$formatted_doenroll)/12)
  
  # select the columns you need
  demographics_adaptive<-demographics_adaptive%>%
    select(bblid,protocol.x, study_coordinator, study_group, study_status, formatted_doenroll, dovisit, dobirth, sex, race, ethnic, age_enroll)
  
  
  
}

{
  adaptive<-read.csv("/data/secure/bbldm/adaptive.csv") 
  adaptive<-adaptive%>%
    select(-c("ORD","REAL_LINK","CONDITION_IDENTIFIER"))%>%
    unite(Order, c("ORDER_NAME", "LINK_NAME"))%>%
    mutate(BBLID = as.numeric(BBLID))%>%
    filter(BBLID>1000)
  #%>%
  # pivot_wider(names_from = Order, values_from = DONE)
  
  names(adaptive)<-tolower(names(adaptive))
  
  Orders_adaptive<-left_join(studyEnroll,adaptive, by="bblid")
  
}

write.csv(demographics_adaptive,"/data/secure/lab/pi_repo/ruben_gur/adaptive_v1/demographics_adaptive.csv", row.names = FALSE)
write.csv(Orders_adaptive,"/data/secure/lab/pi_repo/ruben_gur/adaptive_v1/orders_adaptive.csv", row.names = FALSE)

#combined
demo_orders_adaptive<-left_join(demographics_adaptive,adaptive, by="bblid")
write.csv(demo_orders_adaptive,"/data/secure/lab/pi_repo/ruben_gur/adaptive_v1/demo_orders_adaptive.csv", row.names = FALSE)
