# DISC task scoring (3 methods)

# Akira Di Sandro, 03.15.22



# load necessary packages ----
library(tidyr)
library(dplyr)
library(irtoys)
library(missForest)
library(catR)



# get only DISC data from CNB dump
cnb_dump <- read.csv("data/inputs/cnb_merged_20220311.csv")
cnb_dump <- cnb_dump[,-7]
names(cnb_dump)[6] <- "bblid"

demos <- cnb_dump %>% 
  filter(test_sessions.siteid == "adaptive_v") %>% 
  filter(bblid > 9999) %>% 
  dplyr::select(matches("bblid") | matches("^test_") | matches("^datasetid") | matches("^platform"))

disc <- cnb_dump %>%
  filter(test_sessions.siteid == "adaptive_v") %>% 
  filter(bblid > 9999) %>% 
  dplyr::select(matches("DISC"))

temp <- cbind(demos,disc) %>% 
  subset(rowSums(is.na(.)) < 100)

demos <- temp %>% 
  dplyr::select(!matches("DISC"))

disc <- temp %>% 
  dplyr::select(matches("DISC"))

ddisc <- disc %>% 
  dplyr::select(matches("DDISC")) %>% 
  cbind(demos,.)

edisc <- disc %>% 
  dplyr::select(matches("EDISC")) %>% 
  cbind(demos,.)

rdisc <- disc %>% 
  dplyr::select(matches("RDISC")) %>% 
  cbind(demos,.)

rdisc <- rdisc[,1:104]
# rdisc[178:181,104:185] <- rdisc[178:181,21:102]


# * imputing data ----
# if > 50% missing data in row, delete row
# if < 50% missing data in row, impute missing values

d_qs <- ddisc %>% 
  dplyr::select(matches("KDDISC.q_"))
all(rowSums(is.na(d_qs))==0)     # no imputing necessary

e_qs <- edisc %>% 
  dplyr::select(matches("_resp"))
all(rowSums(is.na(e_qs))==0)     # no imputing necessary

r_qs <- rdisc %>% 
  dplyr::select(matches("KRDISC.q_"))
all(rowSums(is.na(r_qs))==0)     # some NAs

index <- 1:ncol(r_qs)
r_qs[,index] <- lapply(r_qs[,index], as.factor)
r_qs <- missForest(r_qs)$ximp
r_qs[,index] <- lapply(r_qs[,index], as.numeric)



# Method 1: Sum Scores ----

d_qs <- d_qs-1
ddisc$ddisc_sum <- rowSums(d_qs)
# write.csv(ddisc,"data/outputs/DDISC_sum_17Mar22.csv")

e_qs <- e_qs-1
edisc$edisc_sum <- rowSums(e_qs[,101:ncol(e_qs)])
# write.csv(edisc,"data/outputs/EDISC_sum_17Mar22.csv")

r_qs <- r_qs-1                          # last 4 rows are missing q_**.1 data but excluding this, all of these columns match the previous ones
rdisc$rdisc_sum <- rowSums(r_qs)
write.csv(temp_r,"data/outputs/RDISC_17Mar22.csv")


all_same <- c()
for (i in 21:102){
  all_same <- c(all_same,all(temp_r[,i]==temp_r[,i+83]))
}

all(temp_r$KRDISC.q_09==temp_r$KRDISC.q_09.1)







# Method 2: IRT Scoring ----

# use espEst() with model="GRM"







# Method 3: Dan's Matlab method ----








