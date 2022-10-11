# Full vs CAT GOA comparisons for AdaptiveV

# Akira Di Sandro, 10.05.2022


`%notin%` <- Negate(`%in%`)
# load packages ----






# load data ----
allGOA <- read.csv("data/inputs/goa/GOA_merged_221006.csv")

pdf("data/outputs/GOA_test_retest_221006.pdf",height=9,width=12)
pairs.panels(allGOA %>% dplyr::select(matches("ext")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(allGOA %>% dplyr::select(matches("mood")),lm=TRUE,scale=TRUE,ci=TRUE)
# pairs.panels(allGOA %>% dplyr::select(matches("per")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(allGOA %>% dplyr::select(matches("phob|fear")),lm=TRUE,scale=TRUE,ci=TRUE)
pairs.panels(allGOA %>% dplyr::select(matches("psy")),lm=TRUE,scale=TRUE,ci=TRUE)
dev.off()









# Figuring out which BBLIDs to check for active Psychosis ----
library(stats)

`%notin%` <- Negate(`%in%`)

# read in bbl_study_all to determine which BBLIDs to keep
studyAll <- read.csv("data/inputs/cnb/bbl_study_all_221011.csv", comment.char="#")  # 315 rows of unique bblid as of 9/30/22

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

# read in most recent CAT GOA Psychosis scores
psy <- read.csv("data/inputs/goa/gacat_completed_psychosis_v1_20221003_1431.csv")
psy_clean <- psy %>% mutate(BBLID = as.numeric(BBL.ID)) %>% filter(BBLID %in% studyEnroll2$bblid) %>% 
  dplyr::select(BBLID,Domain:Reconstructed.Response.Path) %>% arrange(BBLID) %>% 
  left_join(.,studyEnroll2 %>% dplyr::select(bblid,study_group),by=c("BBLID"="bblid"))

# top quarter of psy scores
top25 <- quantile(psy_clean$Score,0.75,na.rm = T)

# export BBLIDs, study_group, and PSY scores of top 25% of psy scorers
psy_top25 <- psy_clean %>% filter(Score >= top25) %>% 
  dplyr::select(BBLID,Score,study_group) %>% arrange(desc(Score))

write.csv(psy_top25,"data/outputs/cat_goa_top_psy_scorers.csv",row.names = F)



