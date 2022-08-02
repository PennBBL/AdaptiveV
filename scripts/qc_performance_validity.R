# QC/performance validity for Adaptive V CNB data

# Akira Di Sandro, 5.5.22


# load packages ----



# load data ----

fulldat <- read_csv("data/inputs/cnb_merged_20220613.csv")
itemwise <- read_csv("data/inputs/adaptive_cnb_battery_itemlevel_report_335.csv")

itemwise_short <- read_csv("data/inputs/adaptive_cpt_aim_digsym_gng_itemlevel_20220408.csv")
itemwise_short <- itemwise_short[which(rowSums(!is.na(itemwise_short[,17:nrow(itemwise_short)]))),]



# 1. autovalidation codes ----







# 2. performance validity scores ----

common_col <- intersect(colnames(itemwise)[1:35],colnames(itemwise_short)[1:35])
common_bbl <- intersect(unique(itemwise$BBLID),unique(itemwise_short$BBLID))

itemwise_full <- full_join(itemwise,itemwise_short,by=common_col)
itemwise_full <- full_join(itemwise,itemwise_short,by="BBLID")






# 3. assessor comments ----








