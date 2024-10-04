disease <- fread("/public/home/lihaiyun/4-GFAP&NFL/neuropsychiatric disease_50w.csv")
disease <- disease[,-1]

# earliest date of disease diagnosis #
library(lubridate)
date_columns <- lapply(disease, ymd)
date_df <- as.data.frame(date_columns)
earliest_dates <- apply(date_df, 1, function(x) min(x, na.rm = TRUE))
earliest_dates_df <- data.frame(earliest_date_disease = earliest_dates)

baseline_all <-fread("/public/home/lihaiyun/4-GFAP&NFL/baseline_50w.csv")
baseline_all <- cbind(baseline_all,earliest_dates)
all <- cbind(baseline_all, disease)

# merge with protein #
protein_5w <- fread("/public/home/lihaiyun/4-GFAP&NFL/protein_5w.csv")
protein_5w <- protein_5w[,c(1,4,7)]
all_protein_merge <- merge(protein_5w, all, by = "eid")

# step1: exclude outliers*23 #
all_protein_merge_delete <- all_protein_merge %>% filter(earliest_dates == "1900-01-01" | earliest_dates == "1902-02-02")
all_protein_merge <- all_protein_merge %>% filter(!(eid %in% all_protein_merge_delete$eid))
rm(all_protein_merge_delete)

# step2: exclude baseline cases*9k #
all_protein_merge$start <- ymd(all_protein_merge$start)
all_protein_merge_filtered <- all_protein_merge[is.na(all_protein_merge$earliest_dates) | 
                                                  all_protein_merge$earliest_dates > all_protein_merge$start, ]

# step3: exclude follow-up <0 (in control group)*3k #
control <- all_protein_merge_filtered[is.na(all_protein_merge_filtered$earliest_dates), ]
case <- all_protein_merge_filtered[!is.na(all_protein_merge_filtered$earliest_dates), ]

control$end <- pmax(control$death, control$last_inpatient_record, na.rm = TRUE)
control$start <- as.Date(control$start)
control$end <- as.Date(control$end)
control$follow_up <- as.numeric(difftime(control$end, control$start, units = "days"))
control_filtered <- control[is.na(control$follow_up) | control$follow_up >= 0, ]
control_filtered <- control_filtered[,c(1:10,38,39,11:37)]

case$end <- case$earliest_dates
case$start <- as.Date(case$start)
case$end <- as.Date(case$end)
case$follow_up <- as.numeric(difftime(case$end, case$start, units = "days"))
case <- case[,c(1:10,38,39,11:37)]

all_protein_merge_filtered2 <- rbind(case, control_filtered)

# fill the conclusion #
all_protein_merge_filtered2$end[is.na(all_protein_merge_filtered2$end)] <- as.Date("2024-07-01")
all_protein_merge_filtered2$follow_up <- as.numeric(difftime(all_protein_merge_filtered2$end, all_protein_merge_filtered2$start, units = "days"))

# race #
all_protein_merge_filtered2$race <- all_protein_merge_filtered2$`race-0.0`
all_protein_merge_filtered2$race <- ifelse(all_protein_merge_filtered2$race == -3, apply(all_protein_merge_filtered2[, 7:9], 1, function(x) {
  non_na_pos = x[!is.na(x) & x > 0]
  if (length(non_na_pos) > 0) {
    non_na_pos[1]
  } else {
    NA
  }
}), all_protein_merge_filtered2$race)
all_protein_merge_filtered2$race[all_protein_merge_filtered2$race == -1] <- NA
all_protein_merge_filtered2$race[all_protein_merge_filtered2$race %in% c(1, 1001, 1002, 1003)] <- 1  #white
all_protein_merge_filtered2$race[all_protein_merge_filtered2$race %in% c(2, 2004, 6)] <- 2  #mixed/other
all_protein_merge_filtered2$race[all_protein_merge_filtered2$race %in% c(3, 5, 2003, 3001, 3002, 3003, 3004)] <- 3  #Asian
all_protein_merge_filtered2$race[all_protein_merge_filtered2$race %in% c(4, 2001, 2002, 4001, 4002, 4003)] <- 4  #black
all_protein_merge_filtered2$race[is.na(all_protein_merge_filtered2$race)] <- 6
all_protein_merge_filtered2$race[all_protein_merge_filtered2$race == 6] <- 2
table(all_protein_merge_filtered2$race)
#1-white,2-mixed&other,3-Asian,4-black

final <- all_protein_merge_filtered2[,c(1:5,40,10:12,16:39)]
fwrite(final, "/public/home/lihaiyun/4-GFAP&NFL/final+disease_3.6w.csv", row.names = F)
fwrite(all_protein_merge_filtered2, "/public/home/lihaiyun/4-GFAP&NFL/final-all_protein_merge_filtered2.csv", row.names = F)
