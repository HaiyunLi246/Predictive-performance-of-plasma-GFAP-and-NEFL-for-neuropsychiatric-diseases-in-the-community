library(data.table)
library(dplyr)
othercov <- fread("/public/home/lihaiyun/4-GFAP&NFL/environmental covariates_50w.csv")

### Step1: data preparation ###
data <- othercov %>% filter(othercov$eid %in% all_protein_merge_filtered2$eid) #merge with 36,153
#处理education
calculate_education <- function(row) {
  cols <- row[3:8]
  if (any(cols %in% c(1, 2, 5, 6))) {
    return(1)
  }
  if (all(is.na(cols) | cols == -3)) {
    return(NA)
  }
  return(2)
}
data <- data %>%
  rowwise() %>%
  mutate(education = calculate_education(cur_data()))
data <- data[,c(1,2,22,3:21)]
#处理smoking
data$smoking <- data$`20116-0.0`
data$smoking[data$smoking == -3] <- NA
data$smoking[data$smoking == 2] <- 1
#处理alcohol intake
data$alcohol <- data$`1558-0.0`
data$alcohol[data$alcohol == -3] <- NA
data$alcohol <- ifelse(is.na(data$alcohol), NA, 6 - data$alcohol)
#处理exercise
data$exercise_m <- data$`884-0.0`
data$exercise_m[data$exercise_m == -3 | data$exercise_m == -1] <- NA
data$exercise_v <- data$`904-0.0`
data$exercise_v[data$exercise_v == -3 | data$exercise_v == -1] <- NA
#处理diet
#fresh fruit
data$fruit <- data$`1309-0.0`
data$fruit[data$fruit == -3 | data$fruit == -1] <- NA
data$fruit[data$fruit == -10] <- 0
#vegetable
data$vegetable_cook <- data$`1289-0.0`
data$vegetable_cook[data$vegetable_cook == -3 | data$vegetable_cook == -1] <- NA
data$vegetable_cook[data$vegetable_cook == -10] <- 0

data$vegetable_raw <- data$`1299-0.0`
data$vegetable_raw[data$vegetable_raw == -3 | data$vegetable_raw == -1] <- NA
data$vegetable_raw[data$vegetable_raw == -10] <- 0

data$vegetable <- coalesce(data$vegetable_cook + data$vegetable_raw, data$vegetable_cook, data$vegetable_raw)
#salt
data$salt <- data$`1478-0.0`
data$salt[data$salt == -3] <- NA

#coffee&tea&water
data$tea <- data$`1488-0.0`
data$tea[data$tea == -3 | data$tea == -1] <- NA
data$tea[data$tea == -10] <- 0

data$coffee <- data$`1498-0.0`
data$coffee[data$coffee == -3 | data$coffee == -1] <- NA
data$coffee[data$coffee == -10] <- 0

data$water <- data$`1528-0.0`
data$water[data$water == -3 | data$water == -1] <- NA
data$water[data$water == -10] <- 0

### Step2: data_environment for linear regression ###
data_environment <- data[,c(1,2,3,10,22:26,30,27,31:34)]
protein_environment <- all_protein_merge_filtered2[,c(1:5,40)]
data_environment <- merge(protein_environment,data_environment,by="eid")
data_environment[, c("sex", "race","education","smoking","alcohol","salt")] <- lapply(data_environment[, c("sex", "race","education","smoking","alcohol","salt")], as.factor)
#data_environment[, c("exercise_m", "exercise_v","fruit","tea","coffee","water")] <- lapply(data_environment[, c("sex", "race","education","smoking","alcohol","salt")], as.integer)
str(data_environment)
names(data_environment) <- c("eid","result_gfap", "result_nfl","age_baseline","sex","race",
                             "BMI","education","TDI","Environ_score","smoking","alcohol",
                             "exercise_m", "exercise_v","vegetable","fruit","salt","tea","coffee","water")
fwrite(data_environment,"/public/home/lihaiyun/4-GFAP&NFL/data/environment covariates*14_3.6w.csv", row.names = F)

### Step3: Loop!###
data_environment <- fread("/public/home/lihaiyun/4-GFAP&NFL/data/environment covariates*14_3.6w.csv")
data_environment[, c("sex", "race","education","smoking","alcohol","salt")] <- lapply(data_environment[, c("sex", "race","education","smoking","alcohol","salt")], as.factor)
#data_environment[, c("exercise_m", "exercise_v","fruit","tea","coffee","water")] <- lapply(data_environment[, c("sex", "race","education","smoking","alcohol","salt")], as.integer)
str(data_environment)

#normalization?
#data_environment$BMI <- scale(data_environment$BMI)
#data_environment$TDI <- scale(data_environment$TDI)
#data_environment$Environ_score <- scale(data_environment$Environ_score)
#data_environment$result_gfap <- scale(data_environment$result_gfap)
#data_environment$result_nfl <- scale(data_environment$result_nfl)

proteins <- c("result_gfap", "result_nfl")
pheno_columns <- colnames(data_environment)[7:ncol(data_environment)]
Unilm_single <- data.frame()

for (pheno in pheno_columns) {
  pheno <- pheno
  for (protein in proteins) {
    #FML <- as.formula(paste0(protein, '~', pheno)) #model1
    FML <- as.formula(paste0(protein, '~', paste(pheno, 'age_baseline', 'sex', 'race', sep = "+"))) #model2
    tryCatch({
      fit <- lm(FML, data = data_environment)
      sumx <- summary(fit)
      pvalue <- sumx$coefficients[2, 'Pr(>|t|)']
      Estimated_effect <- sumx$coefficients[2, "Estimate"]
      se <- sumx$coefficients[2, "Std. Error"]
      tvalue <- sumx$coefficients[2, 't value']
      
      unilm <- data.frame(
        'Trait' = pheno,
        'Protein' = protein,
        'sample_size' = nrow(data_environment),
        'Estimated_effect' = Estimated_effect,
        'se' = se,
        't_value' = tvalue,
        'lowerCI' = Estimated_effect - 1.96 * se,
        'upperCI' = Estimated_effect + 1.96 * se,
        'pval_raw' = pvalue,
        #'pval_fdr' = p.adjust(pvalue, method = "BH"),
        'pval_bfi' = ifelse(pvalue*21*2 >= 1, 1, pvalue*21*2),
        'pval_significant' = ifelse(pvalue*21*2 >= 0.05, "", ifelse(pvalue*21*2 >= 0.01, "*", ifelse(pvalue*21*2 >= 0.001, "**", "***")))
      )
      Unilm_single <- rbind(Unilm_single, unilm)
    }, error = function(e) {
      cat(paste("ERROR:", pheno, gsub("[\r\n]", "", e), sep = " "))
      info <- paste("ERROR:", gsub("[\r\n]", "", e), sep = " ")
      unilm_err <- data.frame(
        'Trait' = pheno,
        'Protein' = protein,
        'error' = info
      )
      Unilm_single <- rbind(Unilm_single, unilm_err)
    })
  }
}
write.csv(Unilm_single, "/public/home/lihaiyun/4-GFAP&NFL/result/environment/environment covariates-model2.csv", row.names = FALSE)
