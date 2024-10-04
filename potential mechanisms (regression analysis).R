### imaging linear regression ###
# Step1: image for analysis #
library(data.table)
library(utils)

image <- fread("/public/home/lihaiyun/4-GFAP&NFL/image_50w.csv")
cols_to_keep <- c(1, seq(2, ncol(image), by = 2))
new_image <- image[, ..cols_to_keep]
new_image <- new_image %>% filter(new_image$eid %in% all_protein_merge_filtered2$eid)
cov <- all_protein_merge_filtered2[,c(1:5,40)]
cov <- merge(cov,imagesite,by = "eid") #add image scanning site to covs
data <- merge(cov, new_image, by = "eid") #data for linear regression

# Step2: "data_clean" by category (WMH, FA&MD, cortical) #
data_clean <- data[,c(1:7,63:304)] #cortical
data_clean <- data[,c(1:7,9:62)] #FA&MD
data_clean <- data[,c(1:7,8)] #FA&MD
data_clean <- na.omit(data_clean)
data_clean <- as.data.frame(data_clean)

#transform: as.factor for covariates
data_clean$age_baseline <- as.numeric(data_clean$age_baseline)
data_clean$sex <- as.factor(data_clean$sex)
data_clean$race <- as.factor(data_clean$race)
data_clean$site <- as.factor(data_clean$site)

#normalization: scale
data_clean$result_gfap <- scale(data_clean$result_gfap)
data_clean$result_nfl <- scale(data_clean$result_nfl)
data_clean[, 8:ncol(data_clean)] <- scale(data_clean[, 8:ncol(data_clean)])
colnames(data_clean)[2:3] <- c("result_gfap","result_nfl")

# Step3: Loop for linear regression! #
proteins <- c("result_gfap", "result_nfl")
pheno_columns <- colnames(data_clean)[8:ncol(data_clean)]
Unilm_single <- data.frame()

for (pheno_raw in pheno_columns) {
  pheno <- paste0("`", pheno_raw, "`")
  for (protein in proteins) {
    FML <- as.formula(paste0(pheno, '~', paste(protein, 'age_baseline', 'sex', 'race', 'site', sep = "+")))
    tryCatch({
      fit <- lm(FML, data = data_clean)
      sumx <- summary(fit)
      pvalue <- sumx$coefficients[2, 'Pr(>|t|)']
      Estimated_effect <- sumx$coefficients[2, "Estimate"]
      se <- sumx$coefficients[2, "Std. Error"]
      tvalue <- sumx$coefficients[2, 't value']
    
      unilm <- data.frame(
        'Trait' = pheno,
        'Protein' = protein,
        'sample_size' = nrow(data_clean),
        'Estimated_effect' = Estimated_effect,
        'se' = se,
        't_value' = tvalue,
        'lowerCI' = Estimated_effect - 1.96 * se,
        'upperCI' = Estimated_effect + 1.96 * se,
        'pval_raw' = pvalue,
        #'pval_fdr' = p.adjust(pvalue, method = "BH"),
        'pval_bfi' = ifelse(pvalue*297*2 >= 1, 1, pvalue*297*2),
        'pval_significant' = ifelse(pvalue*297*2 >= 0.05, "", ifelse(pvalue*297*2 >= 0.01, "*", ifelse(pvalue*297*2 >= 0.001, "**", "***")))
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
write.csv(Unilm_single, "/public/home/lihaiyun/4-GFAP&NFL/result/imaging/new_imaging_lm.csv", row.names = FALSE)


### cognitive function regression ###
# Step1 #
cognitive_baseline <- fread("/public/home/lihaiyun/4-GFAP&NFL/cognitive test baseline_50w.csv")
cov <- all_protein_merge_filtered2[,c(1:5,40)]
data_cog <- merge(cov, cognitive_baseline, by = "eid") #data for linear regression
data_cog[,c("sex","race")] <- lapply(data_cog[,c("sex","race")], as.factor)

# Step2: linear regression for field 20023 #
proteins <- c("result_gfap", "result_nfl")
pheno_columns <- colnames(data_cog)[9:9]
Unilm_single <- data.frame()

for (pheno_raw in pheno_columns) {
  pheno <- paste0("`", pheno_raw, "`")
  for (protein in proteins) {
    FML <- as.formula(paste0(pheno, '~', paste(protein, 'age_baseline', 'sex', 'race', sep = "+")))
    tryCatch({
      fit <- lm(FML, data = data_cog)
      sumx <- summary(fit)
      pvalue <- sumx$coefficients[2, 'Pr(>|t|)']
      Estimated_effect <- sumx$coefficients[2, "Estimate"]
      se <- sumx$coefficients[2, "Std. Error"]
      tvalue <- sumx$coefficients[2, 't value']
      
      unilm <- data.frame(
        'Trait' = pheno,
        'Protein' = protein,
        'sample_size' = nrow(data_cog)-sum(length(sumx[["na.action"]])),
        'Estimated_effect' = Estimated_effect,
        'se' = se,
        't_value' = tvalue,
        'lowerCI' = Estimated_effect - 1.96 * se,
        'upperCI' = Estimated_effect + 1.96 * se,
        'pval_raw' = pvalue,
        #'pval_fdr' = p.adjust(pvalue, method = "BH"),
        'pval_bfi' = ifelse(pvalue*5*2 >= 1, 1, pvalue*5*2),
        'pval_significant' = ifelse(pvalue*5*2 >= 0.05, "", ifelse(pvalue*5*2 >= 0.01, "*", ifelse(pvalue*5*2 >= 0.001, "**", "***")))
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
write.csv(Unilm_single, "/public/home/lihaiyun/4-GFAP&NFL/result/cognitive test/20023_reaction time_lm.csv", row.names = FALSE)

# Step3: polr for other fields (catord traits) #
#field 4282
data_cog$`4282-0.0` <- sub(-1,NA,data_cog$`4282-0.0`)
#field 399
data_cog$`399-0.0` <- rowMeans(data_cog[,c('399-0.1', '399-0.2')], na.rm = TRUE)
data_cog$`399-0.0`[is.na(data_cog$`399-0.1`) & is.na(data_cog$`399-0.2`)] <- NA
data_cog[,c("20016-0.0","20018-0.0","4282-0.0","399-0.0")] <- lapply(data_cog[,c("20016-0.0","20018-0.0","4282-0.0","399-0.0")], as.factor)

proteins <- c("result_gfap", "result_nfl")
pheno_columns <- colnames(data_cog)[c(7,8,10,14)]
Unilm_single <- data.frame()

library(MASS)
library(lmtest)
for (pheno_raw in pheno_columns) {
  pheno <- paste0("`", pheno_raw, "`")
  for (protein in proteins) {
    FML <- as.formula(paste0(pheno, '~', paste(protein, 'age_baseline', 'sex', 'race', sep = "+")))
    tryCatch({
      fit <- polr(FML, data=data_cog, Hess=TRUE)
      sumx <- summary(fit)
      ctable <- coef(summary(fit))
      ct <- coeftest(fit)
      pvalue <- ct[1,"Pr(>|t|)"]
      Estimated_effect <- ctable[1, "Value"]
      se <- ctable[1, "Std. Error"]
      tvalue <- ctable[1, "t value"]
      
      unilm <- data.frame(
        'Trait' = pheno,
        'Protein' = protein,
        'sample_size' = nrow(data_cog)-sum(length(sumx[["na.action"]])),
        'Estimated_effect' = Estimated_effect,
        'se' = se,
        't_value' = tvalue,
        'lowerCI' = Estimated_effect - 1.96 * se,
        'upperCI' = Estimated_effect + 1.96 * se,
        'pval_raw' = pvalue,
        #'pval_fdr' = p.adjust(pvalue, method = "BH"),
        'pval_bfi' = ifelse(pvalue*5*2 >= 1, 1, pvalue*5*2),
        'pval_significant' = ifelse(pvalue*5*2 >= 0.05, "", ifelse(pvalue*5*2 >= 0.01, "*", ifelse(pvalue*5*2 >= 0.001, "**", "***")))
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
write.csv(Unilm_single, "/public/home/lihaiyun/4-GFAP&NFL/result/cognitive test/other4_polr.csv", row.names = FALSE)
