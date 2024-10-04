### incorporation of genetic factors ###
## cox regression ##
library(survival)
pheno = "AD" #10 diseases
disease_new <- fread(paste0("/public/home/lihaiyun/4-GFAP&NFL/data/",pheno,"_new.csv"))
disease_new$sex <- as.factor(disease_new$sex)
disease_new$race <- as.factor(disease_new$race)

snp <- fread("/public/home/lihaiyun/4-GFAP&NFL/data/snp_cov_3w.csv")
disease_new <- merge(disease_new,snp, by = "eid",all.x = T)
disease_new[,(ncol(disease_new)-2):ncol(disease_new)] <- lapply(disease_new[,(ncol(disease_new)-2):ncol(disease_new)], as.factor)

## ROC, AUC ##
library(pROC)
last_three_cols <- names(disease_new)[(ncol(disease_new)-2):ncol(disease_new)]
filtered_disease_new <- disease_new %>%  filter(!is.na(.[[last_three_cols[1]]]) &  !is.na(.[[last_three_cols[2]]]) &
                                                  !is.na(.[[last_three_cols[3]]]))
model_new <- coxph(Surv(time, status) ~ result_gfap + result_nfl + age_baseline + sex + race + rs429358_C + rs10895991_C + rs59291670_T, data=disease_new)
print(model_new)
predict_scores <- predict(model_new, type = "lp")

time_point <- 6655
status_at_time <- ifelse(filtered_disease_new$time <= time_point & filtered_disease_new$status == 1, 1, 0)

roc_curve <- roc(status_at_time, predict_scores)
plot(roc_curve)
auc_value <- auc(roc_curve)
print(auc_value)
auc_ci <- ci.auc(roc_curve, conf.level = 0.95)

roc_data <- data.frame(
  False_Positive_Rate = roc_curve$specificities,
  True_Positive_Rate = roc_curve$sensitivities,
  Thresholds = roc_curve$thresholds
)
roc_data$AUC <- auc_value
roc_data$AUC_CI_Lower <- auc_ci[1]
roc_data$AUC_CI_Upper <- auc_ci[3]
write.csv(roc_data,paste0("/public/home/lihaiyun/4-GFAP&NFL/result/AUC/new model-genetic/",pheno,"-model4new-roc.csv"),row.names = F)

# compare with old model4
model4_o <- coxph(Surv(time, status) ~ result_gfap + result_nfl + age_baseline + sex + race, data=filtered_disease_new)
predict_scores_o <- predict(model4_o, type = "lp")
time_point <- 6655
status_at_time <- ifelse(filtered_disease_new$time <= time_point & filtered_disease_new$status == 1, 1, 0)

roc_curve_o <- roc(status_at_time, predict_scores_o)
plot(roc_curve_o)
auc_value_o <- auc(roc_curve_o)
print(auc_value_o)
auc_ci_o <- ci.auc(roc_curve_o, conf.level = 0.95)
roc_data_o <- data.frame(
  False_Positive_Rate = roc_curve_o$specificities,
  True_Positive_Rate = roc_curve_o$sensitivities,
  Thresholds = roc_curve_o$thresholds
)
roc_data_o$AUC <- auc_value_o
roc_data_o$AUC_CI_Lower <- auc_ci_o[1]
roc_data_o$AUC_CI_Upper <- auc_ci_o[3]
write.csv(roc_data_o,paste0("/public/home/lihaiyun/4-GFAP&NFL/result/AUC/new model-genetic/",pheno,"-model4old-roc.csv"),row.names = F)

# delong test
roc_test <- roc.test(roc_curve, roc_curve_o, method = "delong")
print(roc_test)
roc_test_results <- data.frame(
  auc_diff = roc_test$estimate,
  p_value = roc_test$p.value,
  conf_low = roc_test$conf.int[1],
  conf_high = roc_test$conf.int[2]
)
write.csv(roc_test_results,paste0("/public/home/lihaiyun/4-GFAP&NFL/result/AUC/new model-genetic/",pheno,"-DeLong.csv"), row.names = T)

### incorporation of environmental factors ###
## cox regression ##
library(survival)
pheno = "AD" #10 diseases
disease_new <- fread(paste0("/public/home/lihaiyun/4-GFAP&NFL/data/",pheno,"_new.csv"))
disease_new$sex <- as.factor(disease_new$sex)
disease_new$race <- as.factor(disease_new$race)

DM_cov <- data_cmbd[,c(1,8)]
CKD_cov <- data_cmbd[,c(1,10)]
AF_cov <- data_cmbd[,c(1,12)]
smoking_cov <- data_environment[,c(1,11)]
disease_new <- merge(disease_new, DM_cov, by = "eid")
disease_new <- merge(disease_new, CKD_cov, by = "eid")
disease_new <- merge(disease_new, AF_cov, by = "eid")
disease_new <- merge(disease_new,smoking_cov, by = "eid")
str(disease_new)

# model2: age, sex, race : interaction analysis#
model2 <- coxph(Surv(time, status) ~ result_gfap + smoking + result_gfap*smoking + age_baseline + sex + race, data=disease_new) #CKD, DM, smoking
model2 <- coxph(Surv(time, status) ~ result_nfl + AF_baseline + result_nfl*AF_baseline + age_baseline + sex + race, data=disease_new) #AF, CKD, DM
print(model2, digits=3)
model_summary2 <- summary(model2, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary2$coefficients),
  HR = exp(model_summary2$coefficients[, "coef"]),
  se = model_summary2$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary2$conf.int[, "lower .95"],
  HR_upper_95 = model_summary2$conf.int[, "upper .95"],
  p_value = model_summary2$coefficients[, "Pr(>|z|)"])

## ROC, AUC ##
library(pROC)
filtered_disease_new <- subset(disease_new, !is.na(smoking))
model4 <- coxph(Surv(time, status) ~ result_gfap + result_nfl + smoking + DM_baseline + CKD_baseline + age_baseline + sex + race, data=disease_new)
print(model4)
predict_scores <- predict(model4, type = "lp")

time_point <- 6655
status_at_time <- ifelse(filtered_disease_new$time <= time_point & filtered_disease_new$status == 1, 1, 0)

roc_curve <- roc(status_at_time, predict_scores)
plot(roc_curve)
auc_value <- auc(roc_curve)
print(auc_value)
auc_ci <- ci.auc(roc_curve, conf.level = 0.95)

roc_data <- data.frame(
  False_Positive_Rate = roc_curve$specificities,
  True_Positive_Rate = roc_curve$sensitivities,
  Thresholds = roc_curve$thresholds
)
roc_data$AUC <- auc_value
roc_data$AUC_CI_Lower <- auc_ci[1]
roc_data$AUC_CI_Upper <- auc_ci[3]
write.csv(roc_data,paste0("/public/home/lihaiyun/4-GFAP&NFL/result/AUC/new model-environment/",pheno,"-model4new-roc.csv"),row.names = F)

# compare with old model4
model4_o <- coxph(Surv(time, status) ~ result_gfap + result_nfl + age_baseline + sex + race, data=filtered_disease_new)
predict_scores_o <- predict(model4_o, type = "lp")
time_point <- 6655
status_at_time <- ifelse(filtered_disease_new$time <= time_point & filtered_disease_new$status == 1, 1, 0)

roc_curve_o <- roc(status_at_time, predict_scores_o)
plot(roc_curve_o)
auc_value_o <- auc(roc_curve_o)
print(auc_value_o)
auc_ci_o <- ci.auc(roc_curve_o, conf.level = 0.95)
roc_data_o <- data.frame(
  False_Positive_Rate = roc_curve_o$specificities,
  True_Positive_Rate = roc_curve_o$sensitivities,
  Thresholds = roc_curve_o$thresholds
)
roc_data_o$AUC <- auc_value_o
roc_data_o$AUC_CI_Lower <- auc_ci_o[1]
roc_data_o$AUC_CI_Upper <- auc_ci_o[3]
write.csv(roc_data_o,paste0("/public/home/lihaiyun/4-GFAP&NFL/result/AUC/new model-environment/",pheno,"-model4old-roc.csv"),row.names = F)

# delong test
roc_test <- roc.test(roc_curve, roc_curve_o, method = "delong")
print(roc_test)
roc_test_results <- data.frame(
  #auc_model = auc(roc_curve),
  #auc_model_o = auc(roc_curve_o),
  auc_diff = roc_test$estimate,
  p_value = roc_test$p.value,
  conf_low = roc_test$conf.int[1],
  conf_high = roc_test$conf.int[2]
)
write.csv(roc_test_results,paste0("/public/home/lihaiyun/4-GFAP&NFL/result/AUC/new model-environment/",pheno,"-DeLong.csv"), row.names = T)
