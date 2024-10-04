## cox regression ##
library(survival)
disease_new <- fread("/public/home/lihaiyun/4-GFAP&NFL/data/AD_new.csv") #10 diseases
disease_new$sex <- as.factor(disease_new$sex)
disease_new$race <- as.factor(disease_new$race)

# model1: single protein, no covariates #
#GFAP
model1 <- coxph(Surv(time, status) ~ result_gfap, data=disease_new)
print(model1, digits=3)
model_summary1 <- summary(model1, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary1$coefficients),
  HR = exp(model_summary1$coefficients[, "coef"]),
  se = model_summary1$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary1$conf.int[, "lower .95"],
  HR_upper_95 = model_summary1$conf.int[, "upper .95"],
  p_value = model_summary1$coefficients[, "Pr(>|z|)"])
#NEFL
model1 <- coxph(Surv(time, status) ~ result_nfl, data=disease_new)
print(model1, digits=3)
model_summary1 <- summary(model1, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary1$coefficients),
  HR = exp(model_summary1$coefficients[, "coef"]),
  se = model_summary1$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary1$conf.int[, "lower .95"],
  HR_upper_95 = model_summary1$conf.int[, "upper .95"],
  p_value = model_summary1$coefficients[, "Pr(>|z|)"])

# model2: single protein + age, sex, race #
#GFAP
model2 <- coxph(Surv(time, status) ~ result_gfap + age_baseline + sex + race, data=disease_new)
print(model2, digits=3)
model_summary2 <- summary(model2, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary2$coefficients),
  HR = exp(model_summary2$coefficients[, "coef"]),
  se = model_summary2$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary2$conf.int[, "lower .95"],
  HR_upper_95 = model_summary2$conf.int[, "upper .95"],
  p_value = model_summary2$coefficients[, "Pr(>|z|)"])
#NEFL
model2 <- coxph(Surv(time, status) ~ result_nfl + age_baseline + sex + race, data=disease_new)
print(model2, digits=3)
model_summary2 <- summary(model2, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary2$coefficients),
  HR = exp(model_summary2$coefficients[, "coef"]),
  se = model_summary2$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary2$conf.int[, "lower .95"],
  HR_upper_95 = model_summary2$conf.int[, "upper .95"],
  p_value = model_summary2$coefficients[, "Pr(>|z|)"])

# model3: 2 proteins, no covariates #
model3 <- coxph(Surv(time, status) ~ result_gfap + result_nfl, data=disease_new)
print(model3, digits=3)
model_summary3 <- summary(model3, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary3$coefficients),
  HR = exp(model_summary3$coefficients[, "coef"]),
  se = model_summary3$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary3$conf.int[, "lower .95"],
  HR_upper_95 = model_summary3$conf.int[, "upper .95"],
  p_value = model_summary3$coefficients[, "Pr(>|z|)"])

# model4: 2 proteins, age, sex, race #
model4 <- coxph(Surv(time, status) ~ result_gfap + result_nfl + age_baseline + sex + race, data=disease_new)
print(model4, digits=3)
model_summary4 <- summary(model4, digits=3)
result_table <- data.frame(
  Variable = rownames(model_summary4$coefficients),
  HR = exp(model_summary4$coefficients[, "coef"]),
  se = model_summary4$coefficients[,"se(coef)"],
  HR_lower_95 = model_summary4$conf.int[, "lower .95"],
  HR_upper_95 = model_summary4$conf.int[, "upper .95"],
  p_value = model_summary4$coefficients[, "Pr(>|z|)"])

## ROC, AUC ##
library(pROC)
#model 1, same for model 2,3,4
predict_scores <- predict(model1, type = "lp")
time_point <- 6655
status_at_time <- ifelse(disease_new$time <= time_point & disease_new$status == 1, 1, 0)

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
write.csv(roc_data,"/public/home/lihaiyun/4-GFAP&NFL/result/AUC/GFAP-disease-model1-roc.csv",row.names = F)