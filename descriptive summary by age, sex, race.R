### 1. age comparison ###
library(dplyr)
library(readr)
library(purrr)
data <- all_protein_merge_filtered2

# define age groups
breaks <- c(39,45,50,55,60,65,70)
labels <- paste(breaks[-length(breaks)], breaks[-1]-1, sep = "-")
labels <- c(labels, paste(tail(breaks, n=1), "+"))

data$age_group <- cut(data$age_baseline, breaks = c(breaks, Inf), right = FALSE, labels = labels)

summary <- data %>%
  group_by(age_group) %>%
  summarise(
    mean_nfl = mean(result_nfl, na.rm = TRUE),
    max_nfl = max(result_nfl, na.rm = TRUE),
    min_nfl = min(result_nfl, na.rm = TRUE)
  )
write.csv(summary, "/public/home/lihaiyun/4-GFAP&NFL/result/description/age_nfl_summary.csv",row.names = F)

summary_male <- data %>% filter(sex == 1) %>% #1 for male, 0 for female
  group_by(age_group) %>%
  summarise(
    mean_nfl = mean(result_gfap, na.rm = TRUE),
    max_nfl = max(result_gfap, na.rm = TRUE),
    min_nfl = min(result_gfap, na.rm = TRUE)
  )
write.csv(summary_male, "/public/home/lihaiyun/4-GFAP&NFL/result/description/age_male_nfl_summary.csv",row.names = F)

### 2. sex comparison ###
t_test_results <- t.test(result_gfap ~ sex, data = data, var.equal = TRUE) #or nfl
# t test results
mean_group0 <- mean(data$result_gfap[data$sex == 0], na.rm = TRUE)
mean_group1 <- mean(data$result_gfap[data$sex == 1], na.rm = TRUE)
p_value <- t_test_results$p.value
t_value <- t_test_results$statistic
freedom <- t_test_results$parameter  # 自由度

results_sex <- data.frame(
  Protein = "GFAP",
  Mean_Group0 = mean_group0,
  Mean_Group1 = mean_group1,
  T_Value = t_value,
  Degrees_of_Freedom = freedom,
  P_Value = p_value)
results_sex2 <- data.frame(
  Protein = "NEFL",
  Mean_Group0 = mean_group0,
  Mean_Group1 = mean_group1,
  T_Value = t_value,
  Degrees_of_Freedom = freedom,
  P_Value = p_value)
results_sex <- rbind(results_sex,results_sex2)
write.csv(results_sex,"/public/home/lihaiyun/4-GFAP&NFL/result/description/sex_gfap&nfl_summary.csv",row.names = F)

### 3. sex by age group ###
agegroup <- c("39-44","45-49","50-54","55-59","60-64","65-69","70 +")
results_sex_by_age <- data.frame(
  Age_Group = character(0),
  Mean_Group0 = numeric(0),       
  Mean_Group1 = numeric(0), 
  T_Value = numeric(0),
  Degrees_of_Freedom = integer(0),
  P_Value = numeric(0))

for (i in agegroup){
  data_sub <- data %>% filter(age_group == i)
  t_test_results <- t.test(result_nfl ~ sex, data = data_sub, var.equal = TRUE)
  # t test results
  mean_group0 <- mean(data_sub$result_nfl[data$sex == 0], na.rm = TRUE)
  mean_group1 <- mean(data_sub$result_nfl[data$sex == 1], na.rm = TRUE)
  p_value <- t_test_results$p.value
  t_value <- t_test_results$statistic
  freedom <- t_test_results$parameter  # 自由度
  
  results_sub <- data.frame(
    Age_group = i,
    Mean_Group0 = mean_group0,
    Mean_Group1 = mean_group1,
    T_Value = t_value,
    Degrees_of_Freedom = freedom,
    P_Value = p_value)
  results_sex_by_age <- rbind(results_sex_by_age,results_sub)
}
write.csv(results_sex_by_age, "/public/home/lihaiyun/4-GFAP&NFL/result/description/sex_nfl_by age.csv",row.names = F)

### 4. race comparison ###
race_summary <- data %>%
  group_by(age_group, race) %>%
  summarise(
    mean_nfl = mean(result_nfl, na.rm = TRUE),
    max_nfl = max(result_nfl, na.rm = TRUE),
    min_nfl = min(result_nfl, na.rm = TRUE)
  ) %>%
  ungroup()
write.csv(race_summary,"/public/home/lihaiyun/4-GFAP&NFL/result/description/race_nfl_summary.csv", row.names = F)

#one-way ANOVA by age groups
data_race <- subset(data, data$age_group == "65-69")
anova_result <- aov(result_nfl ~ race, data = data_race)
summary(anova_result)
anova_summary <- summary(anova_result)

p_value <- anova_summary[[1]][["Pr(>F)"]][1]
anova_output <- data.frame(
  age_group = "70 +",
  p_value = p_value
)

#t-test for significant groups
pairwise_results <- pairwise.t.test(data_race$result_nfl, data_race$race, p.adjust.method = "bonferroni")
print(pairwise_results)
p_values <- pairwise_results$p.value
library(reshape2)
p_values_melt <- melt(p_values, na.rm = TRUE)
colnames(p_values_melt) <- c("group1", "group2", "adjusted_p_value")


### 5. correlation ###
#pearson correlation between protein and age
test_result <- cor.test(data$age_baseline, data$result_nfl, method = "pearson")
results_gfap <- data.frame(
  protein = "GFAP",
  correlation_coefficient = test_result$estimate,
  p_value = test_result$p.value,
  conf_int_lower = test_result$conf.int[1],
  conf_int_upper = test_result$conf.int[2])
results_nfl <- data.frame(
  protein = "NEFL",
  correlation_coefficient = test_result$estimate,
  p_value = test_result$p.value,
  conf_int_lower = test_result$conf.int[1],
  conf_int_upper = test_result$conf.int[2])
results <- rbind(results_gfap,results_nfl)
write.csv(results, "/public/home/lihaiyun/4-GFAP&NFL/result/description/age-protein-pearson correlation.csv", row.names = F)
