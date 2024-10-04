library(ggplot2)
library(dplyr)
library(tidyr)
data_plot <- fread("/public/home/lihaiyun/4-GFAP&NFL/data/AD_new.csv") #10 diseases
data_plot <- data_plot %>% filter(status == 1) %>% mutate(time_year = -ceiling(time / 365),
                                                          time_year = case_when(
                                                            time_year == -16 ~ -15, 
                                                            time_year == -1 ~ -2, 
                                                            TRUE ~ time_year
                                                          )) %>% mutate(time_to_diagnosis = as.numeric(time_year))

selected_years <- seq(-15, -1, by = 1)

#GFAP
new_data <- data.frame(time_to_diagnosis = selected_years)
loess_model <- loess(result_gfap ~ time_to_diagnosis, data = data_plot)
loess_predictions <- predict(loess_model, newdata = new_data, se = TRUE)
new_data$predicted_protein <- loess_predictions$fit
new_data$se <- loess_predictions$se.fit

new_data <- new_data %>%
  mutate(lower_ci = predicted_protein - 1.96 * se,
         upper_ci = predicted_protein + 1.96 * se,
         group = "gfap")

#NEFL
new_data2 <- data.frame(time_to_diagnosis = selected_years)
loess_model <- loess(result_nfl ~ time_to_diagnosis, data = data_plot)
loess_predictions <- predict(loess_model, newdata = new_data2, se = TRUE)
new_data2$predicted_protein <- loess_predictions$fit
new_data2$se <- loess_predictions$se.fit

new_data2 <- new_data2 %>%
  mutate(lower_ci = predicted_protein - 1.96 * se,
         upper_ci = predicted_protein + 1.96 * se,
         group = "nfl")


trajectory <- rbind(new_data,new_data2)
#Drawing
ggplot(trajectory, aes(x = time_to_diagnosis, y = predicted_protein, color = group)) +
  geom_point(size = 2, position = position_dodge(width = 0.3)) +  
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.5, position = position_dodge(width = 0.3)) +
  geom_line(data = trajectory, aes(x = time_to_diagnosis, y = predicted_protein, group=group), size = 1,
            position = position_dodge(width = 0.3)) +
  #geom_ribbon(data = new_data, aes(x = time_to_diagnosis, ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "red") + #CI带
  scale_color_manual(values = c("gfap" = "#ef4968", "nfl" = "#6c7bbb")) +
  labs(title = "AD",
       x = "Time to Diagnosis (years)",
       y = "Plasma protein (NPX)") +
  scale_x_continuous(breaks = seq(-15, -2, by = 1)) +
  theme_minimal()+
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), 
        axis.line = element_line(colour = "black",linewidth = 0.75),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black",size = 12),
        axis.title = element_text(colour = "black",size = 14),
        plot.title = element_text(colour = "black",size = 14),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed", linewidth = 0.5),
        legend.position = "bottom",
        legend.text = element_text(colour = "black",size = 12),
        legend.title = element_text(colour = "black",size = 12))
