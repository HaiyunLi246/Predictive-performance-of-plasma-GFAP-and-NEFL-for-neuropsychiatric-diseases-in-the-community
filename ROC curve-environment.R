pheno = 'anxiety'
roc_data1 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/new model-environment/",pheno,"-model4old-3.5w-roc.csv"))
roc_data1$Label <- "old"
roc_data2 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/new model-environment/",pheno,"-model4new-roc.csv"))
roc_data2$Label <- "new"
roc_data <- rbind(roc_data1,roc_data2)

library(ggplot2)
ggplot(roc_data, aes(x = 1-False_Positive_Rate, y = True_Positive_Rate, color = Label, linetype = Label)) +
  geom_line(linewidth = 0.5) +  # 画出 ROC 曲线
  #geom_point(data = roc_data, size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("old" = "#1b7c3d", "new"= "#f16c23")) +
  scale_linetype_manual(values = c("new"= "solid", "old" = "solid")) +
  labs(
    title = paste0("ROC Curve of ", pheno),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # 确保 x 轴从零开始
  scale_y_continuous(expand = c(0, 0)) +  # 确保 y 轴从零开始
  annotate("text", x = 0.27, y = 0.08, 
           label = paste0("GFAP+NEFL+demo: AUC=",round(roc_data$AUC[1], 3), 
           "[", round(roc_data$AUC_CI_Lower[1], 3), "-", round(roc_data$AUC_CI_Upper[1], 3), "]\nGFAP+NEFL+demo+environment: AUC=",
           round(roc_data$AUC[nrow(roc_data)], 3), "[", round(roc_data$AUC_CI_Lower[nrow(roc_data)], 3),"-", round(roc_data$AUC_CI_Upper[nrow(roc_data)], 4),"]"),
           hjust = 0, size =4) +
  theme_minimal() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5,0.3),
        axis.line = element_line(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 14))
