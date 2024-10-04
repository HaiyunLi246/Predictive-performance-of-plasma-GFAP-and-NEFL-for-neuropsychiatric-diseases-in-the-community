pheno = 'anxiety'
roc_data1 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/GFAP-model1/GFAP-",pheno,"-model1-roc.csv"))
roc_data1$Label <- "GFAP_model1"
roc_data2 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/NEFL-model1/NEFL-",pheno,"-model1-roc.csv"))
roc_data2$Label <- "NEFL_model1"
roc_data3 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/GFAP-model2/GFAP-",pheno,"-model2-roc.csv"))
roc_data3$Label <- "GFAP_model2"
roc_data4 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/NEFL-model2/NEFL-",pheno,"-model2-roc.csv"))
roc_data4$Label <- "NEFL_model2"
roc_data5 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/combine-model3/Combine-",pheno,"-model3-roc.csv"))
roc_data5$Label <- "model3"
roc_data6 <- fread(paste0("D:/Work/3-GFAP, Nfl, brain health/result/AUC/combine-model4/Combine-",pheno,"-model4-roc.csv"))
roc_data6$Label <- "model4"

roc_data <- rbind(roc_data1,roc_data2,roc_data3,roc_data4,roc_data5,roc_data6)

library(ggplot2)
ggplot(roc_data, aes(x = 1-False_Positive_Rate, y = True_Positive_Rate, color = Label)) +
  geom_line(linewidth = 0.5) +  # 画出 ROC 曲线
  #geom_point(data = roc_data, size = 0.1)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("NEFL_model1" = "#99C7F2", "GFAP_model1"= "#F299C7",
                                "NEFL_model2" = "#1663A9", "GFAP_model2"= "#A71661", 
                                "model3" = "#89BF9C", "model4"= "#1b7c3d")) +
  #scale_linetype_manual(values = c("new"= "solid", "old" = "solid")) +
  labs(
    title = paste0("ROC of ", pheno),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # 确保 x 轴从零开始
  scale_y_continuous(expand = c(0, 0)) +  # 确保 y 轴从零开始
  annotate("text", x = 0.5, y = 0.15, 
           label = paste0("GFAP_model1: AUC=",round(roc_data1$AUC[1], 3), "[", round(roc_data1$AUC_CI_Lower[1], 3), "-", round(roc_data1$AUC_CI_Upper[1], 3), 
                          "]\nNEFL_model1: AUC=",round(roc_data2$AUC[1], 3), "[", round(roc_data2$AUC_CI_Lower[1], 3),"-", round(roc_data2$AUC_CI_Upper[1], 3),
                          "]\nGFAP_model2: AUC=",round(roc_data3$AUC[1], 3), "[", round(roc_data3$AUC_CI_Lower[1], 3),"-", round(roc_data3$AUC_CI_Upper[1], 3),
                          "]\nNEFL_model2: AUC=",round(roc_data4$AUC[1], 3), "[", round(roc_data4$AUC_CI_Lower[1], 3),"-", round(roc_data4$AUC_CI_Upper[1], 3),
                          "]\nModel3: AUC=",round(roc_data5$AUC[1], 3), "[", round(roc_data5$AUC_CI_Lower[1], 3),"-", round(roc_data5$AUC_CI_Upper[1], 3),
                          "]\nModel4: AUC=",round(roc_data6$AUC[1], 3), "[", round(roc_data6$AUC_CI_Lower[1], 3),"-", round(roc_data6$AUC_CI_Upper[1], 3),"]"),
           hjust = 0, size =4) +
  theme_minimal() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5,0.3),
        axis.line = element_line(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 14))
