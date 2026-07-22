# 
setwd("E:/PhD_Year1-2/Database_comparison/")


# load package
library(dplyr)  # 
library(ggplot2)  # 
library(cowplot)  # 
library(ggpubr)  # 
library(gridExtra)  #
#####relation of observed and predicted values####
#load csv
data_river<- read.csv("NCBI_BMWP_75TH_TOP_200_Importance_OTU.csv", header = TRUE, row.names = 1)
#calculate R2 and RMSE
r2 <- cor(data_river$Observed_BMWP, data_river$Predicted_BMWP)^2
rmse <- sqrt(mean((data_river$Observed_BMWP - data_river$Predicted_BMWP)^2))

cat("R²:", round(r2, 4), "\n")
cat("RMSE:", round(rmse, 4), "\n")

# plot
NCBI_BMWP_25th <- ggscatter(data_river, x = "Observed_BMWP", y = "Predicted_BMWP",
                             size = 3,
                             color = "#2e2e2e",
                             alpha = 0.75,
                             add = "reg.line",
                             add.params = list(
                               color = "#0072FFFF",
                               fill = "#8EC7FFFF",
                               alpha = 0.25,
                               size = 1.6
                             ),
                             conf.int = TRUE) +
  annotate("text", x = 0, y = 190, 
           label = paste0("R² = ", round(r2, 3), 
                          "\nRMSE = ", round(rmse, 3),
                          "\np < 0.001"),  
           size = 5, hjust = 0) +
  xlab("Measured BMWP") +
  ylab("Predicted BMWP") +
  coord_cartesian(xlim = c(0,250), ylim = c(0, 200)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1.2)) +
  theme(aspect.ratio = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

print(NCBI_BMWP_25th)

# save picture
ggsave("BMWP_prediction_plot.png", data_river_plot, width = 8, height = 8, dpi = 300)
ggsave("NCBI_BMWP_75.pdf",NCBI_BMWP_25th, height = 4, width = 4,dpi = 300)

