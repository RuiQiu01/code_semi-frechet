######################################################################
# This script plots line charts of AMSE for different methods in the 
# real data analysis on the noisy face dataset, showing how performance 
# varies with noise level or contamination proportion.  
# 
# - export figures for inclusion in the main paper / supplementary material.  
######################################################################


# Load the required library
library(ggplot2)

# Create the data matrix 
data <- as.matrix(read.csv('real data_noise_ratio_70.csv'))
data_matrix <- data.frame(
  Method = rep(c("naive-knn", "semi-knn","naive-nw","semi-nw"), each = 6),
  Samples = rep(c(0, 0.2, 0.4, 0.6, 0.8, 1), times = 4),
  Testing_Error = c(data[,2:5])
)

data_matrix$Method <- factor(data_matrix$Method, levels = c(
  "naive-nw", "naive-knn", 
  "semi-nw", "semi-knn"
))

# plot
plot <- ggplot(data_matrix, aes(x = factor(Samples), y = Testing_Error, color = Method, group = Method)) +
  geom_line(size = 1.2) +  
  geom_point(size = 2.2) +  
  scale_color_manual(values = c("naive-knn" = "#7570b3", "naive-nw" = "#e41a1c",
                                "semi-knn" = "#66a61e", "semi-nw"= "#e6ab02"))+  
  theme_bw(base_size = 20) +
  labs(x = "The proportion of noise contamination", y = "Testing error") +
  #labs(x = "Standard deviation of noise", y = "Testing error") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", size = 0.5, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = c(0.99, 0.5),
    legend.justification = c("right", "center"),
    legend.background = element_rect(fill = NA, color = "grey80"),
    legend.key = element_rect(fill = NA)
  )+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))

print(plot)

#ggsave(plot, file='real data_noise_ratio_70.pdf', width=9, height=6)
