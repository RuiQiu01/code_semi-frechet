######################################################################
# This code plots the effect of graph radius on the AMSE of semi-nw and semi-knn methods. It provides a 
# visualization of how performance changes with different choices of graph radius across experimental settings.  
# 
#  - export figures for the supplementary material.  
######################################################################

# Load the required library
library(ggplot2)

# Create the data matrix 
data <- as.matrix(read.csv('spd-200-3000-1000-Swiss2-setting1_snr_2_r_impaction.csv'))
data_matrix <- data.frame(
  Method = rep(c("semi-knn", "semi-nw"), each = 16),
  Samples = rep(seq(0.5,2,0.1), times = 2),
  Testing_Error = c(data[,2:3])
)

data_matrix$Method <- factor(data_matrix$Method, levels = c("semi-nw", "semi-knn"))

# plot
plot <- ggplot(data_matrix, aes(x = factor(Samples), y = Testing_Error, color = Method, group = Method)) +
  geom_line(size = 1.2) +  
  geom_point(size = 2.2) +  
  scale_color_manual(values = c("semi-knn" = "#66a61e", "semi-nw"= "#e6ab02"))+ 
  theme_bw(base_size = 20) +
  labs(x = expression("Graph radius (" * "\u00D7" ~ r[0] * ")"), y = "Testing error") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", size = 0.5, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = c(0.99, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA, color = "grey80"),
    legend.key = element_rect(fill = NA)
  )+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))

print(plot)

#ggsave(plot, file='spd-200-3000-1000-Swiss2-setting1_snr_2_r_impaction.pdf', width=9, height=6)
