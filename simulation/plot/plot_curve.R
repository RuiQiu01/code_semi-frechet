######################################################################
# This code plots line charts of AMSE results for different methods under Setting VI and VII.  
# It is used to visualize and compare the performance of nw, knn, naive-nw, naive-knn, 
# semi-nw, and semi-knn across scenarios. 
#
# - export figures for inclusion in the supplementary material. 
###########################################################################


# Load the required library
library(ggplot2)

# Create the data matrix 
data <- as.matrix(read.csv('spd-100-8000-1000-Curve-snr_2.csv'))
data_matrix <- data.frame(
  Method = rep(c("knn", "nw", "naive-knn", "naive-nw", "semi-knn", "semi-nw"), each = 12),
  Samples = rep(c(0, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000), times = 6),
  Testing_Error = c(data[,2:7])
)

data_matrix$Method <- factor(data_matrix$Method, levels = c(
  "naive-nw", "naive-knn", 
  "semi-nw", "semi-knn", 
  "nw", "knn"
))

# plot
plot <- ggplot(data_matrix, aes(x = factor(Samples), y = Testing_Error, color = Method, group = Method)) +
  #geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = Method), alpha = 0.3, color = NA, show.legend = FALSE) +  #confidence band
  geom_line(size = 1.2) +  
  geom_point(size = 2.2) +  
  #geom_text(aes(label = round(Testing_Error, 3)), vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("knn" = "#377eb8", "nw" = "#8c510a",  
                                "naive-knn" = "#7570b3", "naive-nw" = "#e41a1c",
                                "semi-knn" = "#66a61e", "semi-nw"= "#e6ab02"))+ 
  theme_bw(base_size = 20) +
  labs(x = "Unlabeled sample size", y = "Testing error") +
  coord_cartesian(ylim = c(0, 0.8))+
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

ggsave(plot, file='spd-100-8000-1000-Curve-snr_2.pdf', width=9, height=6)
