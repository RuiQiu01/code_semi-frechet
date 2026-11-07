######################################################################
# This code generates boxplots for the real data analysis. 
######################################################################

library("ggplot2")
library("ggpubr")
data <- matrix(0,800,3)
data = data.frame(data)
temp1 <- read.table(file="real data_face_70label_sphere.csv", header=TRUE,sep=",")
temp2 <- read.table(file="real data_face_140label_sphere.csv", header=TRUE,sep=",")
data[1:100,1] <- temp1$V1[1:100]
data[101:200,1] <- temp1$V2[1:100]
data[201:300,1] <- temp1$V3[1:100]
data[301:400,1] <- temp1$V4[1:100]
data[401:500,1] <- temp2$V1[1:100]
data[501:600,1] <- temp2$V2[1:100]
data[601:700,1] <- temp2$V3[1:100]
data[701:800,1] <- temp2$V4[1:100]

data[1:100,2] <- "naive-knn"
data[101:200,2] <- "semi-knn"
data[201:300,2] <- "naive-nw"
data[301:400,2] <- "semi-nw"
data[401:500,2] <- "naive-knn"
data[501:600,2] <- "semi-knn"
data[601:700,2] <- "naive-nw"
data[701:800,2] <- "semi-nw"

data[1:400,3] <- 70
data[401:800,3] <- 140
names(data)<-c("MSE","Method","n")
data$Method <- as.factor(data$Method)
data$n <- as.factor(data$n)

plot <- ggboxplot(data, 
          
          x = 'Method', 
         
          y = 'MSE', 
          
          fill = 'n', 
          
          bxp.errorbar = T, 
          
          bxp.errorbar.width = 0.2, 
          
          palette = 'npg', 
          
          add = "none"
          
) +
  
  labs(title = 'Up-down angle', 
       
       x = 'Method', 
       
       y = 'Mean absolute error (deg)' 
       
  )+
  
  theme(
    
    plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
    
    plot.subtitle = element_text(color = 'black', size   = 16,hjust = 0.5),
    
    plot.caption  = element_text(color = 'black', size   = 16,face = 'italic', hjust = 1),
    
    axis.text.x   = element_text(color = 'black', size = 16, angle = 0),
    
    axis.text.y   = element_text(color = 'black', size = 16, angle = 0),
    
    axis.title.x  = element_text(color = 'black', size = 16, angle = 0),
    
    axis.title.y  = element_text(color = 'black', size = 16, angle = 90),
    
    legend.title  = element_text(color = 'black', size  = 16),
    
    legend.text   = element_text(color = 'black', size   = 16),
    
    axis.line.y = element_line(color = 'black', linetype = 'solid'), 
    
    axis.line.x = element_line (color = 'black',linetype = 'solid'), 
    
    panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) 
    
  )+
  # Modify the legend order
  scale_x_discrete(limits = c("naive-knn","semi-knn", "naive-nw", "semi-nw"))

print(plot)

ggsave(plot, file='real data boxplot_sphere_ud.pdf', width=8, height=4)
