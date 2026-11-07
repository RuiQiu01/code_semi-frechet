######################################################################
# This code plots AMSE against sample size on a log–log scale, in order to provide a 
# clear visualization of the convergence rate of semi-nw, and semi-knn 
# 
# - export figures for the supplementary material.  
######################################################################



# Load data
data <- as.matrix(read.csv('spd-1200-3000-1000-Swiss2-setting4-snr_2_rate.csv'))

# Prepare sample points and take log base 2
x <- log(seq(500,1200,100))
y <- log(data[,3])

# Set output to PDF file
pdf("spd-1200-3000-1000-Swiss2-setting4-snr_2_rate_nw.pdf", width = 9, height = 6) 

par(mar=c(5,6,4,2)+0.1)
# Plot the line chart (dashed line for data points)
#plotcol = "#66a61e" #semi-knn
plotcol = "#e6ab02" #semi-nw
plot(x, y, type="o", lty=2, lwd=2.7, col=plotcol, pch=16,
     xlab="ln(n)", ylab="ln(AMSE)",
     cex.lab=1.8, cex.axis=1.2,
     col.main="black", col.lab="black", col.axis="black",
     bg="white") 
grid()

# Fit a linear model (least squares)
model <- lm(y ~ x)

# Add fitted line or theoretical line
abline(model, col=plotcol, lwd=2.7, lty=1)
#custom_slope <- -0.387
#intercept <- mean(y) - custom_slope * mean(x)
#abline(a = intercept, b = custom_slope, col=plotcol, lwd=2.7, lty=1)

# Add grid lines
grid(col="grey80", lty=c(8,4), lwd=0.8)
# Draw black axes
box(col="black", lwd=1)

# Calculate R squared and slope
r_squared <- summary(model)$r.squared
slope <- coef(model)[2]

# Print the results
cat("Fitted line slope =", slope, "\n")
cat("R squared =", r_squared, "\n")

# Add legend 
# legend("topright",
#        legend = c(
#          expression(plain("ln(")*AMSE*plain(") vs ln(")*n*plain(")")),
#          expression(plain("ln(")*AMSE*plain(") = a · ln(" * n * plain(")")* " + b ; a") %~~% -0.661 ~ "and" ~ R^2 > 0.99)
#        ),
#        col=plotcol,
#        lty=c(2,1),
#        lwd=2.7,
#        pch=c(16, NA),
#        bty="o",         
#        box.col="grey80",
#        box.lwd=0.8,
#        bg=NA,
#        cex=1.4,        
#        pt.cex=1)

# Close the PDF device
dev.off()

