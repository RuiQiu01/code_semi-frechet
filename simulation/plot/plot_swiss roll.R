######################################################################
# This code generates and visualizes the Swiss roll dataset. It is 
# used for illustrating manifold structures in simulation studies. 
######################################################################

library(scatterplot3d)

# Generate Swiss roll data
set.seed(123)
n <- 2000
t <- (4 * pi) * (1/2 + runif(n))
x <- t * cos(t)/10
y <- 4 * runif(n)
z <- t * sin(t)/10

# Background point grayscale levels
gray_levels <- 0.1 + 0.8 * (rank(z) - 1) / (n - 1)
colors <- gray(gray_levels)

# Special points
t_points <- (c(4/30, 19/30, 12/30)+1/2) * 4 * pi
x_points <- t_points * cos(t_points)/10
y_points <- rep(1, 3)   # y=0
z_points <- t_points * sin(t_points)/10
labels <- c("A", "B", "C")

# Find the grayscale corresponding to the special points
nearest_idx <- sapply(t_points, function(tp) which.min(abs(t - tp)))
gray_for_points <- gray_levels[nearest_idx]-0.2

# Generate “red with grayscale tint” colors
point_colors <- rgb(red = 1, green = gray_for_points, blue = gray_for_points)

png("swiss_roll_points.png", width = 2400, height = 1800, res = 300)

# Plot background points
s3d <- scatterplot3d(x, y, z,
                     color = colors,
                     pch = 16,
                     xlab = expression(X^{(1)}),
                     ylab = expression(X^{(2)}),
                     zlab = expression(X^{(3)}),
                     angle = 45,
                     highlight.3d = FALSE,
                     box = FALSE)

# Add “red with grayscale tint” points
s3d$points3d(x_points, y_points, z_points, col = point_colors, pch = 16)

# Add red labels
coords <- s3d$xyz.convert(x_points, y_points, z_points)
text(coords$x, coords$y, labels = labels, pos = 4, col = "red", cex = 1.4)

dev.off()
