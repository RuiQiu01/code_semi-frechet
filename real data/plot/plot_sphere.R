######################################################################
# This script visualizes the predicted points on the sphere for the 
# real data analysis. 
#
# Methods: naive-knn，semi-knn，naive-nw，semi-nw
######################################################################

library(plot3D)

pp <- c(541,419,210,537,533,339,428,217)
pc <- rep(0,8)
for (i in 1:8) {
  pc[i] <- which(unlabel_index==pp[i])
}

##spherical coordinates############################################################################
M  <- mesh(seq(0, 2*pi, length.out = 100), 
           seq(0,   pi, length.out = 100))
u  <- M$x ; v  <- M$y

xp <- cos(u)*sin(v)
yp <- sin(u)*sin(v)
zp <- cos(v)

par(mfrow = c(1,4), mai=c(0.21,0.21,0.21,0.21))

##knn##########################################################
scatter3D(xp, yp, zp, 
          phi = 1,
          pch = 20, 
          cex = 1, 
          col = 'gray',
          alpha = 0.1,
          ticktype = "detailed",bty = "b", colkey = F, theta = 90)

Y_true <- test$Y[pc,]
xs <- Y_true[,1]
ys <- Y_true[,2]
zs <- Y_true[,3]
scatter3D(xs, ys, zs, pch = 20, cex=1, col = "red", add = TRUE)

Y1 <- knn[pc,]
x1 <- Y1[,1]
y1 <- Y1[,2]
z1 <- Y1[,3]
scatter3D(x1, y1, z1, pch = 20, cex=1, col = "black", add = TRUE)

##semi-knn##########################################################
scatter3D(xp, yp, zp, 
          phi = 1,
          pch = 20, 
          cex = 1, 
          col = 'gray',
          alpha = 0.1,
          ticktype = "detailed",bty = "b", colkey = F, theta = 90)

Y_true <- test$Y[pc,]
xs <- Y_true[,1]
ys <- Y_true[,2]
zs <- Y_true[,3]
scatter3D(xs, ys, zs, pch = 20, cex=1, col = "red", add = TRUE)

Y2 <- knn_semi[pc,]
x2 <- Y2[,1]
y2 <- Y2[,2]
z2 <- Y2[,3]
scatter3D(x2, y2, z2, pch = 20, cex=1, col = "black", add = TRUE)

##nw##########################################################
scatter3D(xp, yp, zp, 
          phi = 1,
          pch = 20, 
          cex = 1, 
          col = 'gray',
          alpha = 0.1,
          ticktype = "detailed",bty = "b", colkey = F, theta = 90)

Y_true <- test$Y[pc,]
xs <- Y_true[,1]
ys <- Y_true[,2]
zs <- Y_true[,3]
scatter3D(xs, ys, zs, pch = 20, cex=1, col = "red", add = TRUE)

Y3 <- nw[pc,]
x3 <- Y3[,1]
y3 <- Y3[,2]
z3 <- Y3[,3]
scatter3D(x3, y3, z3, pch = 20, cex=1, col = "black", add = TRUE)

##semi-nw##########################################################
scatter3D(xp, yp, zp, 
          phi = 1,
          pch = 20, 
          cex = 1, 
          col = 'gray',
          alpha = 0.1,
          ticktype = "detailed",bty = "b", colkey = F, theta = 90)

Y_true <- test$Y[pc,]
xs <- Y_true[,1]
ys <- Y_true[,2]
zs <- Y_true[,3]
scatter3D(xs, ys, zs, pch = 20, cex=1, col = "red", add = TRUE)
Y4 <- nw_semi[pc,]
x4 <- Y4[,1]
y4 <- Y4[,2]
z4 <- Y4[,3]
scatter3D(x4, y4, z4, pch = 20, cex=1, col = "black", add = TRUE)

