######################################################################
# This script evaluates the performance of four Frechet regression 
# methods (naive-nw, naive-knn, semi-nw, semi-knn) on the noisy 
# face dataset. The regression task maps from the 4096-dimensional 
# Euclidean space (64x64 images) to the spherical data space.  
# The script also visualizes how a facial image change under noise contamination.
#
# Adjustable parameters:
#   - n1 (in Get_sphere function): number of labeled samples, 70 or 140  
#   - ratios : noise level or contamination proportion, [0,0.3] or [0,1]  
#
# Outputs:
#   - predicted face orientations (test error measured by spherical geodesic distance)  
######################################################################


library(parallel)
library(pbapply)
library("R.matlab") 
source("function.R")

#----------data Preparation----------
data <- readMat("face_data.mat")
Xdata_clean <- t(data$images) #pixels are taken by image columns
Ydata_sph <- t(data$poses) #the lower right is a positive angle, and the upper left is a negative angle
Ydata <- matrix(0, 698, 3) 
#convert spherical coordinates to rectangular coordinates############
for (i in 1:698) {
  Ydata[i,1] <- cos(Ydata_sph[i,2]*pi/180)*cos(Ydata_sph[i,1]*pi/180)
  Ydata[i,2] <- cos(Ydata_sph[i,2]*pi/180)*sin(Ydata_sph[i,1]*pi/180)
  Ydata[i,3] <- -sin(Ydata_sph[i,2]*pi/180)
}  

#----------Get_sphere function----------
Get_sphere<-function(loop){
  library("igraph")
  library("Matrix")
  library("matrixStats")
  library("Directional")
  library("Riemann")
  library("spherepc")
  source("function.R")
  
  set.seed(loop)
  oldw <- getOption("warn")
  options(warn = -1)
  n1=70; n2=698-n1
  MISE1 = MISE2 = MISE3 <- rep(0,4)
  
  ##############data generation###################################
  label_index <- sample(698,n1)
  unlabel_index <- (1:698)[-label_index]
  train_label <- list(X=Xdata[label_index,],Y=Ydata[label_index,])
  test <- list(X=Xdata[unlabel_index,],Y=Ydata[unlabel_index,])
  test_sph <- Ydata_sph[unlabel_index,]
  Xdist_temp <- Xdist[c(label_index,unlabel_index),c(label_index,unlabel_index)]
  graph_dis_temp <- graph_dis[c(label_index,unlabel_index),c(label_index,unlabel_index)]
  
  ##############supervised learning###################################
  dist <- Xdist_temp[1:n1,1:n1]
  #leave one cross-validation for knn#######
  krange <- 1:10
  k <- knn_validation_sphere(krange=krange, sphere=train_label$Y, dis=dist)$k_opt
  #leave one cross-validation for nw########
  diag(dist) <- Inf
  h_0 <- median(apply(dist,2,min))
  diag(dist) <- 0
  hrange <- h_0*5^(seq(0,1,by=0.1))
  h <- nw_validation_sphere(hrange=hrange, sphere=train_label$Y, dis=dist)$h_opt
  #supervised learning with seleted k,h######################
  knn=nw <- matrix(0,n2,3)
  for (i in 1:n2) {
    E_dist_to_label <- Xdist_temp[1:n1,(n1+i)]
    #knn############################
    knn_index <- order(E_dist_to_label)[1:k]
    if (length(knn_index)==1) {knn[i,] <- train_label$Y[knn_index,]}
    if (length(knn_index)>1) {
      myriem = wrap.sphere(train_label$Y[knn_index,])
      knn[i,] <- c(riem.mean(myriem, weight=rep(1,length(knn_index)), geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
    #nw#############################
    loc_weight <- sapply(E_dist_to_label, function(z) K_Epa(z,h=h))
    if (length(which(loc_weight>0))==0) {nw[i,] <- train_label$Y[n1,]}
    if (length(which(loc_weight>0))>0) {
      myriem = wrap.sphere(train_label$Y)
      nw[i,] <- c(riem.mean(myriem, weight=loc_weight, geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
  }
  
  # ##############semi-supervised learning###############################
  dist <- graph_dis_temp[1:n1,1:n1]
  #leave one cross-validation for semi-knn#######
  k_semi <- knn_validation_sphere(krange=krange, sphere=train_label$Y, dis=dist)$k_opt
  #leave one cross-validation for semi-nw########
  diag(dist) <- Inf
  h_0 <- median(apply(dist,2,min))
  diag(dist) <- 0
  hrange <- h_0*5^(seq(0,1,by=0.1))
  h_semi<- nw_validation_sphere(hrange=hrange, sphere=train_label$Y, dis=dist)$h_opt
  #semi-supervised learning with seleted k_semi, h_semi################################################
  knn_semi=nw_semi <- matrix(0,n2,3)
  for (i in 1:n2) {
    G_dist_to_label <- graph_dis_temp[1:n1,(n1+i)]
    #knn-semi###################
    knn_semi_index <- intersect(order(G_dist_to_label)[1:k_semi],which(G_dist_to_label<10000))
    if (length(knn_semi_index)==0) {knn_semi[i,] <- train_label$Y[n1,]}
    if (length(knn_semi_index)==1) {knn_semi[i,] <- train_label$Y[knn_semi_index,]}
    if (length(knn_semi_index)>1) {
      myriem = wrap.sphere(train_label$Y[knn_semi_index,])
      knn_semi[i,] <- c(riem.mean(myriem, weight=rep(1,length(knn_semi_index)), geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
    #nw-semi###################
    loc_weight_semi <- sapply(G_dist_to_label, function(z) K_Epa(z,h=h_semi))
    if (length(which(loc_weight_semi>0))==0) {nw_semi[i,] <- train_label$Y[n1,]}
    if (length(which(loc_weight_semi>0))>0) {
      myriem = wrap.sphere(train_label$Y)
      nw_semi[i,] <- c(riem.mean(myriem, weight=loc_weight_semi, geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
  }
  
  #########orientation error######################################
  error <- matrix(0,n2,4)
  for (i in 1:n2){
    error[i,1] <- acos(crossprod(knn[i,],test$Y[i,])[1,1])^2
    error[i,2] <- acos(crossprod(knn_semi[i,],test$Y[i,])[1,1])^2
    error[i,3] <- acos(crossprod(nw[i,],test$Y[i,])[1,1])^2
    error[i,4] <- acos(crossprod(nw_semi[i,],test$Y[i,])[1,1])^2
  } 
  MISE1 <- colMeans(error)
  ###########left-right error/up-down error################################
  knn_sph <- matrix(c(atan(knn[,2]/knn[,1])*180/pi,-asin(knn[,3])*180/pi),n2,2)
  knn_semi_sph <- matrix(c(atan(knn_semi[,2]/knn_semi[,1])*180/pi,-asin(knn_semi[,3])*180/pi),n2,2)
  nw_sph <- matrix(c(atan(nw[,2]/nw[,1])*180/pi,-asin(nw[,3])*180/pi),n2,2)
  nw_semi_sph <- matrix(c(atan(nw_semi[,2]/nw_semi[,1])*180/pi,-asin(nw_semi[,3])*180/pi),n2,2)
  MISE2[1] <- colMeans(abs(knn_sph-test_sph))[1]
  MISE2[2] <- colMeans(abs(knn_semi_sph-test_sph))[1]
  MISE2[3] <- colMeans(abs(nw_sph-test_sph))[1]
  MISE2[4] <- colMeans(abs(nw_semi_sph-test_sph))[1]
  MISE3[1] <- colMeans(abs(knn_sph-test_sph))[2]
  MISE3[2] <- colMeans(abs(knn_semi_sph-test_sph))[2]
  MISE3[3] <- colMeans(abs(nw_sph-test_sph))[2]
  MISE3[4] <- colMeans(abs(nw_semi_sph-test_sph))[2]
  
  return(list(MISE1=MISE1, MISE2=MISE2, MISE3=MISE3))
}

#----------Add noise to images----------
add_noise_flat_region <- function(face_matrix, noise_sd = 0.1, pollute_ratio = 0.2, img_size = 64) {
  # face_matrix: n_images × n_pixels (flatten)
  n_images <- nrow(face_matrix)
  n_pixels <- ncol(face_matrix)
  # Calculate the number of rows to be polluted
  n_rows_pollute <- ceiling(img_size * pollute_ratio)
  # Generate a mask: 1 indicates adding noise, 0 indicates no noise
  mask <- matrix(0, nrow = img_size, ncol = img_size)
  if (n_rows_pollute >= 1){mask[1:n_rows_pollute, ] <- 1}
  mask_flat <- as.vector(mask)  # Flatten
  # Generate a noise matrix
  noise <- matrix(rnorm(length(face_matrix), mean = 0, sd = noise_sd),
                  nrow = n_images, ncol = n_pixels)
  # Add noise only at positions where mask = 1
  noisy_matrix <- face_matrix + noise * matrix(mask_flat, nrow = n_images, ncol = n_pixels, byrow = TRUE)
  # Clip values to [0, 1]
  noisy_matrix[noisy_matrix < 0] <- 0
  noisy_matrix[noisy_matrix > 1] <- 1
  
  return(noisy_matrix)
}

#----------------------main code-----------------------------
set.seed(1)
ratios <- seq(0, 1, 0.2)  #noise level
#ratios <- seq(0, 0.3, 0.06)  #pollution ratio
result <- matrix(0,length(ratios),4)
for (kk in 1:length(ratios)) {
  Xdata <- add_noise_flat_region(Xdata_clean, noise_sd = 0.3, pollute_ratio = ratios[kk], img_size = 64) 
  Xdist <- generate_distance(A=Xdata)
  #############graph distance##########################################################
  dist <- Xdist
  poi <- matrix(0,698,698)
  for (i in 1:698) {
    poi[i,order(dist[i,])[1:5]] <- 1
  }
  for (i in 1:698) {
    for (j in 1:698) {
      if (poi[i,j]==0 && poi[j,i]==0){dist[i,j]=dist[j,i]=0}
    }
  }
  dist_sparse <- as(dist, "dgCMatrix")
  g <- graph_from_adjacency_matrix(dist_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
  graph_dis <- distances(g, v = 1:698, to = 1:698)
  
  #-----------parallel run-----------
  simu <- 100
  system.time({
    cores <- min(25, parallel::detectCores() - 1)
    cl <- makeCluster(cores)
    clusterExport(cl, c("Xdata","Ydata_sph","Ydata","Xdist","graph_dis"))
    res <- pblapply(1:simu, Get_sphere, cl = cl)
    stopCluster(cl)
  }) 
  
  #----------final results----------
  MISE <- matrix(0,simu,4)
  for (i in 1:simu) {
    MISE[i,] <- res[[i]]$MISE1
  }

  result[kk,] <- colMeans(MISE)
}

write.csv(x = result,file = "real data_noise_ratio_140.csv")

# ===== plot noisy face images =====
par(mfrow = c(2, length(ratios)), mai=c(0.01,0.01,0.01,0.01))
ratios <- seq(0, 1, 0.2) #noise level
#ratios <- seq(0, 0.3, 0.06) #pollution ratio
for (i in ratios) {
  Xdata <- add_noise_flat_region(Xdata_clean, noise_sd = 0.3, pollute_ratio = i, img_size = 64)
  image(matrix(Xdata[1, ], nrow = 64, byrow = TRUE)[,64:1], 
        col = gray(seq(0, 1, length = 255)), axes = FALSE)
}
