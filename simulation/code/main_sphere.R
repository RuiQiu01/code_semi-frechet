######################################################################
# This code compares six Frechet regression methods (naive-nw, naive-knn, semi-knn, semi-nw)
# for Settings III and V in the paper, where the responses are spherical data. The implementation
# is based on a main function with parallel computation.
#
# Adjustable parameters:
#   - scenario: 
#       "Swiss"  : X on a 2D manifold in R^3 (uniform latent U); 
#       "Swiss2" : X on a 2D manifold in R^6 (uniform U); 
#       "Swiss3" : X on a 2D manifold in R^3 (truncated normal U); 
#       "Swiss4" : X on a 2D manifold in R^6 (truncated normal U); 
#       "Swiss5" : corresponds to Setting V in the supplement, where X lies 
#                  in R^6 without low-dimensional manifold structure.
#
#   - kappa : von Mises-Fisher distribution (vMF), 0 = vMF(·,0); 5 = vMF(·,5).  
#   - sigma : 0.25, noise level used in the paper.  
#   - n1, n2, n3: number of labeled data, unlabeled data, and test samples. 
######################################################################


library(parallel)
library(pbapply)

#----------main function---------------
Get_sphere<-function(loop){
  library("igraph")
  library("Matrix")
  library("matrixStats")
  library("Directional")
  library("tmvtnorm")
  library("Riemann")
  library("spherepc")
  source("function.R")
  
  set.seed(loop)
  oldw <- getOption("warn")
  options(warn = -1)
  scenario="Swiss"; kappa=0; sigma=0.25; n1=100; NN=8000; n3=1000 #n1, NN, n3 : number of labeled data, unlabeled data, and test samples.  
  semigroup <- c(0,100,200,500,1000,2000,3000,4000,5000,6000,7000,NN)
  ISE <- matrix(0,length(semigroup),4)
  
  #----------data generation----------
  train_label <- generate_label_sphere(scenario=scenario, N=n1, kappa=kappa, sigma=sigma)
  train_unlabel <- generate_unlabel_sphere(scenario=scenario, N=NN)
  test <- generate_test_sphere(scenario=scenario, N=n3)
  X_train <- rbind(train_label$X, train_unlabel$X)
  dist_total <- generate_distance(A=X_train)
  dist_test_total <- generate_distance_two(A=test$X, B=X_train)

  
  #----------supervised regression----------
  dist=dist_total[1:n1,1:n1]
  ##leave one cross-validation for knn
  krange <- seq(1,10,1)
  k <- knn_validation_sphere(krange=krange, sphere=train_label$Y, dis=dist)$k_opt
  ##leave one cross-validation for nw
  diag(dist) <- Inf
  h_0 <- median(rowMins(dist))
  diag(dist) <- 0
  hrange <- h_0*5^(seq(0,1,by=0.1))
  h <- nw_validation_sphere(hrange=hrange, sphere=train_label$Y, dis=dist)$h_opt
  ##supervised regression with seleted k,h
  knn=nw <- matrix(0,n3,3)
  for (i in 1:n3) {
    E_dist_to_label <- dist_test_total[i,1:n1]
    ###knn
    knn_index <- order(E_dist_to_label)[1:k]
    if (length(knn_index)==1) {knn[i,] <- train_label$Y[knn_index,]}
    if (length(knn_index)>1) {
      myriem = wrap.sphere(train_label$Y[knn_index,])
      knn[i,] <- c(riem.mean(myriem, weight=rep(1,length(knn_index)), geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
    ###nw
    loc_weight <- sapply(E_dist_to_label, function(z) K_Epa(z,h=h))
    if (length(which(loc_weight>0))==0) {nw[i,] <- train_label$Y[n1,]}
    if (length(which(loc_weight>0))>0) {
      myriem = wrap.sphere(train_label$Y)
      nw[i,] <- c(riem.mean(myriem, weight=loc_weight, geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
  }
  ##testing error for supervised regerssion
  error_super <- matrix(0,n3,2)
  for (i in 1:n3){
    error_super[i,1] <- acos(crossprod(knn[i,],test$Y[i,])[1,1])^2
    error_super[i,2] <- acos(crossprod(nw[i,],test$Y[i,])[1,1])^2
  } 
  ISE_super <- colMeans(error_super) 
  ISE[, 1:2] <- rep(ISE_super, each = nrow(ISE))
  
  #----------semi-supervised regression----------
  b <- 1
  for (n2 in semigroup) {
    dist <- dist_total[1:(n1+n2),1:(n1+n2)]
    dist_test <- dist_test_total[,1:(n1+n2)]
    ##cross-validation for graph radius
    diag(dist) <- Inf
    R <- max(rowMins(dist))
    diag(dist) <- 0
    #dist[dist > 1.2*R] <- Inf
    #dist_test[dist_test > 1.2*R] <- Inf
    rrange <- R*(seq(0.5,2,by=0.1))
    list_knn <- list(); list_nw <- list(); count = 1
    for (r in rrange) {
      dist_temp <- dist
      dist_temp[dist_temp > r] <- Inf
      ##matrix sparse transformation
      dist_temp[is.infinite(dist_temp)] <- 0
      dist_sparse <- as(dist_temp, "dgCMatrix")
      g <- graph_from_adjacency_matrix(dist_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
      graph_dis <- distances(g, v = 1:n1, to = 1:n1)
      ##leave one cross-validation for semi-knn 
      krange <- seq(1,10,1)
      list_knn[[count]] <- knn_validation_sphere(krange=krange, sphere=train_label$Y, dis=graph_dis)
      ##leave one cross-validation for semi-nw 
      diag(graph_dis) <- Inf
      h_0 <- median(rowMins(graph_dis))
      diag(graph_dis) <- 0
      hrange <- h_0*5^(seq(0,1,by=0.1))
      list_nw[[count]] <- nw_validation_sphere(hrange=hrange, sphere=train_label$Y, dis=graph_dis)
      count <- count+1
    }
    rm(dist_temp);gc()
    knn_index <- which.min(sapply(list_knn, function(x) x$err_opt))
    r_knn <- rrange[knn_index]
    k_semi <- list_knn[[knn_index]]$k_opt
    nw_index <- which.min(sapply(list_nw, function(x) x$err_opt))
    r_nw <- rrange[nw_index]
    h_semi <- list_nw[[nw_index]]$h_opt
    ##semi-supervised knn regression with seleted r, k_semi
    dist1 <- dist; dist_test1 <- dist_test
    dist1[dist1 > r_knn] <- Inf
    dist_test1[dist_test1 > r_knn] <- Inf
    ###matrix sparse transformation
    dist1[is.infinite(dist1)] <- 0
    dist_sparse <- as(dist1, "dgCMatrix")
    g <- graph_from_adjacency_matrix(dist_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
    graph_dis_total <- distances(g, v = 1:n1, to = 1:(n1+n2))
    knn_semi <- matrix(0,n3,3)
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test1[i,],n1), nrow = n1, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-knn
      knn_semi_index <- intersect(order(G_dist_to_label)[1:k_semi],which(G_dist_to_label<Inf))
      if (length(knn_semi_index)==0) {knn_semi[i,] <- train_label$Y[n1,]}
      if (length(knn_semi_index)==1) {knn_semi[i,] <- train_label$Y[knn_semi_index,]}
      if (length(knn_semi_index)>1) {
        myriem = wrap.sphere(train_label$Y[knn_semi_index,])
        knn_semi[i,] <- c(riem.mean(myriem, weight=rep(1,length(knn_semi_index)), geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
    }
    rm(list = c("dist1","dist_test1"));gc()
    ##semi-supervised nw regression with seleted r, h_semi
    dist2 <- dist; dist_test2 <- dist_test
    dist2[dist2 > r_nw] <- Inf
    dist_test2[dist_test2 > r_nw] <- Inf
    ###matrix sparse transformation
    dist2[is.infinite(dist2)] <- 0
    dist_sparse <- as(dist2, "dgCMatrix")
    g <- graph_from_adjacency_matrix(dist_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
    graph_dis_total <- distances(g, v = 1:n1, to = 1:(n1+n2))
    nw_semi <- matrix(0,n3,3)
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test2[i,],n1), nrow = n1, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-nw
      loc_weight_semi <- sapply(G_dist_to_label, function(z) K_Epa(z,h=h_semi))
      if (length(which(loc_weight_semi>0))==0) {nw_semi[i,] <- train_label$Y[n1,]}
      if (length(which(loc_weight_semi>0))>0) {
        myriem = wrap.sphere(train_label$Y)
        nw_semi[i,] <- c(riem.mean(myriem, weight=loc_weight_semi, geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
    }
    rm(list = c("dist2","dist_test2"));gc()
    ##testing error for semi-supervised regerssion
    error_semi <- matrix(0,n3,2)
    for (i in 1:n3){
      error_semi[i,1] <- acos(crossprod(knn_semi[i,],test$Y[i,])[1,1])^2
      error_semi[i,2] <- acos(crossprod(nw_semi[i,],test$Y[i,])[1,1])^2
    } 
    ISE_semi <- colMeans(error_semi) 
    ISE[b,3:4] <- ISE_semi 
    b=b+1
  }
  return(ISE)
}

#----------parallel computation----------
system.time({
  cores <- min(25, parallel::detectCores() - 1)
  cl <- makeCluster(cores)
  res <- pblapply(1:100, Get_sphere, cl = cl)
  stopCluster(cl)
}) 

#----------final results----------
MISE <- array(0,c(12,4,100))
for (i in 1:100) {
  MISE[,,i] <- res[[i]]
}
compare <- matrix(0,12,9)
for (i in 1:12){
  compare[i,1:4] <- rowMeans(MISE[i,,]) #mean of MSE (AMSE)
  compare[i,6:9] <- sqrt(rowVars(MISE[i,,])) #variance of MSE (AMSE)
} 
#write.csv(x = compare,file = "sphere-100-8000-1000-Swiss-setting1_0.25.csv")







