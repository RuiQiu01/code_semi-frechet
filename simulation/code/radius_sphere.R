######################################################################
# This code implements the experiments in the supplementary material that study the effect of graph radius on the performance of semi-nw  
# and semi-knn methods, where the responses are spherical data. The implementation is based on a main function with parallel computation.
#  
#  - scenario:
#       "Swiss"  : X on a 2D manifold in R^3 (uniform latent U); 
#       "Swiss2" : X on a 2D manifold in R^6 (uniform U); 
#  - kappa : von Mises-Fisher distribution (vMF), 0 = vMF(·,0); 5 = vMF(·,5).  
#  - sigma : 0.25, noise level used in the paper.  
#  - n1, n2, n3: number of labeled data, unlabeled data, and test samples. 
######################################################################



library(parallel)
library(pbapply)

#----------main function----------
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
  scenario="Swiss2"; kappa=0; sigma=0.25; n1=200; n2=3000; n3=1000
  ISE_semi <- matrix(0,16,2)
  
  #----------data generation----------
  train_label <- generate_label_sphere(scenario=scenario, N=n1, kappa=kappa, sigma=sigma)
  train_unlabel <- generate_unlabel_sphere(scenario=scenario, N=n2)
  test <- generate_test_sphere(scenario=scenario, N=n3)
  X_train <- rbind(train_label$X, train_unlabel$X)
  dist_total <- generate_distance(A=X_train)
  dist_test_total <- generate_distance_two(A=test$X, B=X_train)
  
  ##Impaction of graph radius
  diag(dist_total) <- Inf
  R <- max(rowMins(dist_total))
  diag(dist_total) <- 0
  rrange <- R*(seq(0.5,2,by=0.1))
  #eta <- 1.661; tau <- 0.628
  #rrange <- seq(0.06, min(tau,1/(3*eta)), by = 0.01)
  b <- 1
  for (r in rrange) {
    dist <- dist_total; dist_test <- dist_test_total
    dist[dist > r] <- Inf
    dist_test[dist_test > r] <- Inf
    ##matrix sparse transformation##################
    dist[is.infinite(dist)] <- 0
    dist_sparse <- as(dist, "dgCMatrix")
    g <- graph_from_adjacency_matrix(dist_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
    graph_dis_total <- distances(g, v = 1:n1, to = 1:(n1+n2))
    graph_dis <- graph_dis_total[,1:n1]
    ##leave one cross-validation for semi-knn #######
    krange <- 1:10
    k_semi <- knn_validation_sphere(krange=krange, sphere=train_label$Y, dis=graph_dis)$k_opt
    ##leave one cross-validation for semi-nw ########
    diag(graph_dis) <- Inf
    h_0 <- median(rowMins(graph_dis))
    diag(graph_dis) <- 0
    hrange <- h_0*5^(seq(0,1,by=0.1))
    h_semi <- nw_validation_sphere(hrange=hrange, sphere=train_label$Y, dis=graph_dis)$h_opt
    ##semi-supervised regression with seleted k_semi, h_semi################################################
    knn_semi = nw_semi <- matrix(0,n3,3)
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test[i,],n1), nrow = n1, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-knn###################
      knn_semi_index <- intersect(order(G_dist_to_label)[1:k_semi],which(G_dist_to_label<Inf))
      if (length(knn_semi_index)==0) {knn_semi[i,] <- train_label$Y[n1,]}
      if (length(knn_semi_index)==1) {knn_semi[i,] <- train_label$Y[knn_semi_index,]}
      if (length(knn_semi_index)>1) {
        myriem = wrap.sphere(train_label$Y[knn_semi_index,])
        knn_semi[i,] <- c(riem.mean(myriem, weight=rep(1,length(knn_semi_index)), geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
      ###semi-nw
      loc_weight_semi <- sapply(G_dist_to_label, function(z) K_Epa(z,h=h_semi))
      if (length(which(loc_weight_semi>0))==0) {nw_semi[i,] <- train_label$Y[n1,]}
      if (length(which(loc_weight_semi>0))>0) {
        myriem = wrap.sphere(train_label$Y)
        nw_semi[i,] <- c(riem.mean(myriem, weight=loc_weight_semi, geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
    }
    ##testing error for semi-supervised regerssion
    error_semi <- matrix(0,n3,2)
    for (i in 1:n3){
      error_semi[i,1] <- acos(crossprod(knn_semi[i,],test$Y[i,])[1,1])^2
      error_semi[i,2] <- acos(crossprod(nw_semi[i,],test$Y[i,])[1,1])^2
    } 
    ISE_semi[b,] <- colMeans(error_semi) 
    b <- b+1
  }
  return(ISE_semi)
}

#----------parallel computation----------
system.time({
  cores <- min(25, parallel::detectCores() - 1)
  cl <- makeCluster(cores)
  res <- pblapply(1:100, Get_sphere, cl = cl)
  stopCluster(cl)
}) 

#----------final results----------
MISE <- array(0,c(16,2,100))
for (i in 1:100) {
  MISE[,,i] <- res[[i]]
}
compare <- matrix(0,16,2)
for (i in 1:16){
  compare[i,] <- rowMeans(MISE[i,,]) #AMSE
} 

#write.csv(x = compare,file = "sphere-200-3000-1000-Swiss2-setting1_0.25_r_impaction.csv")












