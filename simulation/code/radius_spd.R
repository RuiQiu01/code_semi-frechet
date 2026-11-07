######################################################################
# This code implements the experiments in the supplementary material that study the effect of graph radius on the performance of semi-nw
# and semi-knn methods, where the responses are spd matrices. The implementation is based on a main function with parallel computation.
#
#  - scenario:
#       "Swiss"  : X on a 2D manifold in R^3 (uniform latent U); 
#       "Swiss2" : X on a 2D manifold in R^6 (uniform U); 
#  - setting: 1 = 2D SPD responses; 2 = 3D SPD responses;
#  - snr: 2 or 4, corresponding to two signal-to-noise ratios.
#  - n1, n2, n3: number of labeled data, unlabeled data, and test samples.  
######################################################################


library(parallel)
library(pbapply)

#----------main results----------
Get_spd<-function(loop){
  library("FNN")
  library("pracma")
  library("igraph")
  library("shapes")
  library("Matrix")
  library("matrixStats")
  library("RiemBase")
  source("function.R")
  
  set.seed(loop)
  scenario="Swiss2"; setting=1; snr=2; n1=200; n2=3000; n3=1000
  method <<- "LogEuclidean"; alpha <- 1
  ISE_semi <- matrix(0,16,2)
  
  #----------data generation----------
  train_label <- generate_label_spd(scenario=scenario, setting=setting, N=n1, snr=snr)
  train_unlabel <- generate_unlabel_spd(scenario=scenario, N=n2)
  test <- generate_test_spd(scenario=scenario, setting=setting, N=n3)
  X_train <- rbind(train_label$X, train_unlabel$X)
  dist_total <- generate_distance(A=X_train)
  dist_test_total <- generate_distance_two(A=test$X, B=X_train)
  m <- dim(train_label$Y)[1]
  
  ##Impaction of graph radius
  diag(dist_total) <- Inf
  R <- max(rowMins(dist_total))
  diag(dist_total) <- 0
  rrange <- R*(seq(0.5,2,by=0.1))
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
    k_semi <- knn_validation_spd(krange=krange, spd=train_label$Y, dis=graph_dis)$k_opt
    ##leave one cross-validation for semi-nw ########
    diag(graph_dis) <- Inf
    h_0 <- median(rowMins(graph_dis))
    diag(graph_dis) <- 0
    hrange <- h_0*5^(seq(0,1,by=0.1))
    h_semi <- nw_validation_spd(hrange=hrange, spd=train_label$Y, dis=graph_dis)$h_opt
    ##semi-supervised regression with seleted k_semi, h_semi################################################
    knn_semi = nw_semi <- array(0,c(m,m,n3))
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test[i,],n1), nrow = n1, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-knn###################
      knn_semi_index <- intersect(order(G_dist_to_label)[1:k_semi],which(G_dist_to_label<Inf))
      if (length(knn_semi_index)==0) {knn_semi[,,i] <- train_label$Y[,,n1]}
      if (length(knn_semi_index)==1) {knn_semi[,,i] <- train_label$Y[,,knn_semi_index]}
      if (length(knn_semi_index)>1) {knn_semi[,,i] <- spd_predict(Y=train_label$Y[,,knn_semi_index], method=method, weight=1, alpha=alpha)}
      ###semi-nw###################
      loc_weight_semi <- sapply(G_dist_to_label, function(z) K_Epa(z,h=h_semi))
      if (length(which(loc_weight_semi>0))==0) {nw_semi[,,i] <- train_label$Y[,,n1]}
      if (length(which(loc_weight_semi>0))>0) {nw_semi[,,i] <- spd_predict(Y=train_label$Y, method=method, weight=loc_weight_semi, alpha=alpha)}
    }
    ##testing error for semi-supervised regerssion######################################
    error_semi <- matrix(0,n3,2)
    for (i in 1:n3){
      error_semi[i,1] <- distcov(knn_semi[,,i], test$Y[,,i], method=method, alpha=alpha)^2
      error_semi[i,2] <- distcov(nw_semi[,,i], test$Y[,,i], method=method, alpha=alpha)^2
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
  res <- pblapply(1:100, Get_spd, cl = cl)
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

#write.csv(x = compare,file = "spd-200-3000-1000-Swiss2-setting1_snr_2_r_impaction.csv")






