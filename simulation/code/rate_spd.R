######################################################################
# This code implements the experiments in the supplementary material that evaluate the empirical convergence rates of semi-nw and semi-knn 
# methods by varying the number of labeled samples, where the responses are spd matrices (Setting I and I-2). The implementation is based on a main function with parallel computation.
#
#  - scenario:
#       "Swiss"  : X on a 2D manifold in R^3 (uniform latent U); 
#       "Swiss2" : X on a 2D manifold in R^6 (uniform U); 
#  - setting: 1 = 2D SPD responses (Setting I); 4 = 2D SPD responses (S.5 in Appendix S.3.4);
#  - snr=2, corresponding to two signal-to-noise ratios.
#  - n1, n2, n3: number of unlabeled data, labeled data, and test samples.   
######################################################################


library(parallel)
library(pbapply)

#----------main function----------
Get_spd<-function(loop){
  library(igraph)
  library(shapes)
  library(Matrix)
  library(matrixStats)
  library(RiemBase)
  library(tmvtnorm)
  source("function.R")
  
  set.seed(loop)
  scenario="Swiss2"; setting=1; snr=2; n1=3000; NN=1500; n3=1000
  method <<- "LogEuclidean"; alpha <<- 1
  semigroup <- seq(500,1200,100)
  ISE <- matrix(0,length(semigroup),2)
  
  #----------data generation----------
  train_label <- generate_label_spd(scenario=scenario, setting=setting, N=NN, snr=snr)
  train_unlabel <- generate_unlabel_spd(scenario=scenario, N=n1)
  test <- generate_test_spd(scenario=scenario, setting=setting, N=n3)
  X_train <- rbind(train_unlabel$X,train_label$X)
  dist_total <- generate_distance(A=X_train)
  dist_test_total <- generate_distance_two(A=test$X, B=X_train)
  m <- dim(train_label$Y)[1]
  
  #----------semi-supervised regression----------
  b <- 1
  for (n2 in semigroup) {
    dist <- dist_total[1:(n1+n2),1:(n1+n2)]
    dist_test <- dist_test_total[,1:(n1+n2)]
    ##cross-validation for graph radius
    diag(dist) <- Inf
    R <- max(rowMins(dist))
    diag(dist) <- 0
    rrange <- R*(seq(0.5,2,by=0.1))
    list_knn <- list(); list_nw <- list(); count = 1
    for (r in rrange) {
      dist_temp <- dist
      dist_temp[dist_temp > r] <- Inf
      ##matrix sparse transformation
      dist_temp[is.infinite(dist_temp)] <- 0
      dist_sparse <- as(dist_temp, "dgCMatrix")
      g <- graph_from_adjacency_matrix(dist_sparse, mode = "undirected", weighted = TRUE, diag = FALSE)
      graph_dis <- distances(g, v = (n1+1):(n1+n2), to = (n1+1):(n1+n2))
      ##leave one cross-validation for semi-knn
      krange <- seq(round(sqrt(n2))-8,round(sqrt(n2))+8,1)
      list_knn[[count]] <- knn_validation_spd(krange=krange, spd=train_label$Y[,,1:n2], dis=graph_dis)
      ##leave one cross-validation for semi-nw
      diag(graph_dis) <- Inf
      h_0 <- median(rowMins(graph_dis))
      diag(graph_dis) <- 0
      hrange <- h_0*seq(1,10,length.out=16)
      list_nw[[count]] <- nw_validation_spd(hrange=hrange, spd=train_label$Y[,,1:n2], dis=graph_dis)
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
    graph_dis_total <- distances(g, v = (n1+1):(n1+n2), to = 1:(n1+n2))
    knn_semi <- array(0,c(m,m,n3))
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test1[i,],n2), nrow = n2, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-knn
      knn_semi_index <- intersect(order(G_dist_to_label)[1:k_semi],which(G_dist_to_label<Inf))
      if (length(knn_semi_index)==0) {knn_semi[,,i] <- train_label$Y[,,n2]}
      if (length(knn_semi_index)==1) {knn_semi[,,i] <- train_label$Y[,,1:n2][,,knn_semi_index]}
      if (length(knn_semi_index)>1) {knn_semi[,,i] <- spd_predict(Y=train_label$Y[,,1:n2][,,knn_semi_index], method=method, weight=1, alpha=alpha)}
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
    graph_dis_total <- distances(g, v = (n1+1):(n1+n2), to = 1:(n1+n2))
    nw_semi <- array(0,c(m,m,n3))
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test2[i,],n2), nrow = n2, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-nw
      loc_weight_semi <- sapply(G_dist_to_label, function(z) K_Epa(z,h=h_semi))
      if (length(which(loc_weight_semi>0))==0) {nw_semi[,,i] <- train_label$Y[,,n2]}
      if (length(which(loc_weight_semi>0))>0) {nw_semi[,,i] <- spd_predict(Y=train_label$Y[,,1:n2], method=method, weight=loc_weight_semi, alpha=alpha)}
    }
    rm(list = c("dist2","dist_test2"));gc()
    ##testing error for semi-supervised regerssion
    error_semi <- matrix(0,n3,2)
    for (i in 1:n3){
      error_semi[i,1] <- distcov(knn_semi[,,i], test$Y[,,i], method=method, alpha=alpha)^2
      error_semi[i,2] <- distcov(nw_semi[,,i], test$Y[,,i], method=method, alpha=alpha)^2
    } 
    ISE[b,] <- colMeans(error_semi)
    b=b+1
  }
  return(ISE)
}

#----------parallel computation----------
system.time({
  cores <- min(25, parallel::detectCores() - 1)
  cl <- makeCluster(cores)
  res <- pblapply(1:100, Get_spd, cl = cl)
  stopCluster(cl)
}) 

#----------final results----------
MISE <- array(0,c(8,2,100))
for (i in 1:100) {
  MISE[,,i] <- res[[i]]
}
compare <- matrix(0,8,2)
for (i in 1:8){
  compare[i,] <- rowMeans(MISE[i,,]) #AMSE
} 

#write.csv(x = compare,file = "spd-1200-3000-1000-Swiss2-setting1-snr_2_rate.csv")






