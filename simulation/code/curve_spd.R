######################################################################
# This code compares six Fréchet regression methods (oracle-nw, oracle-knn, 
# naive-nw, naive-knn, semi-knn, semi-nw) for Settings VI in the 
# paper, where responses are SPD matrices and features are on a 1D manifold in R^3.
#
# Adjustable parameters:
#   - snr:     2 or 4, corresponding to two signal-to-noise ratios.  
#   - method:  "LogEuclidean" = log-Euclidean metric; "LogCholesky" = log-Cholesky metric.
#   - n1, n2, n3: number of labeled data, unlabeled data, and test samples.   
######################################################################


library(parallel)
library(pbapply)

#----------main function---------------
Get_spd<-function(loop){
  library(igraph)
  library(shapes)
  library(Matrix)
  library(matrixStats)
  library(RiemBase)
  library(tmvtnorm)
  source("function.R")
  
  set.seed(loop)
  snr=2; n1=100; NN=8000; n3=1000 #n1, NN, n3 : number of labeled data, unlabeled data, and test samples.  
  method <<- "LogEuclidean"; alpha <- 1
  semigroup <- c(6000,7000,NN)
  ISE <- matrix(0,length(semigroup),6)
  
  #----------data generation----------
  train_label <- generate_label_spd_curve(N=n1, snr=snr)
  train_unlabel <- generate_unlabel_spd_curve(N=NN)
  test <- generate_test_spd_curve(N=n3)
  X_train <- rbind(train_label$X, train_unlabel$X)
  dist_total <- generate_distance(A=X_train)
  dist_test_total <- generate_distance_two(A=test$X, B=X_train)
  dist_oral <- abs(outer(train_label$U, train_label$U, "-"))*sqrt(36*pi^2+0.25)
  dist_test_oral <- abs(outer(test$U, train_label$U, "-"))*sqrt(36*pi^2+0.25)
  m <- dim(train_label$Y)[1]
  
  
  #----------oracle regression----------
  dist=dist_oral
  ##leave one cross-validation for oracle-knn
  krange <- seq(1,10,1)
  k_oral <- knn_validation_spd(krange=krange, spd=train_label$Y, dis=dist)$k_opt
  ##leave one cross-validation for oracle-nw
  diag(dist) <- Inf
  h_0 <- median(rowMins(dist))
  diag(dist) <- 0
  hrange <- h_0*(seq(10,30,by=2))
  h_oral <- nw_validation_spd(hrange=hrange, spd=train_label$Y, dis=dist)$h_opt
  ##oracle regression with seleted k,h
  knn_oral=nw_oral <- array(0,c(m,m,n3))
  for (i in 1:n3) {
    oracle_dist_to_label <- dist_test_oral[i,]
    ###knn
    knn_index <- order(oracle_dist_to_label)[1:k_oral]
    if (length(knn_index)==1) {knn_oral[,,i] <- train_label$Y[,,knn_index]}
    if (length(knn_index)>1) {knn_oral[,,i] <- spd_predict(Y=train_label$Y[,,knn_index], method=method, weight=1, alpha=alpha)}
    ###nw
    loc_weight <- sapply(oracle_dist_to_label, function(z) K_Epa(z,h=h_oral))
    if (length(which(loc_weight>0))==0) {nw_oral[,,i] <- train_label$Y[,,n1]}
    if (length(which(loc_weight>0))>0) {nw_oral[,,i] <- spd_predict(Y=train_label$Y, method=method, weight=loc_weight, alpha=alpha)}
  }
  ##testing error for oracle regerssion
  error_oral <- matrix(0,n3,2)
  for (i in 1:n3){
    error_oral[i,1] <- distcov(knn_oral[,,i], test$Y[,,i], method=method, alpha=alpha)^2
    error_oral[i,2] <- distcov(nw_oral[,,i], test$Y[,,i], method=method, alpha=alpha)^2
  } 
  ISE_oral <- colMeans(error_oral) 
  ISE[, 1:2] <- rep(ISE_oral, each = nrow(ISE))
  
  #----------supervised regression----------
  dist=dist_total[1:n1,1:n1]
  ##leave one cross-validation for knn
  krange <- seq(1,10,1)
  k <- knn_validation_spd(krange=krange, spd=train_label$Y, dis=dist)$k_opt
  ##leave one cross-validation for nw
  diag(dist) <- Inf
  h_0 <- median(rowMins(dist))
  diag(dist) <- 0
  hrange <- h_0*5^(seq(0,1,by=0.1))
  h <- nw_validation_spd(hrange=hrange, spd=train_label$Y, dis=dist)$h_opt
  ##supervised regression with seleted k,h
  knn=nw <- array(0,c(m,m,n3))
  for (i in 1:n3) {
    E_dist_to_label <- dist_test_total[i,1:n1]
    ###knn
    knn_index <- order(E_dist_to_label)[1:k]
    if (length(knn_index)==1) {knn[,,i] <- train_label$Y[,,knn_index]}
    if (length(knn_index)>1) {knn[,,i] <- spd_predict(Y=train_label$Y[,,knn_index], method=method, weight=1, alpha=alpha)}
    ###nw
    loc_weight <- sapply(E_dist_to_label, function(z) K_Epa(z,h=h))
    if (length(which(loc_weight>0))==0) {nw[,,i] <- train_label$Y[,,n1]}
    if (length(which(loc_weight>0))>0) {nw[,,i] <- spd_predict(Y=train_label$Y, method=method, weight=loc_weight, alpha=alpha)}
  }
  ##testing error for supervised regerssion
  error_super <- matrix(0,n3,2)
  for (i in 1:n3){
    error_super[i,1] <- distcov(knn[,,i], test$Y[,,i], method=method, alpha=alpha)^2
    error_super[i,2] <- distcov(nw[,,i], test$Y[,,i], method=method, alpha=alpha)^2
  } 
  ISE_super <- colMeans(error_super) 
  ISE[, 3:4] <- rep(ISE_super, each = nrow(ISE))
  
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
      list_knn[[count]] <- knn_validation_spd(krange=krange, spd=train_label$Y, dis=graph_dis)
      ##leave one cross-validation for semi-nw
      diag(graph_dis) <- Inf
      h_0 <- median(rowMins(graph_dis))
      diag(graph_dis) <- 0
      hrange <- h_0*(seq(10,30,by=2))
      list_nw[[count]] <- nw_validation_spd(hrange=hrange, spd=train_label$Y, dis=graph_dis)
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
    knn_semi <- array(0,c(m,m,n3))
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test1[i,],n1), nrow = n1, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-knn
      knn_semi_index <- intersect(order(G_dist_to_label)[1:k_semi],which(G_dist_to_label<Inf))
      if (length(knn_semi_index)==0) {knn_semi[,,i] <- train_label$Y[,,n1]}
      if (length(knn_semi_index)==1) {knn_semi[,,i] <- train_label$Y[,,knn_semi_index]}
      if (length(knn_semi_index)>1) {knn_semi[,,i] <- spd_predict(Y=train_label$Y[,,knn_semi_index], method=method, weight=1, alpha=alpha)}
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
    nw_semi <- array(0,c(m,m,n3))
    for (i in 1:n3) {
      mat <- matrix(rep(dist_test2[i,],n1), nrow = n1, byrow = TRUE)
      G_dist_to_label <- rowMins(graph_dis_total+mat)
      ###semi-nw
      loc_weight_semi <- sapply(G_dist_to_label, function(z) K_Epa(z,h=h_semi))
      if (length(which(loc_weight_semi>0))==0) {nw_semi[,,i] <- train_label$Y[,,n1]}
      if (length(which(loc_weight_semi>0))>0) {nw_semi[,,i] <- spd_predict(Y=train_label$Y, method=method, weight=loc_weight_semi, alpha=alpha)}
    }
    rm(list = c("dist2","dist_test2"));gc()
    ##testing error for semi-supervised regerssion
    error_semi <- matrix(0,n3,2)
    for (i in 1:n3){
      error_semi[i,1] <- distcov(knn_semi[,,i], test$Y[,,i], method=method, alpha=alpha)^2
      error_semi[i,2] <- distcov(nw_semi[,,i], test$Y[,,i], method=method, alpha=alpha)^2
    } 
    ISE_semi <- colMeans(error_semi) 
    ISE[b,5:6] <- ISE_semi 
    b=b+1
  }
  return(ISE)
}

#----------parallel computation----------
system.time({
  cores <- min(25, parallel::detectCores() - 1)
  cl <- makeCluster(cores)
  res <- pblapply(1:20, Get_spd, cl = cl)
  stopCluster(cl)
}) 

#----------final results----------
MISE <- array(0,c(12,6,20))
for (i in 1:20) {
  MISE[,,i] <- res[[i]]
}
compare <- matrix(0,12,13)
for (i in 1:12){
  compare[i,1:6] <- rowMeans(MISE[i,,]) #mean of MSE (AMSE)
  compare[i,8:13] <- sqrt(rowVars(MISE[i,,])) #variance of MSE
} 
write.csv(x = compare,file = "spd-100-8000-1000-Curve-snr_2—2.csv")





