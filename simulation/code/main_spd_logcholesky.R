######################################################################
# This code compares naive-knn, naive-nw, semi-knn and semi-nw for
# Setting II with SPD responses under the log-Cholesky metric.
# Select one scenario per run: "Swiss" for Figure S.3 or "Swiss2" for Figure S.4.
#   - snr: 2 and 4 for each scenario.
#   - reference_pairs: number of noise pairs used for the test means.
# Run all lines from simulation/code to save two CSVs in the working directory.
######################################################################

library(parallel)
library(pbapply)

scenario <- "Swiss"  # Figure S.3; change to "Swiss2" for Figure S.4.

#----------main function---------------
Get_spd_logcholesky <- function(loop, scenario, n1, snr, unlabeled_sizes,
                                test_size, reference_pairs) {
  source("function.R")

  set.seed(loop)
  setting <- 2
  semigroup <- unlabeled_sizes
  NN <- max(semigroup)
  n3 <- test_size
  ISE <- matrix(0, length(semigroup), 4)

  #----------data generation----------
  train_label <- generate_label_spd(scenario=scenario, setting=setting, N=n1, snr=snr)
  train_unlabel <- generate_unlabel_spd(scenario=scenario, N=NN)
  test <- generate_test_spd(scenario=scenario, setting=setting, N=n3)
  X_train <- rbind(train_label$X, train_unlabel$X)
  dist_total <- dist(X_train)
  dist_test_total <- generate_distance_two(A=test$X, B=X_train)
  coords <- t(vapply(seq_len(n1), function(i) lc_coordinates(train_label$Y[,,i]), numeric(6)))

  #----------test means under the log-Cholesky metric----------
  # test$Y contains exp(D(x)). Average the noisy responses in log-Cholesky
  # coordinates using the same noise scale as the labeled sample.
  set.seed(1000000L + loop)
  sigma <- lc_setting2_sigma(train_label$U, snr)
  reference <- lc_reference(test$Y, sigma, reference_pairs)
  testing_error <- function(prediction) {
    # Subtract the sampling variance of the Monte Carlo reference mean.
    rowSums((prediction - reference$center)^2) - reference$mc_variance
  }

  #----------supervised regression----------
  ##leave-one-out cross-validation for naive-knn and naive-nw
  dist_label <- generate_distance(A=train_label$X)
  fit <- lc_tune(dist_label, coords)
  k <- fit$k
  h <- fit$h

  ##prediction with the selected k and h
  labels <- seq_len(n1)
  E_dist_to_label <- dist_test_total[,labels,drop=FALSE]
  fallback <- rep.int(n1, n3)
  knn <- lc_knn_predictions(E_dist_to_label, coords, k, fallback)[[1]]
  nw <- lc_nw_prediction(E_dist_to_label, coords, h, fallback)
  ISE_super <- c(mean(testing_error(knn)), mean(testing_error(nw)))
  ISE[,1:2] <- rep(ISE_super, each=nrow(ISE))

  #----------semi-supervised regression----------
  for (b in seq_along(semigroup)) {
    n2 <- semigroup[b]
    n <- n1+n2

    ##cross-validation for graph radius
    candidates <- lc_radius_candidates(dist_total, n)
    rrange <- candidates$rrange
    fits <- vector("list", length(rrange))
    for (j in seq_along(rrange)) {
      g <- lc_radius_graph(candidates$pairs, n, rrange[j])
      graph_dis <- igraph::distances(g, v=labels, to=labels)
      fits[[j]] <- lc_tune(graph_dis, coords)
    }
    knn_index <- which.min(sapply(fits, function(fit) fit$k_error))
    r_knn <- rrange[knn_index]
    k_semi <- fits[[knn_index]]$k
    nw_index <- which.min(sapply(fits, function(fit) fit$h_error))
    r_nw <- rrange[nw_index]
    h_semi <- fits[[nw_index]]$h

    ##semi-supervised knn regression with the selected r and k
    g <- lc_radius_graph(candidates$pairs, n, r_knn)
    graph_dis_total <- igraph::distances(g, v=labels)
    G_dist_to_label <- lc_query_distances(graph_dis_total, dist_test_total, r_knn)
    knn_semi <- lc_knn_predictions(G_dist_to_label, coords, k_semi, fallback)[[1]]

    ##semi-supervised nw regression with the selected r and h
    if (r_nw != r_knn) {
      g <- lc_radius_graph(candidates$pairs, n, r_nw)
      graph_dis_total <- igraph::distances(g, v=labels)
      G_dist_to_label <- lc_query_distances(graph_dis_total, dist_test_total, r_nw)
    }
    nw_semi <- lc_nw_prediction(G_dist_to_label, coords, h_semi, fallback)
    ISE[b,3:4] <- c(mean(testing_error(knn_semi)), mean(testing_error(nw_semi)))
  }
  return(ISE)
}

#----------parallel computation----------
run_logcholesky <- function(configurations, unlabeled_sizes, test_size, replicates,
                            reference_pairs, workers) {
  cl <- NULL
  if (workers > 1) {
    cl <- makeCluster(min(workers, replicates))
    on.exit(stopCluster(cl))
  }

  for (j in seq_len(nrow(configurations))) {
    scenario <- configurations$scenario[j]
    n1 <- configurations$n1[j]
    snr <- configurations$snr[j]
    filename <- sprintf("spd-%d-%d-%d-%s-setting2-snr_%g_logcholesky.csv",
                        n1, max(unlabeled_sizes), test_size, scenario, snr)
    message(filename)
    res <- pblapply(seq_len(replicates), Get_spd_logcholesky, cl=cl,
                    scenario=scenario, n1=n1, snr=snr, unlabeled_sizes=unlabeled_sizes,
                    test_size=test_size, reference_pairs=reference_pairs)

    #----------final results----------
    MISE <- array(0, c(length(unlabeled_sizes), 4, replicates))
    for (i in seq_len(replicates)) {
      MISE[,,i] <- res[[i]]
    }
    compare <- matrix(0, length(unlabeled_sizes), 9)
    for (i in seq_along(unlabeled_sizes)) {
      compare[i,1:4] <- rowMeans(MISE[i,,])                  # AMSE
      compare[i,6:9] <- sqrt(matrixStats::rowVars(MISE[i,,])) # SD across replicates
    }
    pending <- tempfile("logcholesky-", tmpdir=".")
    write.csv(compare, pending)
    if (!file.rename(pending, filename)) stop("Could not save ", filename, ". Results are in ", pending)
    message("Saved ", filename)
  }
}

#----------settings for the selected scenario----------
run_logcholesky(
  configurations = data.frame(scenario = scenario, n1 = switch(scenario, Swiss=100L, Swiss2=200L),
                              snr = c(2, 4)),
  unlabeled_sizes = c(0L, 100L, 200L, 500L, 1000L, 2000L, 3000L, 4000L,
                      5000L, 6000L, 7000L, 8000L),
  test_size = 1000L,
  replicates = 100L,
  reference_pairs = 64L,  # 128 conditional noise draws per test point.
  workers = 50L           # Set to 1 for serial execution.
)
