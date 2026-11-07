######################################################################
# This script provides all utility functions required for both simulation 
# studies and real data analysis. It includes:
#   - data generation functions,  
#   - hyperparameter selection routines,  
#   - distance and metric computations,  
#   - other supporting tools.  
#
# The functions here are designed as a utility module and should be 
# imported into the main scripts for different experiments.  
######################################################################



#######################matrix exp############################
matrix.exp <- function(A){
  eig <- eigen(A)
  EA=eig$vectors%*%diag(exp(eig$values))%*%t(eig$vectors)
  return((EA+t(EA))/2)
}

#######################matrix log############################
matrix.log <- function(A){
  eig <- eigen(A)
  LA=eig$vectors%*%diag(log(eig$values))%*%t(eig$vectors)
  return((LA+t(LA))/2)
}

########################Log_p(q)###############################
log_map <- function(p, q) {
  dot <- sum(p * q)
  dot <- max(min(dot, 1), -1)  
  theta <- acos(dot)
  if (abs(theta) < 1e-10) {
    return(c(0, 0, 0))  
  }
  if (abs(pi - theta) < 1e-10) {
    v <- rnorm(3)
    v <- v - sum(v * p) * p  
    v <- v / sqrt(sum(v^2)) * pi
    return(v)
  }
  return(theta * (q - cos(theta) * p) / sin(theta))
}

########################Log_p(v)###############################
exp_map <- function(p, v) {
  norm_v <- sqrt(sum(v^2))
  if (norm_v < 1e-10) {
    return(p) 
  }
  return(cos(norm_v) * p + sin(norm_v) * (v / norm_v))
}

#######################Weierstrass function##########################
weierstrass_2d <- function(point, n_terms = 100) {
  x <- point[1]
  y <- point[2]
  sum_x <- 0
  sum_y <- 0
  for (n in 0:n_terms) {
    sum_x <- sum_x + cos(3^n * pi * x) / 2^n
    sum_y <- sum_y + cos(3^n * pi * y) / 2^n
  }
  return(sum_x + sum_y)
}


######################Epanechnikov kernel####################
K_Epa <- function(z, h = 1) {
  ifelse(is.infinite(z), 0, 3 / (4 * h) * (1 - (z / h)^2 + 10^(-5)) * (abs(z) <= h))
}

#######################calculate distance of one matrix############################
generate_distance <- function(A, method = "euclidean") {
  distance <- as.matrix(dist(A, method=method, diag=T,upper=T))
  return(distance)
}

#######################calculate cosine distance of one matrix############################
generate_nonEdistance <- function(A) {
  sim <- A / sqrt(rowSums(A * A))
  sim <- sim %*% t(sim)
  cos_distance <- 1 - sim
  return(cos_distance^2)
}

#######################calculate Euclidean distance between two matrices############################
generate_distance_two <- function(A, B) {
  n <- nrow(A); m <- nrow(B)
  distance <- matrix(rep(rowSums(A*A),m), nrow = n)+matrix(rep(rowSums(B*B),n),nrow = n, byrow=TRUE)-2*A%*%t(B)
  distance <- sqrt(distance)
  return(distance)
}

#######################generate unlabeled spd data##############################
generate_unlabel_spd <- function(scenario,N) {
  if (scenario == "Swiss") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss2") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss3") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss4") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss5") {
    U_data <- matrix(runif(N*6, 0, 1), N, 6) #U
    X_data <- U_data #X
  }
  return(list(U=U_data, X=X_data))
}

#######################generate unlabeled spherical data############################
generate_unlabel_sphere <- function(scenario,N) {
  if (scenario == "Swiss") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss2") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss3") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss4") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss5") {
    U_data <- matrix(runif(N*6, 0, 1), N, 6) #U
    X_data <- U_data #X
  }
  return(list(U=U_data, X=X_data))
}

#######################generate noised labeled spd data##############################
generate_label_spd <- function(scenario,setting,N,snr) {
  if (scenario == "Swiss") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss2") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss3") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss4") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss5") {
    U_data <- matrix(runif(N*6, 0, 1), N, 6) #U
    X_data <- U_data #X
  }
  ##generation of Y######
  if (setting ==1) {
    beta <- c(0.75, 0.25)
    M_var <- matrix(0,N,3)
    for (j in 1:N) {
      pho <- cos((as.numeric(t(beta)%*%U_data[j,]))*4*pi)
      M_var[j,] <- c(1,pho,1)
    }
    sigma <- sqrt(sum(diag(cov(M_var)))/(2.5*snr))
    Y_data=array(0,c(2,2,N))
    for(j in 1:N){
      pho <- cos((as.numeric(t(beta)%*%U_data[j,]))*4*pi)
      D <- matrix(c(1,pho,pho,1),2,2)
      Z <- matrix(0,2,2)
      Z[1,1] <- rnorm(1,0,1)
      Z[2,2] <- rnorm(1,0,1)
      Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
      logY <- sigma*Z+D
      Y_data[,,j] <-matrix.exp(logY)}}
  if (setting ==2) {
    beta_1 <- c(0.75, 0.25); beta_2 <- c(0.25, 0.75)
    M_var <- matrix(0,N,6)
    for (j in 1:N) {
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%U_data[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%U_data[j,])*4*pi)
      M_var[j,] <- c(1,pho_1,pho_2,1,pho_1,1)
    }
    sigma <- sqrt(sum(diag(cov(M_var)))/(4.5*snr))
    Y_data=array(0,c(3,3,N))
    for(j in 1:N){
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%U_data[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%U_data[j,])*4*pi)
      D <- matrix(c(1,pho_1,pho_2,pho_1,1,pho_1,pho_2,pho_1,1),3,3)
      Z <- matrix(0,3,3)
      Z[1,1] <- rnorm(1,0,1)
      Z[2,2] <- rnorm(1,0,1)
      Z[3,3] <- rnorm(1,0,1)
      Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
      Z[1,3]=Z[3,1] <- rnorm(1,0,1/sqrt(2))
      Z[2,3]=Z[3,2] <- rnorm(1,0,1/sqrt(2))
      logY <- sigma*Z+D
      Y_data[,,j] <-matrix.exp(logY)}} 
  if (setting == 3) {
    beta_1 <- c(0.1,0.2,0.3,0.4,0,0); beta_2 <- c(0,0,0.1,0.2,0.3,0.4)
    M_var <- matrix(0,N,6)
    for (j in 1:N) {
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%U_data[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%U_data[j,])*4*pi)
      M_var[j,] <- c(1,pho_1,pho_2,1,pho_1,1)
    }
    sigma <- sqrt(sum(diag(cov(M_var)))/(4.5*snr))
    Y_data=array(0,c(3,3,N))
    for(j in 1:N){
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%U_data[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%U_data[j,])*4*pi)
      D <- matrix(c(1,pho_1,pho_2,pho_1,1,pho_1,pho_2,pho_1,1),3,3)
      Z <- matrix(0,3,3)
      Z[1,1] <- rnorm(1,0,1)
      Z[2,2] <- rnorm(1,0,1)
      Z[3,3] <- rnorm(1,0,1)
      Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
      Z[1,3]=Z[3,1] <- rnorm(1,0,1/sqrt(2))
      Z[2,3]=Z[3,2] <- rnorm(1,0,1/sqrt(2))
      logY <- sigma*Z+D
      Y_data[,,j] <-matrix.exp(logY)}}
  if (setting ==4) {
    M_var <- matrix(0,N,3)
    for (j in 1:N) {
      pho <- weierstrass_2d(U_data[j,])
      M_var[j,] <- c(1,pho,1)
    }
    sigma <- sqrt(sum(diag(cov(M_var)))/(2.5*snr))
    Y_data=array(0,c(2,2,N))
    for(j in 1:N){
      pho <- weierstrass_2d(U_data[j,])
      D <- matrix(c(1,pho,pho,1),2,2)
      Z <- matrix(0,2,2)
      Z[1,1] <- rnorm(1,0,1)
      Z[2,2] <- rnorm(1,0,1)
      Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
      logY <- sigma*Z+D
      Y_data[,,j] <-matrix.exp(logY)}}
  train_spd <- list(U=U_data, X=X_data, Y=Y_data)
  return(train_spd)
}

#######################generate noised labeled spherical data############################
generate_label_sphere <- function(scenario, N, kappa, sigma) {
  if (scenario == "Swiss") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss2") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss3") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss4") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss5") {
    U_data <- matrix(runif(N*6, 0, 1), N, 6) #U
    X_data <- U_data #X
  }
  ##generation of Y###############
  Y_true <- matrix(0,N,3)
  if (scenario != "Swiss5") {
    for (i in 1:N){
      tmpx <- sin(pi/2 * (1 - U_data[i,1])) * cos(pi * U_data[i,2])
      tmpy <- sin(pi/2 * (1 - U_data[i,1])) * sin(pi * U_data[i,2])
      tmpz <- cos(pi/2 * (1 - U_data[i,1]))
      Y_true[i,] <- c(tmpx,tmpy,tmpz)}
  }else{
    beta_1 <- c(0.1,0.2,0.3,0.4,0,0); beta_2 <- c(0,0,0.1,0.2,0.3,0.4)
    for (i in 1:N){
      tmpx <- sin(pi/2 * t(beta_1)%*%U_data[i,])*cos(pi * t(beta_2)%*%U_data[i,])
      tmpy <- sin(pi/2 * t(beta_1)%*%U_data[i,])*sin(pi * t(beta_2)%*%U_data[i,])
      tmpz <- cos(pi/2 * t(beta_1)%*%U_data[i,])
      Y_true[i,] <- c(tmpx,tmpy,tmpz)}
  }
  Y_orig <- t(sapply(1:N, function(i) {
    rvmf(1, mu = Y_true[i,], k = kappa)
  }))
  Y_data <- matrix(0, N, 3)
  if (kappa==0) {alpha <- sqrt(2*sigma^2/(pi^2-4))}
  if (kappa>0) {
    norconst <- kappa/(4*pi*sinh(kappa))
    integrand <- function(theta) {
      theta^2 * sin(theta) * exp(kappa * cos(theta))
    }
    intconst <- integrate(integrand, lower = 0, upper = pi)$value
    alpha <- sqrt(sigma^2/(2*pi*norconst*intconst))
  }
  for (i in 1:N) {
    v <- log_map(Y_true[i,], Y_orig[i, ])
    Y_data[i, ] <- exp_map(Y_true[i,], alpha * v)
  }
  train_sphere <- list(U=U_data, X=X_data, Y=Y_data)
  return(train_sphere)
}

#######################generate true labeled spd (testing data)##############################
generate_test_spd <- function(scenario,setting,N) {
  if (scenario == "Swiss") {  
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss2") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss3") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss4") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss5") {
    U_data <- matrix(runif(N*6, 0, 1), N, 6) #U
    X_data <- U_data #X
  }
  ##Y
  if (setting ==1) {
    Y_data=array(0,c(2,2,N))
    beta <- c(0.75, 0.25)
    for(j in 1:N){
      pho <- cos((as.numeric(t(beta)%*%U_data[j,]))*4*pi)
      Y_data[,,j] <- matrix.exp(matrix(c(1,pho,pho,1),2,2))}}
  if (setting ==2) {
    Y_data=array(0,c(3,3,N))
    beta_1 <- c(0.75, 0.25); beta_2 <- c(0.25, 0.75)
    for(j in 1:N){
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%U_data[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%U_data[j,])*4*pi)
      Y_data[,,j] <-matrix.exp(matrix(c(1,pho_1,pho_2,pho_1,1,pho_1,pho_2,pho_1,1),3,3))}} 
  if (setting ==3) {
    Y_data=array(0,c(3,3,N))
    beta_1 <- c(0.1,0.2,0.3,0.4,0,0); beta_2 <- c(0,0,0.1,0.2,0.3,0.4)
    for(j in 1:N){
      pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%U_data[j,])*4*pi)
      pho_2 <- 0.4*cos(as.numeric(t(beta_2)%*%U_data[j,])*4*pi)
      Y_data[,,j] <-matrix.exp(matrix(c(1,pho_1,pho_2,pho_1,1,pho_1,pho_2,pho_1,1),3,3))}} 
  if (setting ==4) {
    Y_data=array(0,c(2,2,N))
    for(j in 1:N){
      pho <- weierstrass_2d(U_data[j,])
      Y_data[,,j] <- matrix.exp(matrix(c(1,pho,pho,1),2,2))}}
  test_spd <- list(U=U_data, X=X_data, Y=Y_data)
  return(test_spd)
}

#######################generate true labeled spherical data (testing data)############################
generate_test_sphere <- function(scenario, N) {
  if (scenario == "Swiss") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss2") {
    U_data <- matrix(runif(N*2, 0, 1), N, 2) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss3") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta <- 4*pi*(U_data[,1]+1/2)
    X_data <- matrix(0, N, 3) #X
    X_data[,1] <- 0.1*theta*cos(theta); X_data[,2] <- 4*U_data[,2]; X_data[,3] <- 0.1*theta*sin(theta)}
  if (scenario == "Swiss4") {
    mu <- c(0.5,0.5)
    sig <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
    U_data <- rtmvnorm(N, mean = mu, sigma = sig, lower = c(0,0), upper = c(1,1)) #U
    theta1 <- 4*pi*(U_data[,1]+1/2)
    theta2 <- 4*pi*(U_data[,2]+1/2)
    X_data <- matrix(0, N, 6) #X
    X_data[,1] <- 0.1*theta1*cos(theta1); X_data[,2] <- U_data[,2]; X_data[,3] <- 0.1*theta1*sin(theta1)
    X_data[,4] <- 0.1*theta2*cos(theta2); X_data[,5] <- -U_data[,1]; X_data[,6] <- 0.1*theta2*sin(theta2)
  }
  if (scenario == "Swiss5") {
    U_data <- matrix(runif(N*6, 0, 1), N, 6) #U
    X_data <- U_data #X
  }
  ##Y
  Y_data <- matrix(0,N,3)
  if (scenario != "Swiss5") {
    for (i in 1:N){
      tmpx <- sin(pi/2 * (1 - U_data[i,1])) * cos(pi * U_data[i,2])
      tmpy <- sin(pi/2 * (1 - U_data[i,1])) * sin(pi * U_data[i,2])
      tmpz <- cos(pi/2 * (1 - U_data[i,1]))
      Y_data[i,] <- c(tmpx,tmpy,tmpz)}
  }else{
    beta_1 <- c(0.1,0.2,0.3,0.4,0,0); beta_2 <- c(0,0,0.1,0.2,0.3,0.4)
    for (i in 1:N){
      tmpx <- sin(pi/2 * t(beta_1)%*%U_data[i,])*cos(pi * t(beta_2)%*%U_data[i,])
      tmpy <- sin(pi/2 * t(beta_1)%*%U_data[i,])*sin(pi * t(beta_2)%*%U_data[i,])
      tmpz <- cos(pi/2 * t(beta_1)%*%U_data[i,])
      Y_data[i,] <- c(tmpx,tmpy,tmpz)}
  }
  text_sphere <- list(U=U_data, X=X_data, Y=Y_data)
  return(text_sphere)
}

########spd prediction function##########
spd_predict <- function(Y, method=method, weight, alpha=alpha){
  result <- estcov(Y, method=method, weights=weight, alpha=alpha)$mean
  return(result)
}


#########cross validation for knn method of spd ########################
knn_validation_spd <- function(krange, spd, dis){
  n1 <- min(dim(dis)[1],200)
  num <- length(krange)
  record <- rep(0,num)
  t <- 1
  for (k in krange) {
    err <- 0
    for (i in 1:n1) {
      E_dis_to_label <- dis[i,]
      knn_index <- intersect(order(E_dis_to_label)[2:(k+1)],which(E_dis_to_label<Inf))
      if (length(knn_index)==0){knn_i <- spd[,,-i][,,1]}
      if (length(knn_index)==1){knn_i <- spd[,,knn_index]}
      if (length(knn_index)>1){knn_i <- spd_predict(Y=spd[,,knn_index], method=method, weight=1, alpha=alpha)}
      err <- err+distcov(knn_i, spd[,,i], method=method, alpha=alpha)^2
    }
    record[t] <- err
    t=t+1
  }
  return(list(err_opt=min(record),k_opt=krange[which.min(record)]))
}

#########cross validation for nw method of spd ########################
nw_validation_spd <- function(hrange, spd, dis){
  n1 <- min(dim(dis)[1],200)
  num <- length(hrange)
  record <- rep(0,num)
  t <- 1
  for (h in hrange) {
    err <- 0
    for (i in 1:n1) {
      E_dis_to_label <- dis[i,]
      loc_weight <- sapply(E_dis_to_label, function(z) K_Epa(z,h=h))
      if (length(which(loc_weight>0))<=1) {nw_i <- spd[,,-i][,,1]}
      if (length(which(loc_weight>0))>1) {nw_i <- spd_predict(Y=spd[,,-i], method=method, weight=loc_weight[-i], alpha=alpha)}
      err <- err+distcov(nw_i, spd[,,i], method=method, alpha=alpha)^2
    }
    record[t] <- err
    t=t+1
  }
  return(list(err_opt=min(record),h_opt=hrange[which.min(record)]))
}


#########cross validation for knn method of sphere ########################
knn_validation_sphere <- function(krange, sphere, dis){
  n1 <- min(dim(dis)[1],200)
  num <- length(krange)
  record <- rep(0,num)
  t <- 1
  for (k in krange) {
    err <- 0
    for (i in 1:n1) {
      E_dis_to_label <- dis[i,]
      knn_index <- intersect(order(E_dis_to_label)[2:(k+1)],which(E_dis_to_label<Inf))
      if (length(knn_index)==0){knn_i <- sphere[-i,][1,]}
      if (length(knn_index)==1){knn_i <- sphere[knn_index,]}
      if (length(knn_index)>1){
        myriem = wrap.sphere(sphere[knn_index,])
        knn_i <- c(riem.mean(myriem, weight=rep(1,length(knn_index)), geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
      cross <- crossprod(knn_i,sphere[i,])[1,1]
      err <- err+acos(cross)^2
    }
    record[t] <- err
    t=t+1
  }
  return(list(err_opt=min(record),k_opt=krange[which.min(record)]))
}

#########cross validation for nw method of sphere ########################
nw_validation_sphere <- function(hrange, sphere, dis){
  n1 <- min(dim(dis)[1],200)
  num <- length(hrange)
  record <- rep(0,num)
  t <- 1
  for (h in hrange) {
    err <- 0
    for (i in 1:n1) {
      E_dis_to_label <- dis[i,]
      loc_weight <- sapply(E_dis_to_label, function(z) K_Epa(z,h=h))
      if (length(which(loc_weight>0))<=1) {nw_i <- sphere[-i,][1,]}
      if (length(which(loc_weight>0))>1) {
        myriem = wrap.sphere(sphere[-i,])
        nw_i <- c(riem.mean(myriem, weight=loc_weight[-i], geometry="intrinsic", maxiter=100, eps=1e-5)$mean)}
      cross <- crossprod(nw_i,sphere[i,])[1,1]
      err <- err+acos(cross)^2
    }
    record[t] <- err
    t=t+1
  }
  return(list(err_opt=min(record),h_opt=hrange[which.min(record)]))
}

#########cross validation for knn method of Euclidean data ########################
knn_validation_Euc <- function(krange, Euc, dis){
  n1 <- min(dim(dis)[1],200)
  num <- length(krange)
  record <- rep(0,num)
  t <- 1
  for (k in krange) {
    err <- 0
    for (i in 1:n1) {
      E_dis_to_label <- dis[i,]
      knn_index <- intersect(order(E_dis_to_label)[2:(k+1)],which(E_dis_to_label<Inf))
      if (length(knn_index)==0){knn_i <- Euc[-i][1]}
      if (length(knn_index)==1){knn_i <- Euc[knn_index]}
      if (length(knn_index)>1){knn_i <- mean(Euc[knn_index])}
      err <- err+abs(knn_i-Euc[i])
    }
    record[t] <- err
    t <- t + 1
  }
  return(list(err_opt=min(record),k_opt=krange[which.min(record)]))
}

#########cross validation for nw method of Euclidean data ########################
nw_validation_Euc <- function(hrange, Euc, dis){
  n1 <- min(dim(dis)[1],200)
  num <- length(hrange)
  record <- rep(0,num)
  t <- 1
  for (h in hrange) {
    err <- 0
    for (i in 1:n1) {
      E_dis_to_label <- dis[i,]
      loc_weight <- sapply(E_dis_to_label, function(z) K_Epa(z,h=h))
      if (length(which(loc_weight>0))==1) {nw_i <- Euc[-i][1]}
      if (length(which(loc_weight>0))>1) {nw_i <- sum(Euc[-i]*loc_weight[-i])/sum(loc_weight[-i])}
      err <- err+abs(nw_i-Euc[i])
    }
    record[t] <- err
    t=t+1
  }
  return(list(err_opt=min(record),h_opt=hrange[which.min(record)]))
}

######### Estimate eta and local curvatures ###########################
generate_quadratic_terms <- function(Y) {
  n <- nrow(Y)
  d <- ncol(Y)
  terms <- NULL
  for (i in 1:d) {
    terms <- cbind(terms, Y[, i]^2)
  }
  if (d > 1) {
    for (i in 1:(d-1)) {
      for (j in (i+1):d) {
        terms <- cbind(terms, 2 * Y[, i] * Y[, j])
      }
    }
  }
  return(terms)
}

estimate_eta_with_curvature <- function(X, d, k = 30) {
  n <- nrow(X)
  D <- ncol(X)
  eta_vals <- numeric(n)
  
  for (i in 1:n) {
    nn_idx <- get.knnx(X, X[i, , drop = FALSE], k = k)$nn.index
    neighbors <- X[nn_idx, , drop = FALSE]
    
    center <- colMeans(neighbors)
    neighbors_centered <- sweep(neighbors, 2, center)
    
    svd_res <- svd(neighbors_centered)
    tangent_basis <- svd_res$v[, 1:d, drop = FALSE]
    normal_basis <- svd_res$v[, (d + 1):D, drop = FALSE]
    
    Y_tangent <- neighbors_centered %*% tangent_basis
    Z_normal <- neighbors_centered %*% normal_basis
    
    curvatures <- c()
    for (j in 1:ncol(Z_normal)) {
      z <- Z_normal[, j]
      Phi <- generate_quadratic_terms(Y_tangent)
      beta <- lm(z ~ Phi - 1)$coefficients
      
      H <- matrix(0, d, d)
      diag(H) <- beta[1:d]
      idx <- d + 1
      for (p1 in 1:(d-1)) {
        for (p2 in (p1+1):d) {
          H[p1, p2] <- beta[idx]
          H[p2, p1] <- beta[idx]
          idx <- idx + 1
        }
      }
      eigvals <- svd(H)$d
      curvatures <- c(curvatures, eigvals)
    }
    eta_vals[i] <- max(abs(curvatures))
  }
  return(list(global_eta = max(eta_vals), local_eta = eta_vals))
}


######### log-cholesky metric ###########################
distLogCholesky <- function(A, B) {
  M.a <- chol(A)
  D.a <- diag(M.a)
  L.a <- M.a - diag(D.a)
  M.b <- chol(B)
  D.b <- diag(M.b)
  L.b <- M.b - diag(D.b) 
  sqrt(sum((L.a-L.b)^2)+sum((log(D.a)-log(D.b))^2)) 
}

######### weighted freshet mean for spd ###########################
estLogCholesky <- function(S, weights = 1) {
  M <- dim(S)[3]
  if (length(weights) == 1) {
    weights <- rep(1, times = M)
  }
  sum <- S[, , 1] * 0
  for (j in 1:M) {
    C <- chol(S[, , j])
    D <- diag(C)
    L <- C - diag(D)
    sum <- sum + t(diag(log(D))+L) * weights[j] / sum(weights)
  }
  cc <- sum
  D <- diag(cc)
  L <- cc - diag(D)
  cc <- diag(exp(D))+L
  cc %*% t(cc)
}

estcov <-
  function (S,
            method = "Riemannian",
            weights = 1,
            alpha = 1 / 2,
            MDSk = 2)
  {
    out <- list(
      mean = 0,
      sd = 0,
      pco = 0,
      eig = 0,
      dist = 0
    )
    M <- dim(S)[3]
    if (length(weights) == 1) {
      weights <- rep(1, times = M)
    }
    if (method == "Procrustes") {
      dd <- estSS(S, weights)
    }
    if (method == "ProcrustesShape") {
      dd <- estShape(S, weights)
    }
    if (method == "Riemannian") {
      dd <- estLogRiem2(S, weights)
    }
    if (method == "Cholesky") {
      dd <- estCholesky(S, weights)
    }
    if (method == "LogCholesky") {
      dd <- estLogCholesky(S, weights)
    }
    if (method == "Power") {
      dd <- estPowerEuclid(S, weights, alpha)
    }
    if (method == "Euclidean") {
      dd <- estEuclid(S, weights)
    }
    if (method == "LogEuclidean") {
      dd <- estLogEuclid(S, weights)
    }
    if (method == "RiemannianLe") {
      dd <- estRiemLe(S, weights)
    }
    out$mean <- dd
    out
  }

distcov <- function(S1,
                    S2 ,
                    method = "Riemannian",
                    alpha = 1 / 2) {
  if (method == "Procrustes") {
    dd <- distProcrustesSizeShape(S1, S2)
  }
  if (method == "ProcrustesShape") {
    dd <- distProcrustesFull(S1, S2)
  }
  if (method == "Riemannian") {
    dd <- distRiemPennec(S1, S2)
  }
  if (method == "Cholesky") {
    dd <- distCholesky(S1, S2)
  }
  if (method == "LogCholesky") {
    dd <- distLogCholesky(S1, S2)
  }
  if (method == "Power") {
    dd <- distPowerEuclidean(S1, S2, alpha)
  }
  if (method == "Euclidean") {
    dd <- distEuclidean(S1, S2)
  }
  if (method == "LogEuclidean") {
    dd <- distLogEuclidean(S1, S2)
  }
  if (method == "RiemannianLe") {
    dd <- distRiemannianLe(S1, S2)
  }
  dd
}


#######################generate unlabeled spd data##############################
generate_unlabel_spd_curve <- function(N, a=6*pi, b=0.5) {
  U_data <- runif(N, 0, 1) #U
  X_data <- matrix(0, N, 3) #X
  X_data[,1] <- cos(a*U_data); X_data[,2] <- sin(a*U_data); X_data[,3] <- b*U_data
  return(list(U=U_data, X=X_data))
}

#######################generate noised labeled spd data##############################
generate_label_spd_curve <- function(N, snr) {
  data <- generate_unlabel_spd_curve(N)
  U_data <- data$U
  X_data <- data$X 
  ##generation of Y######
  M_var <- matrix(0,N,3)
  for (j in 1:N) {
    pho <- cos(U_data[j]*4*pi)
    M_var[j,] <- c(1,pho,1)
  }
  sigma <- sqrt(sum(diag(cov(M_var)))/(2.5*snr))
  Y_data=array(0,c(2,2,N))
  for(j in 1:N){
    pho <- cos(U_data[j]*4*pi)
    D <- matrix(c(1,pho,pho,1),2,2)
    Z <- matrix(0,2,2)
    Z[1,1] <- rnorm(1,0,1)
    Z[2,2] <- rnorm(1,0,1)
    Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
    logY <- sigma*Z+D
    Y_data[,,j] <-matrix.exp(logY)}
  train_spd <- list(U=U_data, X=X_data, Y=Y_data)
  return(train_spd)
}

#######################generate true labeled spd (testing data)##############################
generate_test_spd_curve <- function(N) {
  data <- generate_unlabel_spd_curve(N)
  U_data <- data$U
  X_data <- data$X 
  ##Y
  Y_data=array(0,c(2,2,N))
  for(j in 1:N){
    pho <- cos(U_data[j]*4*pi)
    Y_data[,,j] <- matrix.exp(matrix(c(1,pho,pho,1),2,2))}
  test_spd <- list(U=U_data, X=X_data, Y=Y_data)
  return(test_spd)
}

#######################generate unlabeled spherical data############################
generate_unlabel_sphere_curve <- function(N, a=6*pi, b=0.5) {
  U_data <- runif(N, 0, 1) #U
  X_data <- matrix(0, N, 3) #X
  X_data[,1] <- cos(a*U_data); X_data[,2] <- sin(a*U_data); X_data[,3] <- b*U_data
  return(list(U=U_data, X=X_data))
}

#######################generate noised labeled spherical data############################
generate_label_sphere_curve <- function(N, kappa, sigma) {
  data <- generate_unlabel_sphere_curve(N)
  U_data <- data$U
  X_data <- data$X 
  ##generation of Y###############
  Y_true <- matrix(0,N,3)
  for (i in 1:N){
    tmpx <- sin(pi/2 * U_data[i]) * cos(pi * U_data[i])
    tmpy <- sin(pi/2 * U_data[i]) * sin(pi * U_data[i])
    tmpz <- cos(pi/2 * U_data[i])
    Y_true[i,] <- c(tmpx,tmpy,tmpz)}
  Y_orig <- t(sapply(1:N, function(i) {
    rvmf(1, mu = Y_true[i,], k = kappa)
  }))
  Y_data <- matrix(0, N, 3)
  if (kappa==0) {alpha <- sqrt(2*sigma^2/(pi^2-4))}
  if (kappa>0) {
    norconst <- kappa/(4*pi*sinh(kappa))
    integrand <- function(theta) {
      theta^2 * sin(theta) * exp(kappa * cos(theta))
    }
    intconst <- integrate(integrand, lower = 0, upper = pi)$value
    alpha <- sqrt(sigma^2/(2*pi*norconst*intconst))
  }
  for (i in 1:N) {
    v <- log_map(Y_true[i,], Y_orig[i, ])
    Y_data[i, ] <- exp_map(Y_true[i,], alpha * v)
  }
  train_sphere <- list(U=U_data, X=X_data, Y=Y_data)
  return(train_sphere)
}

#######################generate true labeled spherical data (testing data)############################
generate_test_sphere_curve <- function(N) {
  data <- generate_unlabel_sphere_curve(N)
  U_data <- data$U
  X_data <- data$X 
  ##Y
  Y_data <- matrix(0,N,3)
  for (i in 1:N){
    tmpx <- sin(pi/2 * U_data[i]) * cos(pi * U_data[i])
    tmpy <- sin(pi/2 * U_data[i]) * sin(pi * U_data[i])
    tmpz <- cos(pi/2 * U_data[i])
    Y_data[i,] <- c(tmpx,tmpy,tmpz)}
  text_sphere <- list(U=U_data, X=X_data, Y=Y_data)
  return(text_sphere)
}