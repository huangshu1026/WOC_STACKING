####################################################################################
########## This file includes all previous individual weighting models  ############
########## including CWM, SS, TOPn, RP, OW, CLASSIC LASSO, RIDGE REGRESSION ########
####################################################################################

########## Input: (X, y) ###########
########## Output: weights #########

# loading the package and library
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('mgcv')) install.packages('mgcv'); library('mgcv')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('MBESS')) install.packages('MBESS'); library('MBESS')
if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc')
if (!require('quadprog')) install.packages('quadprog'); library('quadprog')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('tseries')) install.packages('tseries'); library('tseries')
if (!require('PACLasso')) install.packages('PACLasso'); library('PACLasso')

#################################################
########## Crowds Selection Methods #############
#################################################
###### Contribution Weighted Model (CWM) ######
score <- function(X, y){
  outcome <- y
  predict <- rowMeans(X)
  s <- sum((outcome - predict)^2)
  return(s)
}
contribution <- function(X, y){
  n <- nrow(X)
  s <- score(X, y)
  sj <- c()
  for(j in 1:ncol(X)){
    sj[j] <- score(X[ , -j], y)
  }
  c <- -(s - sj)/n
  return(c)
}
CWM <- function(X, y){
  M <- ncol(X)
  cont <- contribution(X, y)
  ## CEWM: equal weighting people with positive contribution
  CEWM.select <- which(cont >= 0)
  CEWM.weight <- rep(0, M)
  CEWM.weight[CEWM.select] <- 1/length(CEWM.select)
  ## CWM: weighted average of people with positive contribution
  CWM.select <- which(cont >= 0)
  CWM.weight <- rep(0, M)
  CWM.weight[CWM.select] <- cont[CWM.select]/sum(cont[CWM.select])
  
  return(list(weight.CEWM = CEWM.weight, weight.CWM = CWM.weight))
  # weight is a M-length vector
}
###### Sequential Searching Model ######
SeqSearch <- function(X, y){
  # set-up
  M <- ncol(X)
  MSE.person <- sapply(1:M, function(i) mean((X[ , i] - y)^2))
  rank.person <- rank(MSE.person) # from lowest mse to the highest mse
  ## increasing way 
  SSin.select <- c(which.min(rank.person))
  mse.train.ssin <- c(mean((X[ , SSin.select] - y)^2))
  for(i in 2:M){
    candidate <- c(1:M)[-SSin.select]
    new.mse <- sapply(1:length(candidate), function(j) mean((rowMeans(X[ , c(SSin.select , candidate[j])]) - y)^2))
    SSin.select[i] <- candidate[which.min(new.mse)]
    mse.train.ssin[i] <- min(new.mse)
  }
  SSin.select <- SSin.select[1:which.min(mse.train.ssin)]
  SSin.weight <- rep(0, M)
  SSin.weight[SSin.select] <- 1/length(SSin.select)

  ## decreasing way 
  SSde.select <- c(1:M)
  mse.train.ssde <- c()
  remove <- c()
  for(i in 1:(M-2)){
    new.mse <- sapply(1:length(SSde.select), function(j) mean((rowMeans(X[ , SSde.select[-which(SSde.select == SSde.select[j])]]) - y)^2))
    remove[i] <- SSde.select[which.min(new.mse)]
    mse.train.ssde[i] <- mean((rowMeans(X[ , -remove]) - y)^2)
    SSde.select <- c(1:M)[-remove]
  }
  new.mse1 <- mean((X[ , SSde.select[1]] - y)^2) 
  new.mse2 <- mean((X[ , SSde.select[2]] - y)^2) 
  remove[M-1] <- SSde.select[3 - which.min(c(new.mse1, new.mse2))]
  mse.train.ssde[M-1] <- min(c(new.mse1, new.mse2))
  remove <- remove[1:which.min(mse.train.ssde)]
  SSde.select <- c(1:M)[-remove]
  SSde.weight <- rep(0, M)
  SSde.weight[SSde.select] <- 1/length(SSde.select)
  
  return(list(weight.SSIN = SSin.weight, weight.SSDE = SSde.weight))
}
###### Top N Models ######
TopNModel <- function(X, y, N){
  # set-up
  M <- ncol(X)
  # selecting the top N model
  MSE.train.AllModel <- sapply(1:M, function(i) mean((X[ , i] - y)^2))
  MSEranking <- rank(MSE.train.AllModel, ties.method = "random")
  choose <- which(MSEranking <= N)
  weight <- rep(0, M)
  weight[choose] <- 1/length(choose)
  return(weight)
} 
###### Rank Performance ######
RankPerformance <- function(X, y){
  # set-up
  M <- ncol(X)
  # rank performance model
  MSE.person <- sapply(1:M, function(i) mean((X[ , i] - y)^2))
  rank.person <- rank(MSE.person, ties.method = "random") # from lowest mse to the highest mse
  MSE.RANK <- sapply(2:M, function(i) mean((rowMeans(X[ , which(rank.person <= i)]) - y)^2))
  select.num <- which.min(MSE.RANK)
  RANK.select <- which(rank.person <= select.num + 1)
  weight <- rep(0, M)
  weight[RANK.select] <- 1/length(RANK.select)
  return(weight)
} # only return the selected model's index 
###### ApproxGOpt ######
PermuteMatrix <- function(rank, Sigma, mux, sigmaxy){
  M <- length(mux)
  ##### step1: permute the matrix: Sigma + mux %*% t(mux) #####
  mat <- Sigma + mux %*% t(mux)
  index.order <- c()
  summation <- c()
  upperbound <- c()
  lowerbound <- c()
  upperbound[1] <- -2 * min(sigmaxy)
  lowerbound[1] <- -2 * max(sigmaxy)
  
  # row 1
  index.order[1] <- rank
  summation[1] <- diag(mat)[rank]
  upperbound[2] <- -2 * min(sigmaxy[-index.order])
  lowerbound[2] <- -2 * max(sigmaxy[-index.order])
  # row 2
  addition <- diag(mat) + 2 * mat[index.order, ]
  addition[index.order] <- 1000000000
  index.order[2] <- which.min(addition)
  summation[2] <- summation[1] + min(addition)
  upperbound[3] <- -2 * min(sigmaxy[-index.order])
  lowerbound[3] <- -2 * max(sigmaxy[-index.order])
  # row 3 to M-1
  for(i in 3:(M-1)){
    addition <- diag(mat) + 2 * colSums(mat[index.order, ])
    addition[index.order] <- 1000000000
    index.order[i] <- which.min(addition)
    summation[i] <- summation[i-1] + min(addition)
    upperbound[i+1] <- -2 * min(sigmaxy[-index.order])
    lowerbound[i+1] <- -2 * max(sigmaxy[-index.order])
  }
  # row M
  index.order[M] <- c(1:M)[-index.order]
  summation[M] <- sum(mat)
  diff.crit <- upperbound - lowerbound # ith element for ith person being selected 
  return(list(index.order = index.order, summation = summation, diff.crit = diff.crit))
} # return the index.order, summation, and diff.crit
ApproxGopt <- function(X, y){
  M <- ncol(X)
  Sigma <- cov(X)
  mux <- colMeans(X - matrix(rep(y, M), nrow = length(y), ncol = M))
  sigmaxy <- cor(X, y)*sqrt(diag(Sigma))*sd(y)
  sigma2y <- var(y)
  mat <- Sigma + mux %*% t(mux)
  #index.order.mat <- matrix(0, nrow = M, ncol = M) # each row is for each person 
  #summation.mat <- matrix(0, nrow = M, ncol = M)
  #diff.crit.mat <- matrix(0, nrow = M, ncol = M)
  MSE <- 100000
  SELECT <- c()
  for(m in 1:M){
    ##### step1: permute the matrix: Sigma + mux %*% t(mux) for each person being the first row #####
    permute <- PermuteMatrix(m, Sigma, mux, sigmaxy)
    index.order <- permute$index.order
    summation <- permute$summation
    diff.crit <- permute$diff.crit
    ##### step2: select crowds when first term is dominant ##### 
    select <- c()
    select[1] <- index.order[1]
    for(i in 2:M){ # adding ith person
      diff <- summation[i-1]/((i-1)^2) - summation[i]/(i^2)
      if(diff >= diff.crit[i]){
        select[i] <- index.order[i]
      } else {
        break
      }
    } # leave select to step 3
    ##### step3: adding one by one by considering all terms #####
    num <- length(select)
    if(num == 1){
      existnum <- length(select)
      candidates <- c(1:M)[-select]
      oldmse.term1 <- sum(mat[select, select])
      oldmse <- oldmse.term1/(existnum^2) - 2*sum(sigmaxy[select])/existnum
      tempmse <- (oldmse.term1 + diag(mat)[candidates] + 2 * mat[select, candidates])/((existnum+1)^2) - 2*(sum(sigmaxy[select]) + sigmaxy[candidates])/(existnum+1)
      reduction <- oldmse - tempmse
      if(max(reduction) > 0){
        select <- c(select, candidates[which.max(reduction)])
        for(j in (num + 2):(M-1)){
          existnum <- length(select)
          candidates <- c(1:M)[-select]
          oldmse.term1 <- sum(mat[select, select])
          oldmse <- oldmse.term1/(existnum^2) - 2*sum(sigmaxy[select])/existnum
          tempmse <- (oldmse.term1 + diag(mat)[candidates] + 2 * colSums(mat[select, candidates]))/((existnum+1)^2) - 2*(sum(sigmaxy[select]) + sigmaxy[candidates])/(existnum+1)
          reduction <- oldmse - tempmse
          if(max(reduction) > 0){
            select <- c(select, candidates[which.max(reduction)])
          } else {
            break
          }
        }
        MSE.new <- sum(mat[select, select])/(length(select)^2) - 2*sum(sigmaxy[select])/(length(select))
        if(MSE.new < MSE){
          SELECT <- select
          MSE <- MSE.new
        }
      } else {
        MSE.new <- sum(mat[select, select])/(length(select)^2) - 2*sum(sigmaxy[select])/(length(select))
        if(MSE.new < MSE){
          SELECT <- select
          MSE <- MSE.new
        }
      }
    } else {
      for(j in (num + 1):(M-1)){
        existnum <- length(select)
        candidates <- c(1:M)[-select]
        oldmse.term1 <- sum(mat[select, select])
        oldmse <- oldmse.term1/(existnum^2) - 2*sum(sigmaxy[select])/existnum
        tempmse <- (oldmse.term1 + diag(mat)[candidates] + 2 * colSums(mat[select, candidates]))/((existnum+1)^2) - 2*(sum(sigmaxy[select]) + sigmaxy[candidates])/(existnum+1)
        reduction <- oldmse - tempmse
        if(max(reduction) > 0){
          select <- c(select, candidates[which.max(reduction)])
        } else {
          break
        }
      }
      MSE.new <- sum(mat[select, select])/(length(select)^2) - 2*sum(sigmaxy[select])/(length(select))
      if(MSE.new < MSE){
        SELECT <- select
        MSE <- MSE.new
      }
    }
  }
  weight <- rep(0, M)
  weight[SELECT] <- 1/length(SELECT)
  return(weight)
} # for M >= 3


#################################################
####### Theoretical True Optimal Weights#########
#################################################
###### Optimal weights: with sum to one & non-negative constraints ######
optw <- function(X, y){
  M <- ncol(X)
  mat <- t(X) %*% X
  if(!is.positive.definite(mat)){
    mat <- nearPD(mat)$mat
  } # change the optimal weights 
  Rinv <- solve(chol(mat))
  C <- cbind(rep(1, M), diag(M))
  b <- c(1, rep(0, M))
  d <- t(y) %*% X
  weights <- rep(NA, M)
  tryCatch({
    qp.model <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    weights <- qp.model$solution
  }, error = function(e){
    cat("No solution for OPTW.SUM1POS!")
  })
  return(weights)
}

#################################################
#### Zero-weight prior: regularized weights #####
#################################################
### Constrained Elastic Net Regression with Zero Weights as Prior 
CalGCV.lasso <- function(X, y, lambda1){
  ## gamma = 1
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  D <- rbind(cbind(XTX, -XTX), cbind(-XTX, XTX))
  if(!is.positive.definite(D)){
    D <- nearPD(D)$mat
  }
  A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
  b0 <- c(1, rep(0, 2*M))
  d <- c(XTY, -XTY) - lambda1*rep(1, 2*M)
  result <- solve.QP(D, d, A, b0, meq = 1)$solution
  #result <- solve.QP(D, d, A, b0, meq = 1)$solution
  weights <- result[1:M] - result[(M+1):(2*M)]
  ## step 2: p(lambda1)
  # remove zero weights
  remove.index <- which(abs(weights) <= 10^(-8))
  n0 <- length(remove.index)
  if(n0 == 0){
    invW <- diag(1/(2*abs(weights)))
    B <- t(X)%*%X + lambda1*invW
    BB <- X %*% solve(B) %*% t(X)
    plambda <- sum(diag(BB)) - n0
  }
  else if(n0 >= M-1){
    plambda <- M - length(remove.index)
  }
  else if(n0 > 0 && n0 < M-1){
    invW <- diag(1/(2*abs(weights))[-remove.index])
    B <- t(X[ , -remove.index])%*%X[ , -remove.index] + lambda1*invW
    BB <- X[ , -remove.index] %*% solve(B) %*% t(X[ , -remove.index])
    plambda <- sum(diag(BB)) - n0
  }
  ## step 3: GCV
  GCV <- (t(y-X%*%weights) %*% (y-X%*%weights))/(n * (1 - plambda/n)^2)
  return(GCV)
}
CalGCV.ridge <- function(X, y, lambda2){
  ## gamma = 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  mat <- cbind(XTX + lambda2*diag(M), rep(1, M))
  mat <- rbind(mat, t(c(rep(1, M), 0)))
  result <- solve(mat) %*% c(XTY, 1)
  weights <- result[1:M]
  ## step 2: p(lambda2)
  B <- X %*% solve(t(X)%*%X + lambda2*diag(M)) %*% t(X)
  plambda <- sum(diag(B))
  ## step 3: GCV
  GCV <- (t(y-X%*%weights) %*% (y-X%*%weights))/(n * (1 - plambda/n)^2)
  return(GCV)
}
CalGCV.enet <- function(X, y, lambda1, lambda2){
  ## gamma = 1 & 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
  b0 <- c(1, rep(0, 2*M))
  d <- c(XTY, -XTY) - lambda1*rep(1, 2*M)
  D <- rbind(cbind(XTX + lambda2*diag(M), -XTX - lambda2*diag(M)), 
             cbind(-XTX - lambda2*diag(M), XTX + lambda2*diag(M)))
  if(!is.positive.definite(D)){
    D <- nearPD(D)$mat
  }
  result <- solve.QP(D, d, A, b0, meq = 1)$solution
  weights <- result[1:M] - result[(M+1):(2*M)]
  ## step 2: p(lambda1, lambda2)
  # remove zero weights
  remove.index <- which(abs(weights) <= 10^(-8))
  n0 <- length(remove.index)
  if(n0 == 0){
    invW1 <- diag(1/(2*abs(weights)))
    B <- t(X)%*%X + lambda1*invW1 + lambda2*diag(M)
    BB <- X %*% solve(B) %*% t(X)
    plambda <- sum(diag(BB)) - n0
  }
  else if(n0 >= M-1){
    plambda <- M - length(remove.index)
  }
  else if(n0 > 0 && n0 < M-1){
    invW1 <- diag(1/(2*abs(weights))[-remove.index])
    B <- t(X[ , -remove.index])%*%X[ , -remove.index] + lambda1*invW1 + lambda2*diag(M-n0)
    BB <- X[ , -remove.index] %*% solve(B) %*% t(X[ , -remove.index])
    plambda <- sum(diag(BB)) - n0
  }
  ## step 3: GCV
  GCV <- (t(y-X%*%weights) %*% (y-X%*%weights))/(n * (1 - plambda/n)^2)
  return(GCV)
}
ConstrENetZERO <- function(X, y, addl1, addl2, l1.seq, l2.seq){
  bestl1 <- NA
  bestl2 <- NA
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  if(addl1 == FALSE && addl2 == FALSE){ 
    #mat <- cbind(XTX, rep(1, M))
    #mat <- rbind(mat, t(c(rep(1, M), 0)))
    #if(!is.positive.definite(mat)){
    #  mat <- nearPD(mat)$mat
    #}
    #result <- solve(mat) %*% c(XTY, 1)
    #weights <- result[1:M]
    weights <- optweight(Sigma = cov(X)*((n-1)/n), mux = colMeans(X), sigmaxy = cov(X, y))
  } ### estimated optimal weights
  else if (addl1 == TRUE && addl2 == FALSE){ 
    gcv <- sapply(1:length(l1.seq), function(i) CalGCV.lasso(X, y, l1.seq[i]))
    bestl1 <- l1.seq[which.min(gcv)]
    D <- rbind(cbind(XTX, -XTX), cbind(-XTX, XTX))
    if(!is.positive.definite(D)){
      D <- nearPD(D)$mat
    }
    A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
    b0 <- c(1, rep(0, 2*M))
    d <- c(XTY, -XTY) - bestl1*rep(1, 2*M)
    result <- solve.QP(D, d, A, b0, meq = 1)$solution
    weights <- result[1:M] - result[(M+1):(2*M)]
  }### Lasso regression
  else if (addl1 == FALSE && addl2 == TRUE){
    gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ridge(X, y, l2.seq[i]))
    bestl2 <- l2.seq[which.min(gcv)]
    mat <- cbind(XTX + bestl2*diag(M), rep(1, M))
    mat <- rbind(mat, t(c(rep(1, M), 0)))
    result <- solve(mat) %*% c(XTY, 1)
    weights <- result[1:M]
  } ### Ridge Regression 
  else if (addl1 == TRUE && addl2 == TRUE){
    lambda.set <- cbind(sort(rep(l1.seq, length(l2.seq))), rep(l2.seq, length(l1.seq)))
    gcv <- sapply(1:(length(l1.seq)*length(l2.seq)), function(i) CalGCV.enet(X, y, lambda.set[i, 1], lambda.set[i, 2]))
    bestl1 <- lambda.set[which.min(gcv), 1]
    bestl2 <- lambda.set[which.min(gcv), 2]
    A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
    b0 <- c(1, rep(0, 2*M))
    d <- c(XTY, -XTY) - bestl1*rep(1, 2*M)
    D <- rbind(cbind(XTX + bestl2*diag(M), -XTX - bestl2*diag(M)), 
               cbind(-XTX - bestl2*diag(M), XTX + bestl2*diag(M)))
    if(!is.positive.definite(D)){
      D <- nearPD(D)$mat
    }
    result <- solve.QP(D, d, A, b0, meq = 1)$solution
    weights <- result[1:M] - result[(M+1):(2*M)]
  } ### Elastic Net regression 
  return(list(weights = weights, bestl1 = bestl1, bestl2 = bestl2))
}


#################################################
#### Equal-weight prior: regularized weights ####
#################################################
### Our Constrained Elastic Net Regression with Equal Weights as Prior 
CalGCV.ourlasso <- function(X, y, lambda1, wprior){
  ## gamma = 1
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  D <- rbind(cbind(XTX, -XTX), cbind(-XTX, XTX))
  if(!is.positive.definite(D)){
    D <- nearPD(D)$mat
  }
  A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
  b0 <- c(0, rep(0, 2*M)) # different part 
  d <- c(XTY - XTX%*%wprior, -XTY + XTX%*%wprior) - lambda1*rep(1, 2*M) # different part 
  result <- solve.QP(D, d, A, b0, meq = 1)$solution
  v <- result[1:M] - result[(M+1):(2*M)]
  ## step 2: p(lambda1)
  # remove zero v
  remove.index <- which(abs(v) <= 10^(-6))
  n0 <- length(remove.index)
  if(n0 == 0){
    invW <- diag(1/(2*abs(v)))
    B <- t(X)%*%X + lambda1*invW
    BB <- X %*% solve(B) %*% t(X)
    plambda <- sum(diag(BB)) - n0
  }
  else if(n0 >= M-1){
    plambda <- M - length(remove.index)
  }
  else if(n0 > 0 && n0 < M-1){
    invW <- diag(1/(2*abs(v))[-remove.index])
    B <- t(X[ , -remove.index])%*%X[ , -remove.index] + lambda1*invW
    BB <- X[ , -remove.index] %*% solve(B) %*% t(X[ , -remove.index])
    plambda <- sum(diag(BB)) - n0
  }
  ## step 3: GCV
  GCV <- (t(y-X%*%(v + wprior)) %*% (y-X%*%(v + wprior)))/(n * (1 - plambda/n)^2)
  return(GCV)
}
CalGCV.ourridge <- function(X, y, lambda2, wprior){
  ## gamma = 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  mat <- cbind(XTX + lambda2*diag(M), rep(1, M))
  mat <- rbind(mat, t(c(rep(1, M), 0)))
  result <- solve(mat) %*% c(XTY - XTX%*%wprior, 0) # different part 
  v <- result[1:M]  
  ## step 2: p(lambda2)
  B <- X %*% solve(t(X)%*%X + lambda2*diag(M)) %*% t(X)
  plambda <- sum(diag(B))
  ## step 3: GCV
  GCV <- (t(y-X%*%(v + wprior)) %*% (y-X%*%(v + wprior)))/(n * (1 - plambda/n)^2)
  return(GCV)
}
CalGCV.ourenet <- function(X, y, lambda1, lambda2, wprior){
  ## gamma = 1 & 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
  b0 <- c(0, rep(0, 2*M)) # different part
  d <- c(XTY - XTX%*%wprior, -XTY + XTX%*%wprior) - lambda1*rep(1, 2*M) # different part
  D <- rbind(cbind(XTX + lambda2*diag(M), -XTX - lambda2*diag(M)), 
             cbind(-XTX - lambda2*diag(M), XTX + lambda2*diag(M)))
  if(!is.positive.definite(D)){
    D <- nearPD(D)$mat
  }
  result <- solve.QP(D, d, A, b0, meq = 1)$solution
  v <- result[1:M] - result[(M+1):(2*M)]
  ## step 2: p(lambda1, lambda2)
  # remove zero v
  remove.index <- which(abs(v) <= 10^(-6))
  n0 <- length(remove.index)
  if(n0 == 0){
    invW1 <- diag(1/(2*abs(v)))
    B <- t(X)%*%X + lambda1*invW1 + lambda2*diag(M)
    BB <- X %*% solve(B) %*% t(X)
    plambda <- sum(diag(BB)) - n0
  }
  else if(n0 >= M-1){
    plambda <- M - length(remove.index)
  }
  else if(n0 > 0 && n0 < M-1){
    invW1 <- diag(1/(2*abs(v))[-remove.index])
    B <- t(X[ , -remove.index])%*%X[ , -remove.index] + lambda1*invW1 + lambda2*diag(M-n0)
    BB <- X[ , -remove.index] %*% solve(B) %*% t(X[ , -remove.index])
    plambda <- sum(diag(BB)) - n0
  }
  ## step 3: GCV
  GCV <- (t(y-X%*%(v + wprior)) %*% (y-X%*%(v + wprior)))/(n * (1 - plambda/n)^2)
  return(GCV)
}
ConstrENetEQAUL <- function(X, y, addl1, addl2, l1.seq, l2.seq, wprior){
  bestl1 <- NA
  bestl2 <- NA
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  if(addl1 == FALSE && addl2 == FALSE){
    #mat <- cbind(XTX, rep(1, M))
    #mat <- rbind(mat, t(c(rep(1, M), 0)))
    #if(!is.positive.definite(mat)){
    #  mat <- nearPD(mat)$mat
    #}
    #result <- solve(mat) %*% c(XTY, 1)
    #weights <- result[1:M]
    weights <- optweight(Sigma = cov(X)*((n-1)/n), mux = colMeans(X), sigmaxy = cov(X, y))
  } ## estimated optimal weights 
  else if (addl1 == TRUE && addl2 == FALSE){
    gcv <- sapply(1:length(l1.seq), function(i) CalGCV.ourlasso(X, y, l1.seq[i], wprior))
    bestl1 <- l1.seq[which.min(gcv)]
    D <- rbind(cbind(XTX, -XTX), cbind(-XTX, XTX))
    if(!is.positive.definite(D)){
      D <- nearPD(D)$mat
    }
    A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
    b0 <- c(0, rep(0, 2*M)) # different part 
    d <- c(XTY - XTX%*%wprior, -XTY + XTX%*%wprior) - bestl1*rep(1, 2*M) # different part 
    result <- solve.QP(D, d, A, b0, meq = 1)$solution
    weights <- result[1:M] - result[(M+1):(2*M)] + wprior # different part
  } ## our lasso regression 
  else if (addl1 == FALSE && addl2 == TRUE){
    gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ourridge(X, y, l2.seq[i], wprior))
    bestl2 <- l2.seq[which.min(gcv)]
    mat <- cbind(XTX + bestl2*diag(M), rep(1, M))
    mat <- rbind(mat, t(c(rep(1, M), 0)))
    result <- solve(mat) %*% c(XTY-XTX%*%wprior, 0) # different part 
    weights <- result[1:M] + wprior # different part 
  } ## our ridge regression 
  else if (addl1 == TRUE && addl2 == TRUE){
    lambda.set <- cbind(sort(rep(l1.seq, length(l2.seq))), rep(l2.seq, length(l1.seq)))
    gcv <- sapply(1:(length(l1.seq)*length(l2.seq)), function(i) CalGCV.ourenet(X, y, lambda.set[i, 1], lambda.set[i, 2], wprior))
    bestl1 <- lambda.set[which.min(gcv), 1]
    bestl2 <- lambda.set[which.min(gcv), 2]
    A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
    b0 <- c(0, rep(0, 2*M)) # different part
    d <- c(XTY - XTX%*%wprior, -XTY + XTX%*%wprior) - bestl1*rep(1, 2*M) # different part
    D <- rbind(cbind(XTX + bestl2*diag(M), -XTX - bestl2*diag(M)), 
               cbind(-XTX - bestl2*diag(M), XTX + bestl2*diag(M)))
    if(!is.positive.definite(D)){
      D <- nearPD(D)$mat
    }
    result <- solve.QP(D, d, A, b0, meq = 1)$solution
    weights <- result[1:M] - result[(M+1):(2*M)] + wprior # different part
  }
  return(list(weights = weights, bestl1 = bestl1, bestl2 = bestl2))
}



