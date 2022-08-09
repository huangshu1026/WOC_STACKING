##############################################################################################
########## This file includes function to get all weights and the stacking method  ###########
##############################################################################################
#### function to calculate different weights ####
WeightsCalculation <- function(X, y, weight_index, w5.select = c(5:8), w6.select = c(9:10)){
  Allweights <- c()
  if(1 %in% weight_index){
    w1.sa <- eqw(X, y)
    Allweights <- cbind(Allweights, w1.sa)
  }
  if(2 %in% weight_index){
    w2.optw_bias <- optw_bias(X, y)
    Allweights <- cbind(Allweights, w2.optw_bias)
  }
  if(3 %in% weight_index){
    w3.model <- optw_cov(X, y)
    w3.optw_cov  <- w3.model$weight
    bias <- w3.model$bias
    Allweights <- cbind(Allweights, w3.optw_cov)
  }
  if(4 %in% weight_index){
    w4.model <- optw_var(X, y)
    w4.optw_var  <- w4.model$weight
    bias <- w4.model$bias
    Allweights <- cbind(Allweights, w4.optw_var)
  }
  if(5 %in% weight_index & 5 %in% w5.select){
    w5.model.exo <- optw_var_1cor_exo(X, y)
    w5.optw_var_1cor_exo <- w5.model.exo$weight
    Allweights <- cbind(Allweights, w5.optw_var_1cor_exo)
  }
  if(5 %in% weight_index & 6 %in% w5.select){
    w5.model.bme <- optw_var_1cor_BayesMean(X, y)
    w5.optw_var_1cor_BME <- w5.model.bme$weight
    Allweights <- cbind(Allweights, w5.optw_var_1cor_BME)
  }
  if(5 %in% weight_index & 7 %in% w5.select){
    w5.model.bmap <- optw_var_1cor_BayesMAP(X, y)
    w5.optw_var_1cor_BMAP <- w5.model.bmap$weight
    Allweights <- cbind(Allweights, w5.optw_var_1cor_BMAP)
  }
  if(5 %in% weight_index & 8 %in% w5.select){
    w5.model.mean <- optw_var_1cor_mean(X, y)
    w5.optw_var_1cor_mean <- w5.model.mean$weight
    Allweights <- cbind(Allweights, w5.optw_var_1cor_mean)
  }
  if(6 %in% weight_index & 9 %in% w6.select){
    w6.OurConRdg <- rep(NA, ncol(X))
    tryCatch({
      w6.OurConRdg <- ConsRidge(X, y)
    }, error = function(e){
      #cat("no consridge.")
    })
    Allweights <- cbind(Allweights, w6.OurConRdg)
  }
  if(6 %in% weight_index & 10 %in% w6.select){
    w6.OurConLas <- rep(NA, ncol(X))
    tryCatch({
      w6.OurConLas <- ConsLas(X, y)
    }, error = function(e){
      #cat("no conslas.")
    })
    Allweights <- cbind(Allweights, w6.OurConLas)
  }
  if(7 %in% weight_index){
    w7.cwm <- CWM(X, y)$weight.CWM 
    Allweights <- cbind(Allweights, w7.cwm)
  }
  if(8 %in% weight_index){
    w8.ssin <- SeqSearch(X, y)$weight.SSIN 
    Allweights <- cbind(Allweights, w8.ssin)
  }
  if(9 %in% weight_index){
    w9.ssde <- SeqSearch(X, y)$weight.SSDE
    Allweights <- cbind(Allweights, w9.ssde) 
  }
  if(10 %in% weight_index){
    w10.rp <- RankPerformance(X, y)
    Allweights <- cbind(Allweights, w10.rp)
  }
  if(11 %in% weight_index){
    w11.top3 <- TopNModel(X, y, 3)
    #w11.top5 <- TopNModel(X, y, 5)
    Allweights <- cbind(Allweights, w11.top3)
  }
  if(12 %in% weight_index){
    w12.vanilla <- VanillaStacking(X, y)
    Allweights <- cbind(Allweights, w12.vanilla)
  }
  return(list(weights = Allweights, bias = bias))
}
#### prediction from individual models ####
MakePrediction <- function(X.test, weights, debias){
  # with bias 
  pred1 <- X.test %*% weights
  # without bias
  debias.mat <- matrix(rep(debias, each = nrow(X.test)), ncol = ncol(X.test), nrow = nrow(X.test))
  pred2 <- (X.test - debias.mat) %*% weights
  pred.all <- cbind(pred1[ , c(1:2)], pred2[ , 3:5], pred1[ , 6:12])
  return(pred.all)
}
#### meta learners ####
# for the second meta learner
ConstrRidge_sum1 <- function(X, y, l.min = -3, l.max = 6, step = 0.1){
  l2.seq <- 10^seq(l.min, l.max, by = step)
  bestl2 <- NA
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  
  gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ridge_sum1(X, y, l2.seq[i]))
  if(all(is.na(gcv) == TRUE)){
    l2.seq.new <- 10^seq(l.min-3, l.min, by = 0.1)
    gcv <- sapply(1:length(l2.seq.new), function(i) CalGCV.ridge_sum1(X, y, l2.seq.new[i]))
    bestl2 <- l2.seq.new[which.min(gcv)]
  } else {
    if(all(is.na(gcv) == TRUE)){
      bestl2 <- 0
    } else {
      bestl2 <- l2.seq[which.min(gcv)]
    }
  }
  H <- cbind(XTX + bestl2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(1)
  results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq)), bvec = c(b.eq), meq = 1)
  weights <- results$solution
  return(weights)
}
CalGCV.ridge_sum1 <- function(X, y, lambda2){
  ## gamma = 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  H <- cbind(XTX + lambda2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(1)
  GCV <- NA
  tryCatch({
    results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq)), bvec = c(b.eq), meq = 1)
    weights <- results$solution
    ## step 2: p(lambda2)
    B <- X %*% solve(t(X)%*%X + lambda2*diag(M)) %*% t(X)
    plambda <- sum(diag(B))
    ## step 3: GCV
    GCV <- (t(y-X%*%weights) %*% (y-X%*%weights))/(n * (1 - plambda/n)^2)
  }, error = function(e){
    NULL
  })
  return(GCV)
}
# for the third meta learner
ConstrRidge <- function(X, y, l.min = -3, l.max = 6, step = 0.1){
  l2.seq <- 10^seq(l.min, l.max, by = step)
  bestl2 <- NA
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  
  gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ridge(X, y, l2.seq[i]))
  if(all(is.na(gcv) == TRUE)){
    l2.seq.new <- 10^seq(l.min-3, l.min, by = 0.1)
    gcv <- sapply(1:length(l2.seq.new), function(i) CalGCV.ridge(X, y, l2.seq.new[i]))
    bestl2 <- l2.seq.new[which.min(gcv)]
  } else {
    if(all(is.na(gcv) == TRUE)){
      bestl2 <- 0
    } else {
      bestl2 <- l2.seq[which.min(gcv)]
    }
  }
  H <- cbind(XTX + bestl2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(1)
  # inequalities
  A.ge <- diag(M)
  b.ge <- rep(0, M)
  results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq, A.ge)), bvec = c(b.eq, b.ge), meq = 1)
  weights <- results$solution
  return(weights)
}
CalGCV.ridge <- function(X, y, lambda2){
  ## gamma = 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  H <- cbind(XTX + lambda2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(1)
  # inequalities
  A.ge <- diag(M)
  b.ge <- rep(0, M)
  GCV <- NA
  tryCatch({
    results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq, A.ge)), bvec = c(b.eq, b.ge), meq = 1)
    weights <- results$solution
    ## step 2: p(lambda2)
    B <- X %*% solve(t(X)%*%X + lambda2*diag(M)) %*% t(X)
    plambda <- sum(diag(B))
    ## step 3: GCV
    GCV <- (t(y-X%*%weights) %*% (y-X%*%weights))/(n * (1 - plambda/n)^2)
  }, error = function(e){
    NULL
  })
  return(GCV)
}
# for the last meta leaner
ConstrRidge.priorw <- function(X, y, l1.seq, wprior){
  M <- ncol(X)
  w.prime <- wprior
  y.new <- y - X %*% w.prime
  model <- ConsRidge.priorw(X, y.new)
  weight <- model + w.prime
  return(weight)
}
ConsRidge.priorw <- function(X, y, l.min = -3, l.max = 6, step = 0.1){
  l2.seq <- 10^seq(l.min, l.max, by = step)
  bestl2 <- NA
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  
  gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ridge.priorw(X, y, l2.seq[i]))
  if(all(is.na(gcv) == TRUE)){
    l2.seq.new <- 10^seq(l.min-3, l.min, by = 0.1)
    gcv <- sapply(1:length(l2.seq.new), function(i) CalGCV.ridge.priorw(X, y, l2.seq.new[i]))
    bestl2 <- l2.seq.new[which.min(gcv)]
  } else {
    if(all(is.na(gcv) == TRUE)){
      bestl2 <- 0
    } else {
      bestl2 <- l2.seq[which.min(gcv)]
    }
  }
  H <- cbind(XTX + bestl2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(0)
  # inequalities
  A.ge <- diag(M)
  b.ge <- c(-1, rep(0, M-1))
  results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq, A.ge)), bvec = c(b.eq, b.ge), meq = 1)
  weights <- results$solution
  return(weights)
}
CalGCV.ridge.priorw <- function(X, y, lambda2){
  ## gamma = 2
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  ## step 1: calculate beta hat
  H <- cbind(XTX + lambda2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(0)
  # inequalities
  A.ge <- diag(M)
  b.ge <- c(-1, rep(0, M-1))
  GCV <- NA
  tryCatch({
    results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq, A.ge)), bvec = c(b.eq, b.ge), meq = 1)
    weights <- results$solution
    ## step 2: p(lambda2)
    B <- X %*% solve(t(X)%*%X + lambda2*diag(M)) %*% t(X)
    plambda <- sum(diag(B))
    ## step 3: GCV
    GCV <- (t(y-X%*%weights) %*% (y-X%*%weights))/(n * (1 - plambda/n)^2)
  }, error = function(e){
    NULL
  })
  return(GCV)
}


MetaAlgorithm <- function(Z.mat, y.vec, meta_index){
  all_coef <- c()
  if(1 %in% meta_index){
    ## conventional Lasso 
    coef.stacking1 <- rep(NA, 1 + ncol(Z.mat))
    tryCatch({
      lambda_seq <- 10^seq(-3, 6, by = 0.1)
      stacking <- cv.glmnet(Z.mat, y.vec, alpha = 1, lambda = lambda_seq)
      best_lambda <- stacking$lambda.min
      stacking <- glmnet(Z.mat, y.vec, alpha = 1, lambda = best_lambda)
      coef.stacking1 <- as.vector(coef(stacking))
    }, error = function(e){
      cat("No Stacking 1!")
    })
    all_coef <- cbind(all_coef, coef.stacking1)
  }
  if(2 %in% meta_index){
    ## Constrained Ridge: sum to one + no intercept 
    coef.stacking2 <- rep(NA, 1 + ncol(Z.mat))
    tryCatch({
      M <- ncol(Z.mat)
      model <- ConstrRidge_sum1(Z.mat, y.vec)
      coef.stacking2 <- c(0, model)
      #model <- ConsLasso_eq(Z.mat, y.vec, C.full = t(rep(1, M)), b = 1, l.min = -3, l.max = 6, step = 0.1, intercept = F)
      #error_MSE <- colMeans((Z.mat %*% model$coefs - y.vec)^2)
      #coef.stacking2 <- c(0, model$coefs[ , which.min(error_MSE)])
    }, error = function(e){
      cat("No Stacking 2!")
    })
    all_coef <- cbind(all_coef, coef.stacking2)
  }
  if(3 %in% meta_index){
    ## Constrained Ridge: sum to one + positive + no intercept 
    coef.stacking3 <- rep(NA, 1 + ncol(Z.mat))
    tryCatch({
      M <- ncol(Z.mat)
      model <- ConstrRidge(Z.mat, y.vec, l.min = -3, l.max = 6, step = 0.1)
      coef.stacking3 <- c(0, model)
    }, error = function(e){
      cat("No Stacking 3!")
    })
    all_coef <- cbind(all_coef, coef.stacking3)
  }
  if(4 %in% meta_index){
    ## Constrained Ridge: sum to one + positive + no intercept + prior weights 
    coef.stacking4 <- rep(NA, 1 + ncol(Z.mat))
    tryCatch({
      M <- ncol(Z.mat)
      model <- ConstrRidge.priorw(X = Z.mat, y = y.vec, l1.seq = 10^seq(-3, 6, 0.1), wprior = c(1, rep(0, M-1)))
      coef.stacking4 <- c(0, model)
    }, error = function(e){
      cat("No Stacking 4!")
    })
    all_coef <- cbind(all_coef, coef.stacking4)
  }
  return(all_coef)
}
###### stacking method #######
SuperLearner <- function(X, y, weights, meta_index, w5.select, w6.select){
  # find NA weights 
  remove.index <- which(is.na(weights[1, ]) == TRUE)
  # Preparing the leave-one-out data
  if(length(remove.index) > 0){
    #Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights) + 1 - length(remove.index)) # change the number of models
    Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights) - length(remove.index)) # change the number of models
    y.vec <- c()
    num_eachfold <- floor(nrow(X)/10)
    start_n <- seq(1, nrow(X), by = num_eachfold)
    end_n <- c(start_n[2:length(start_n)] - 1, nrow(X))
    for(n in 1:length(start_n)){
      train.X <- X[-c(start_n[n]:end_n[n]), ]
      test.X <- X[start_n[n]:end_n[n], ]
      train.y <- y[-c(start_n[n]:end_n[n])]
      test.y <- y[start_n[n]:end_n[n]]
      X.cv <- as.matrix(train.X)
      y.cv <- as.vector(train.y)
      # train all models 
      Allweights.cv <- WeightsCalculation(X.cv, y.cv, weight_index = c(1:12)[-remove.index], w5.select, w6.select)
      weights.cv <- Allweights.cv$weights
      pred1 <- as.matrix(test.X) %*% weights.cv
      pred2 <- as.matrix(test.X - Allweights.cv$bias) %*% weights.cv
      #pred3 <- apply(test.X , 1, function(x) median(na.omit(x)))
      #Z.mat[start_n[n]:end_n[n], ] <- cbind(pred1[, 1:2], pred2[, 3:5], pred1[, 6:ncol(weights.cv)], pred3)
      Z.mat[start_n[n]:end_n[n], ] <- cbind(pred1[, 1:2], pred2[, 3:5], pred1[, 6:ncol(weights.cv)])
      y.vec[start_n[n]:end_n[n]] <- test.y
    }
  } else {
    #Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights) + 1) # change the number of models
    Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights)) # change the number of models
    y.vec <- c()
    num_eachfold <- floor(nrow(X)/10)
    start_n <- seq(1, nrow(X), by = num_eachfold)
    end_n <- c(start_n[2:length(start_n)] - 1, nrow(X))
    for(n in 1:length(start_n)){
      train.X <- X[-c(start_n[n]:end_n[n]), ]
      test.X <- X[start_n[n]:end_n[n], ]
      train.y <- y[-c(start_n[n]:end_n[n])]
      test.y <- y[start_n[n]:end_n[n]]
      X.cv <- as.matrix(train.X)
      y.cv <- as.vector(train.y)
      # train all models 
      Allweights.cv <- WeightsCalculation(X.cv, y.cv, weight_index = c(1:12), w5.select, w6.select)
      weights.cv <- Allweights.cv$weights
      if(start_n[n] != end_n[n]){
        pred1 <- as.matrix(test.X) %*% weights.cv
        pred2 <- as.matrix(test.X - Allweights.cv$bias) %*% weights.cv
        #pred3 <- apply(test.X , 1, function(x) median(na.omit(x)))
        #Z.mat[start_n[n]:end_n[n], ] <- cbind(pred1[, 1:2], pred2[, 3:5], pred1[, 6:ncol(weights.cv)], pred3)
        Z.mat[start_n[n]:end_n[n], ] <- cbind(pred1[, 1:2], pred2[, 3:5], pred1[, 6:ncol(weights.cv)])
        y.vec[start_n[n]:end_n[n]] <- test.y
      } else {
        pred1 <- t(as.matrix(test.X)) %*% weights.cv
        pred2 <- t(as.matrix(test.X - Allweights.cv$bias)) %*% weights.cv
        #pred3 <- median(test.X)
        #Z.mat[start_n[n], ] <- c(pred1[1:2], pred2[3:5], pred1[6:ncol(weights.cv)], pred3)
        Z.mat[start_n[n], ] <- c(pred1[1:2], pred2[3:5], pred1[6:ncol(weights.cv)])
        y.vec[start_n[n]] <- test.y
      }
    }
    na.index <- which(is.na(Z.mat[ , 6]))
    if(length(na.index) > 0){
      Z.mat <- Z.mat[-na.index, ]
      y.vec <- y.vec[-na.index]
    }
  }
  # stacking algorithm 
  coef.stacking <- MetaAlgorithm(Z.mat, y.vec, meta_index = meta_index)
  return(coef.stacking)
}
SuperLearner_onebyone <- function(X, y, weights, meta_index, w5.select, w6.select){
  # find NA weights 
  remove.index <- which(is.na(weights[1, ]) == TRUE)
  # Preparing the leave-one-out data
  if(length(remove.index) > 0){
    #Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights) + 1 - length(remove.index)) # change the number of models
    Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights) - length(remove.index)) # change the number of models
    y.vec <- c()
    #num_eachfold <- floor(nrow(X)/10)
    #start_n <- seq(1, nrow(X), by = num_eachfold)
    #end_n <- c(start_n[2:length(start_n)] - 1, nrow(X))
    for(n in 1:nrow(X)){
      train.X <- X[-n, ]
      test.X <- X[n, ]
      train.y <- y[-n]
      test.y <- y[n]
      X.cv <- as.matrix(train.X)
      y.cv <- as.vector(train.y)
      # train all models 
      Allweights.cv <- WeightsCalculation(X.cv, y.cv, weight_index = c(1:12)[-remove.index], w5.select, w6.select)
      weights.cv <- Allweights.cv$weights
      pred1 <- t(as.matrix(test.X)) %*% weights.cv
      pred2 <- t(as.matrix(test.X - Allweights.cv$bias)) %*% weights.cv
      #pred3 <- apply(test.X , 1, function(x) median(na.omit(x)))
      #Z.mat[start_n[n]:end_n[n], ] <- cbind(pred1[, 1:2], pred2[, 3:5], pred1[, 6:ncol(weights.cv)], pred3)
      Z.mat[n, ] <- cbind(pred1[, 1:2], pred2[, 3:5], pred1[, 6:ncol(weights.cv)])
      y.vec[n] <- test.y
    }
  } else {
    #Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights) + 1) # change the number of models
    Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights)) # change the number of models
    y.vec <- c()
    #num_eachfold <- floor(nrow(X)/10)
    #start_n <- seq(1, nrow(X), by = num_eachfold)
    #end_n <- c(start_n[2:length(start_n)] - 1, nrow(X))
    for(n in 1:nrow(X)){
      train.X <- X[-n, ]
      test.X <- X[n, ]
      train.y <- y[-n]
      test.y <- y[n]
      X.cv <- as.matrix(train.X)
      y.cv <- as.vector(train.y)
      # train all models 
      Allweights.cv <- WeightsCalculation(X.cv, y.cv, weight_index = c(1:12), w5.select, w6.select)
      weights.cv <- Allweights.cv$weights
      pred1 <- t(as.matrix(test.X)) %*% weights.cv
      pred2 <- t(as.matrix(test.X - Allweights.cv$bias)) %*% weights.cv
      #pred3 <- median(test.X)
      #Z.mat[start_n[n], ] <- c(pred1[1:2], pred2[3:5], pred1[6:ncol(weights.cv)], pred3)
      Z.mat[n, ] <- c(pred1[1:2], pred2[3:5], pred1[6:ncol(weights.cv)])
      y.vec[n] <- test.y
    }
    na.index <- which(is.na(Z.mat[ , 6]))
    if(length(na.index) > 0){
      Z.mat <- Z.mat[-na.index, ]
      y.vec <- y.vec[-na.index]
    }
  }
  # stacking algorithm 
  coef.stacking <- MetaAlgorithm(Z.mat, y.vec, meta_index = meta_index)
  return(coef.stacking)
}
