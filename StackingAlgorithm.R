##############################################################################################
########## This file includes function to get all weights and the stacking method  ###########
##############################################################################################
#### function to calculate different weights ####
WeightsCalculation <- function(X, y){
  # w1: simple average 
  w1.sa <- rep(1/ncol(X), ncol(X))
  # w2: optimal weights 
  w2.optw <- optw(X, y)
  # w3: CEWM
  w3.cwm <- CWM(X, y)$weight.CWM 
  # w4: sequential search 
  w4.ssin <- SeqSearch(X, y)$weight.SSIN 
  w4.ssde <- SeqSearch(X, y)$weight.SSDE
  # w5: ranked performance 
  w5.rp <- RankPerformance(X, y)
  # w6: top 3 and 5 models
  w6.top5 <- TopNModel(X, y, 5)
  # w7: AppoxGopt
  w7.ApproxGopt <- ApproxGopt(X, y)
  # w8: constrained lasso regression with zero prior 
  w8.ConLas <- rep(NA, ncol(X))
  tryCatch({
    w8.ConLas <- ConstrENetZERO(X, y, addl1 = TRUE, addl2 = FALSE, l1.seq = 10^(seq(-3, 4, by = 0.1)), l2.seq = NULL)$weights
  }, error=function(e){
    cat("ERROR: No solution for ConLas(zero)!")
  })
  # w9: constrained lasso regression with equal-weight prior 
  w9.OurConLas <- rep(NA, ncol(X))
  tryCatch({
    w9.OurConLas <- ConstrENetEQAUL(X, y, addl1 = TRUE, addl2 = FALSE, l1.seq = 10^(seq(-3, 4, by = 0.1)), l2.seq = NULL, wprior = rep(1/ncol(X), ncol(X)))$weights
  }, error=function(e){
    cat("ERROR: No solution for OurConLas!")
  })
  # record all weights for testing data set 
  weights <- cbind(w1.sa, w2.optw, w3.cwm, 
                   w4.ssin, w4.ssde, w5.rp, w6.top5, w7.ApproxGopt, 
                   w8.ConLas, w9.OurConLas)
  return(weights)
}

###### stacking method #######
SuperLearner <- function(X, y, weights){
  # Preparing the leave-one-out data
  Z.mat <- matrix(NA, nrow = nrow(X), ncol = ncol(weights)) # change the number of models
  y.vec <- c()
  for(n in 1:nrow(X)){
    train.X <- X[-n, ]
    test.X <- X[n, ]
    train.y <- y[-n]
    test.y <- y[n]
    X.cv <- as.matrix(train.X)
    y.cv <- as.vector(train.y)
    # train all models 
    weights.cv <- WeightsCalculation(X.cv, y.cv)
    Z.mat[n, ] <- t(weights.cv) %*% as.matrix(test.X)
    y.vec[n] <- test.y
  }
  # Identify weighting models without solution
  model.index <- which(is.na(Z.mat[1, ]) == FALSE)
  Z.mat <- Z.mat[ , model.index]
  # stacking algorithm 
  coef.stacking <- rep(NA, 1 + ncol(Z.mat))
  tryCatch({
    lambda_seq <- 10^seq(-3, 4, by = 0.1)
    stacking <- cv.glmnet(Z.mat, y.vec, alpha = 1, lambda = lambda_seq)
    best_lambda <- stacking$lambda.min
    stacking <- glmnet(Z.mat, y.vec, alpha = 1, lambda = best_lambda)
    coef.stacking <- as.vector(coef(stacking))
  }, error = function(e){
    cat("No solution for Stacking Algorithm.")
  })
  return(list(coef.stacking = coef.stacking, model.index = model.index))
}
