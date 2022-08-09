##################################################################################################
########## This file includes all previous individual WOC models (i.e., base leaners) ############
##################################################################################################

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

##########################################
########## Simple Approaches #############
##########################################
###### Simple Average (SA) ######
eqw <- function(X, y){
  M <- ncol(X)
  weight <- rep(1/M, M)
  return(weight)
}

###################################################
######### All-crowd weighting methods #############
###################################################
###### Optimal weighting ######
## Davis-Stover (2014, 2015)
optw_bias <- function(X, y){
  mu.obs <- colMeans(X) - mean(y)
  sigma.obs <- cov(X)
  mat <- sigma.obs + mu.obs%*%t(mu.obs)
  sigma.xy <- cov(X, y)
  A <- t(rep(1, ncol(X)))
  QPmodel <- solve.QP(Dmat = nearPD(mat)$mat, 
                      dvec = sigma.xy, 
                      Amat = t(A), 
                      bvec = c(1), 
                      meq = 1) 
  weight <- QPmodel$solution
  return(weight)
}
## Winkler (1986): require debias
optw_cov <- function(X, y){
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  sigma <- cov(X.new)
  sigma <- nearPD(sigma)$mat
  weight <- colSums(solve(sigma))/sum(solve(sigma))
  return(list(weight = weight, bias = mu.obs))
}
## Bunn (1985): require debias
optw_var <- function(X, y){
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  sigma <- diag(apply(X.new, 2, var))
  sigma <- nearPD(sigma)$mat
  weight <- colSums(solve(sigma))/sum(solve(sigma))
  return(list(weight = weight, bias = mu.obs))
}
## Soule et al., (2020): require debias
Calweight_rho_c <- function(rho_c, sd_vec){
  k <- length(sd_vec)
  cor.mat <- matrix(rep(rho_c, k*k), k, k)
  diag(cor.mat) <- 1
  sigma <- cor2cov(cor.mat, sd_vec)
  sigma <- nearPD(sigma)$mat 
  weight <- colSums(solve(sigma))/sum(solve(sigma))
  return(weight)
}
optw_var_1cor_exo <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model
  k <- ncol(X.new)
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  rho_c_seq <- seq(1/(1-k)+0.001, 0.999, length.out = 1000)
  Weights <- c()
  for(i in 1:length(rho_c_seq)){
    tempdata <- Calweight_rho_c(rho_c_seq[i], sd_vec)
    Weights <- rbind(Weights, tempdata)
  }
  error_MSE <- colMeans((X.new %*% t(Weights) - y)^2)
  best_rho_c <- rho_c_seq[which.min(error_MSE)]
  weight <- Calweight_rho_c(best_rho_c, sd_vec)
  return(list(rho_c = best_rho_c, weight = weight, bias = mu.obs))
}
optw_var_1cor_BayesMean <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model
  k <- ncol(X.new)
  delta0 <- k + 1
  rho_c_0 <- 0
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  cor_mat <- cor(X.new)
  if(length(which(is.na(cor_mat) == TRUE)) > 0){
    k.na <- length(which(is.na(cor_mat) == TRUE))
    rho_bar <- (sum(na.omit(cor_mat)) - k)/(k*(k-1) - k.na)
  } else {
    rho_bar <- (sum(cor_mat)-k)/(k*(k-1))
  }
  rho_c <- ((delta0-k-1)*rho_c_0 + n*rho_bar) / (delta0+n-k-1)
  weight <- Calweight_rho_c(rho_c, sd_vec)
  return(list(rho_c = rho_c, weight = weight, bias = mu.obs))
}
optw_var_1cor_BayesMAP <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model 
  k <- ncol(X.new)
  delta0 <- k + 1
  rho_c_0 <- 0
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  cor_mat <- cor(X.new)
  if(length(which(is.na(cor_mat) == TRUE)) > 0){
    k.na <- length(which(is.na(cor_mat) == TRUE))
    rho_bar <- (sum(na.omit(cor_mat)) - k)/(k*(k-1) - k.na)
  } else {
    rho_bar <- (sum(cor_mat)-k)/(k*(k-1))
  }
  rho_c <- ((delta0-k-1)*rho_c_0 + n*rho_bar) / (delta0+n+k+2)
  weight <- Calweight_rho_c(rho_c, sd_vec)
  return(list(rho_c = rho_c, weight = weight, bias = mu.obs))
}
optw_var_1cor_mean <- function(X, y){
  # obtain unbiased judgments 
  n <- nrow(X)
  mu.obs <- colMeans(X - y)
  debias.mat <- matrix(rep(mu.obs, each = n), ncol = ncol(X), nrow = nrow(X))
  X.new <- X - debias.mat
  # model 
  k <- ncol(X.new)
  sd_vec <- apply(X.new, 2, sd)
  if(length(which(sd_vec <= 0)) > 0){
    sd_vec[which(sd_vec <= 0)] <- min(sd_vec[which(sd_vec > 0)])*0.001
  }
  cor_mat <- cor(X.new)
  if(length(which(is.na(cor_mat) == TRUE)) > 0){
    k.na <- length(which(is.na(cor_mat) == TRUE))
    rho_c <- (sum(na.omit(cor_mat)) - k)/(k*(k-1) - k.na)
  } else {
    rho_c <- (sum(cor_mat)-k)/(k*(k-1))
  }
  weight <- Calweight_rho_c(rho_c, sd_vec)
  return(list(rho_c = rho_c, weight = weight, bias = mu.obs))
}

###### Regularized weighting ######
ConsLas <- function(X, y){
  M <- ncol(X)
  w.prime <- rep(1/M, M)
  y.new <- y - X %*% w.prime
  model <- ConsLasso_eq(X, y.new, C.full = t(rep(1, M)), b = 0, l.min = -3, l.max = 6, step = 0.1, intercept = F)
  # select the best lambda 
  error_MSE <- colMeans((X %*% (model$coefs + w.prime) - y)^2)
  #best.labmda <- model$lambda[which.min(error_MSE)] 
  weight <- model$coefs[ , which.min(error_MSE)] + w.prime
  return(weight)
}
ConsRidge <- function(X, y){
  M <- ncol(X)
  w.prime <- rep(1/M, M)
  y.new <- y - X %*% w.prime
  w <- ConstrRidge_sum0(X, y.new)
  weight <- w + w.prime
  return(weight)
}
ConstrRidge_sum0 <- function(X, y, l.min = -3, l.max = 6, step = 0.1){
  l2.seq <- 10^seq(l.min, l.max, by = step)
  bestl2 <- NA
  M <- ncol(X)
  n <- nrow(X)
  XTX <- t(X) %*% X
  XTY <- t(X) %*% y
  
  gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ridge_sum0(X, y, l2.seq[i]))
  if(all(is.na(gcv) == TRUE)){
    l2.seq.new <- 10^seq(l.min-3, l.min, by = 0.1)
    gcv <- sapply(1:length(l2.seq.new), function(i) CalGCV.ridge_sum0(X, y, l2.seq.new[i]))
    if(all(is.na(gcv) == TRUE)){
      bestl2 <- 0
    } else {
      bestl2 <- l2.seq.new[which.min(gcv)]
    }
  } else {
    bestl2 <- l2.seq[which.min(gcv)]
  }
  H <- cbind(XTX + bestl2*diag(M))
  f <- XTY
  # equalities
  A.eq <- rbind(rep(1, M))
  b.eq <- c(0)
  results <- solve.QP(Dmat = H, dvec = f, Amat = t(rbind(A.eq)), bvec = c(b.eq), meq = 1)
  weights <- results$solution
  return(weights)
}
CalGCV.ridge_sum0 <- function(X, y, lambda2){
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
ConsLasso_ineq <- function (x, y, C.full, b, l.min = -2, l.max = 6, step = 0.2, 
                            beta0 = NULL, verbose = F, max.it = 12, intercept = T, normalize = T, 
                            backwards = F){
  if (!backwards) {
    fit = ConsLars_ineq(x, y, C.full, b, l.min = l.min, l.max = l.max, 
                        step = step, beta0 = beta0, verbose = verbose, max.it = max.it, 
                        intercept = intercept, normalize = normalize, forwards = T)
    if (fit$error | is.null(fit$coefs)) {
      if (is.null(fit$lambda)) 
        fit$lambda = 10^l.min
      fit2 = ConsLars_ineq(x, y, C.full, b, l.min = log10(max(fit$lambda)), 
                           l.max = l.max, step = step, beta0 = beta0, verbose = verbose, 
                           max.it = max.it, intercept = intercept, normalize = normalize, 
                           forwards = F)
      if (is.null(fit$coefs)) 
        fit$lambda = NULL
      fit$coefs = cbind(fit$coefs, fit2$coefs)
      fit$lambda = c(fit$lambda, fit2$lambda)
      fit$intercept = c(fit$intercept, fit$intercept)
    }
  }
  else fit = ConsLars_ineq(x, y, C.full, b, l.min = l.min, l.max = l.max, 
                           step = step, beta0 = beta0, verbose = verbose, max.it = max.it, 
                           intercept = intercept, normalize = normalize, forwards = F)
  fit
}
ConsLars_ineq <- function (x, y, C.full, b, l.min = -2, l.max = 6, step = 0.2, 
                           beta0 = NULL, verbose = F, max.it = 12, intercept = T, normalize = T, 
                           forwards = T){
  p = ncol(x)
  n = nrow(x)
  M = nrow(C.full)
  one <- rep(1, n)
  if (intercept) {
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- drop(y - mu)
  }
  else {
    meanx <- rep(0, p)
    mu <- 0
    y <- drop(y)
  }
  normx <- rep(1, p)
  if (normalize) {
    normx <- sqrt(n) * apply(x, 2, stats::sd, na.rm = T)
    x <- scale(x, FALSE, normx)
  }
  C.full = t(t(C.full)/normx)
  if (!forwards) {
    lambda = lambda.old = 10^l.max
    step = -step
    if (is.null(beta0)) 
      beta0 = lin.int.ineq(C.full, b)
  }
  else {
    lambda = lambda.old = 10^l.min
    if (is.null(beta0)) 
      beta0 = quad.int.ineq(x, y, C.full, b, lambda)
  }
  A.old = C.full
  x.old = x
  C.full = cbind(C.full, -diag(M))
  x = cbind(x, matrix(0, n, M))
  p = p + M
  beta.new = rep(0, p)
  step.orig = step
  coefs = grid = b2index = NULL
  t.data = transformed.ineq(x, y, C.full, b, lambda, beta0)
  beta1.old = t.data$beta1
  beta2.old = t.data$beta2
  iterations = 1
  end.loop = F
  while (!end.loop & (iterations <= max.it)) {
    iterations = 1
    loop = T
    while (loop & (iterations <= max.it)) {
      t.data$y = t.data$y + (lambda - lambda.old) * t.data$C
      beta1.new = rep(0, length(beta1.old))
      act = length(t.data$active)
      lambda.pen = rep(lambda, act)
      lambda.pen[t.data$delta1.index] = 0
      positive.pen = rep(F, act)
      if (act > p) {
        positive.pen[t.data$delta1.index] = T
      }
      fit.pen = Penalized_update(t.data$y, t.data$x[ ,
                                                     t.data$active], ~0, lambda1 = lambda.pen, lambda2 = 0.001, trace = F, positive = positive.pen, 
                                 maxiter = 10000, startbeta = beta1.old, epsilon = 10^-6)
      beta1.new[t.data$active] = coef(fit.pen, "all")
      # if (!attr(fit.pen, "converged")) 
      #   print(paste("Warning: Maximum iterations exceeded in penalized fit."))
      beta2.new = beta2.old + t.data$a2 %*% (beta1.old - 
                                               beta1.new)
      bad.beta2 = (sum(abs(sign(t.data$beta2) - sign(beta2.new))) != 
                     0)
      X_star = t.data$x
      derivs = abs(as.vector(t(X_star) %*% (X_star %*% 
                                              beta1.new)) - t(X_star) %*% t.data$Y_star - 
                     lambda * t.data$C2)
      bad.active = F
      if (n < (p - M)) 
        bad.active = (max(derivs[-t.data$active]) > 
                        lambda)
      if (bad.beta2 | bad.active) {
        t.data$y = t.data$y - (lambda - lambda.old) * 
          t.data$C
        step = step/2
        lambda = lambda.old * 10^step
        iterations = iterations + 1
      }
      else loop = F
    }
    if (iterations <= max.it) {
      # if (verbose == T) {
      #   print(paste("Lambda =", round(lambda, 3)))
      #   if (abs(step) < abs(step.orig)) 
      #     print(paste("Step size reduced to ", step))
      # }
      step = step.orig
      beta.new[t.data$beta2.index] = beta2.new
      beta.new[-t.data$beta2.index] = beta1.new
      coefs = cbind(coefs, beta.new)
      change.beta = (min(abs(beta2.new)) < max(abs(beta1.new)))
      change.active = F
      if (n < (p - M)) 
        change.active = (min(derivs[t.data$active]) < 
                           max(derivs[-t.data$active]))
      if (change.beta | change.active) {
        t.data = transformed.ineq(x, y, C.full, b, lambda, 
                                  beta.new)
        beta1.new = t.data$beta1
        beta2.new = t.data$beta2
      }
      beta1.old = beta1.new
      beta2.old = beta2.new
      grid = c(grid, lambda)
      lambda.old = lambda
      lambda = lambda * 10^step
    }
    if ((forwards & (lambda > 10^l.max)) | (!forwards & 
                                            (lambda < 10^l.min))) 
      end.loop = T
  }
  # if (iterations > max.it) 
  # print("Exceed.")
  # print(paste("Warning: Algorithm terminated at lambda =", 
  #             round(lambda.old, 1), ": Maximum iterations exceeded."))
  colnames(coefs) = intercept = NULL
  if (!is.null(grid)) {
    normx = c(normx, rep(1, M))
    coefs = coefs/normx
    coefs = coefs[, order(grid)]
    grid = sort(grid)
    meanx = c(meanx, rep(0, M))
    intercept = mu - drop(t(coefs) %*% meanx)
  }
  coefs = coefs[1:(p - M), ]
  list(coefs = coefs, lambda = grid, intercept = intercept, 
       error = (iterations > max.it))
}
Penalized_update <- function (response, penalized, unpenalized, lambda1 = 0, lambda2 = 0, 
                              positive = FALSE, data, fusedl = FALSE, model = c("cox", "logistic", "linear", "poisson"), startbeta, startgamma, 
                              steps = 1, epsilon = 1e-10, maxiter, standardize = FALSE, 
                              trace = TRUE){
  if (missing(maxiter)) 
    maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive) 
      25
  else Inf
  if (steps == "Park" || steps == "park") {
    steps <- 1
    park <- TRUE
  }
  else park <- FALSE
  prep <- penalized:::.checkinput(match.call(), parent.frame())
  if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) && 
      all(lambda2 == 0) && !any(prep$positive)) 
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", 
         call. = FALSE)
  fit <- penalized:::.modelswitch(prep$model, prep$response, prep$offset, 
                                  prep$strata)$fit
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  nr <- nrow(prep$X)
  fusedl <- prep$fusedl
  if (length(lambda1) == pp && (!all(lambda1 == 0))) {
    wl1 <- c(numeric(pu), lambda1)
    lambda1 <- 1
  }
  else {
    wl1 <- 1
  }
  if (length(lambda2) == pp) 
    lambda2 <- c(numeric(pu), lambda2)
  if (park || steps > 1 && fusedl == FALSE) {
    if (pu > 0) 
      lp <- drop(prep$X[, 1:pu, drop = FALSE] %*% prep$nullgamma)
    else lp <- numeric(n)
    chck <- (wl1 > 0) & c(rep(FALSE, pu), rep(TRUE, pp))
    gradient <- drop(crossprod(prep$X[, chck, drop = FALSE], 
                               fit(lp)$residuals))
    if (length(wl1) > 1) {
      rel <- gradient/(wl1[chck] * prep$baselambda1[chck])
    }
    else {
      rel <- gradient/(wl1 * prep$baselambda1[chck])
    }
    from <- max(ifelse(prep$positive[chck], rel, abs(rel)))
    if (from < lambda1) {
      warning("Chosen lambda1 greater than maximal lambda1: \"steps\" argument ignored")
      steps <- 1
      park <- FALSE
      from <- lambda1
    }
  }
  else {
    from <- lambda1
  }
  lambda1s <- sapply(1:length(lambda1), function(x) seq(from[x], lambda1[x], length.out = steps))
  beta <- prep$beta
  louts <- if (park) 
    4 * pp
  else length(lambda1s)
  outs <- vector("list", louts)
  rellambda1 <- lambda1s[1]
  ready <- FALSE
  i <- 0
  while (!ready) {
    ready <- all(rellambda1 == lambda1)
    i <- i + 1
    if (!fusedl) {
      if (rellambda1 != 0 || any(prep$positive)) {
        if (all(lambda2 == 0)) {
          out <- penalized:::.steplasso(beta = beta, lambda = rellambda1 * 
                                          wl1 * prep$baselambda1, lambda2 = 0, positive = prep$positive, 
                                        X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
                                        maxiter = maxiter)
        }
        else {
          out <- penalized:::.lasso(beta = beta, lambda = rellambda1 * 
                                      wl1 * prep$baselambda1, lambda2 = lambda2 * 
                                      prep$baselambda2, positive = prep$positive, 
                                    X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
                                    maxiter = maxiter)
        }
      }
      else {
        if (pp > n) {
          P <- penalized:::.makeP(prep$X, lambda2 * prep$baselambda2)
          gams <- .solve(crossprod(t(P)), P %*% beta)
          PX <- P %*% t(prep$X)
          Pl <- P * matrix(sqrt(lambda2 * prep$baselambda2), 
                           nrow(P), ncol(P), byrow = TRUE)
          PlP <- crossprod(t(Pl))
          out <- penalized:::.ridge(beta = gams, Lambda = PlP, X = t(PX), 
                                    fit = fit, trace = trace, epsilon = epsilon, 
                                    maxiter = maxiter)
          out$beta <- drop(crossprod(P, out$beta))
        }
        else {
          out <- penalized:::.ridge(beta = beta, Lambda = lambda2 * 
                                      prep$baselambda2, X = prep$X, fit = fit, 
                                    trace = trace, epsilon = epsilon, maxiter = maxiter)
        }
      }
    }
    if (fusedl) {
      out <- penalized:::.flasso(beta = beta, lambda1 = rellambda1 * 
                                   wl1 * prep$baselambda1, lambda2 = lambda2 * 
                                   prep$baselambda2, chr = prep$chr, positive = prep$positive, 
                                 X = prep$X, fit = fit, trace = trace, epsilon = epsilon, 
                                 maxiter = maxiter)
    }
    if (trace) 
      cat("\n")
    beta <- out$beta
    if (!ready) {
      if (!fusedl) {
        if (park) {
          newpark <- penalized:::.park(beta = beta, lambda = rellambda1 * 
                                         wl1 * prep$baselambda1, lambda2 = 0, positive = prep$positive, 
                                       X = prep$X, fit = out$fit)
          rellambda1 <- rellambda1 * (1 - newpark$hh)
          if (rellambda1 < lambda1 || rellambda1 == 
              Inf) {
            rellambda1 <- lambda1
            beta <- out$beta
          }
          else {
            beta <- newpark$beta
          }
          lambda1s <- c(lambda1s, rellambda1)
        }
        else {
          rellambda1 <- lambda1s[i + 1]
          beta <- out$beta
        }
      }
      else {
        rellambda1 <- lambda1s[i + 1]
        beta <- out$beta
      }
    }
    outs[[i]] <- out
  }
  if (length(lambda2) > 1) 
    lambda2 <- lambda2[pu + 1:pp]
  outs <- sapply(1:i, function(nr) {
    thislambda1 <- lambda1s[[nr]] * ifelse(length(wl1) > 
                                             1, wl1[pu + 1:pp], wl1)
    penalized:::.makepenfit(outs[[nr]], pu, fusedl = fusedl, prep$model, 
                            thislambda1, lambda2, prep$orthogonalizer, prep$weights, 
                            prep$formula, rownames(prep$X))
  })
  if (length(outs) == 1) 
    outs <- outs[[1]]
  outs
}
ConsLasso_eq <- function (x, y, C.full, b, l.min = -2, l.max = 6, step = 0.2, 
            beta0 = NULL, verbose = F, max.it = 12, intercept = T, normalize = T, 
            backwards = F){
    if (!backwards) {
      fit = ConsLars_eq(x, y, C.full, b, l.min = l.min, l.max = l.max, 
                        step = step, beta0 = beta0, verbose = verbose, max.it = max.it, 
                        intercept = intercept, normalize = normalize, forwards = T)
      if (fit$error | is.null(fit$coefs)) {
        if (is.null(fit$lambda)) 
          fit$lambda = 10^l.min
        fit2 = ConsLars_eq(x, y, C.full, b, l.min = log10(max(fit$lambda)), 
                           l.max = l.max, step = step, beta0 = beta0, verbose = verbose, 
                           max.it = max.it, intercept = intercept, normalize = normalize, 
                           forwards = F)
        if (is.null(fit$coefs)) 
          fit$lambda = NULL
        fit$coefs = cbind(fit$coefs, fit2$coefs)
        fit$lambda = c(fit$lambda, fit2$lambda)
        fit$intercept = c(fit$intercept, fit2$intercept)
        fit$b2index = cbind(fit$b2index, fit2$b2index)
      }
    }
    else fit = ConsLars_eq(x, y, C.full, b, l.min = l.min, l.max = l.max, 
                           step = step, beta0 = beta0, verbose = verbose, max.it = max.it, 
                           intercept = intercept, normalize = normalize, forwards = F)
    fit
  }
ConsLars_eq <- function (x, y, C.full, b, l.min = -2, l.max = 6, step = 0.2, 
                         beta0 = NULL, verbose = F, max.it = 12, intercept = T, normalize = T, 
                         forwards = T){
  p = ncol(x)
  n = nrow(x)
  m = nrow(C.full)
  beta.new = rep(0, p)
  one <- rep(1, n)
  if (intercept) {
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- drop(y - mu)
  }
  else {
    meanx <- rep(0, p)
    mu <- 0
    y <- drop(y)
  }
  normx <- rep(1, p)
  if (normalize) {
    normx <- sqrt(n) * apply(x, 2, stats::sd, na.rm = T)
    x <- scale(x, FALSE, normx)
  }
  C.full = t(t(C.full)/normx)
  if (!forwards) {
    lambda = lambda.old = 10^l.max
    step = -step
    if (is.null(beta0)) 
      beta0 = lin.int(C.full, b)
  }
  else {
    lambda = lambda.old = 10^l.min
    if (is.null(beta0)) 
      beta0 = quad.int(x, y, C.full, b, lambda)
  }
  step.orig = step
  coefs = grid = b2index = NULL
  t.data = transformed(x, y, C.full, b, lambda, beta0)
  beta1.old = t.data$beta1
  beta2.old = t.data$beta2
  iterations = 1
  end.loop = F
  while (!end.loop & (iterations <= max.it)) {
    iterations = 1
    loop = T
    while (loop & (iterations <= max.it)) {
      t.data$y = t.data$y + (lambda - lambda.old) * t.data$C
      beta1.new = rep(0, length(beta1.old))
      fit = lars::lars(t.data$x[, t.data$active], t.data$y, 
                       normalize = F, intercept = F)
      beta1.new[t.data$active] = stats::predict(fit, s = lambda, 
                                                type = "coefficients", mode = "lambda")$coef
      beta2.new = beta2.old + t.data$a2 %*% (beta1.old - 
                                               beta1.new)
      bad.beta2 = (sum(abs(sign(t.data$beta2) - sign(beta2.new))) != 
                     0)
      X_star = t.data$x
      derivs = abs(as.vector(t(X_star) %*% (X_star %*% 
                                              beta1.new)) - t(X_star) %*% t.data$Y_star - 
                     lambda * t.data$C2)
      bad.active = F
      if (n < (p - m)) 
        bad.active = (max(derivs[-t.data$active]) > 
                        lambda)
      if (bad.beta2 | bad.active) {
        t.data$y = t.data$y - (lambda - lambda.old) * 
          t.data$C
        step = step/2
        lambda = lambda.old * 10^step
        iterations = iterations + 1
      }
      else loop = F
    }
    if (iterations <= max.it) {
      # if (verbose == T) {
      #   print(paste("Lambda =", round(lambda, 3)))
      #   if (abs(step) < abs(step.orig)) 
      #     print(paste("Step size reduced to ", step))
      # }
      step = step.orig
      beta.new[t.data$beta2.index] = beta2.new
      beta.new[-t.data$beta2.index] = beta1.new
      coefs = cbind(coefs, beta.new)
      b2index = cbind(b2index, t.data$beta2.index)
      change.beta = (min(abs(beta2.new)) < max(abs(beta1.new)))
      change.active = F
      if (n < (p - m)) 
        change.active = (min(derivs[t.data$active]) < 
                           max(derivs[-t.data$active]))
      if (change.beta | change.active) {
        t.data = transformed(x, y, C.full, b, lambda, 
                             beta.new)
        beta1.new = t.data$beta1
        beta2.new = t.data$beta2
      }
      beta1.old = beta1.new
      beta2.old = beta2.new
      grid = c(grid, lambda)
      lambda.old = lambda
      lambda = lambda * 10^step
    }
    if ((forwards & (lambda > 10^l.max)) | (!forwards & 
                                            (lambda < 10^l.min))) 
      end.loop = T
  }
  # if (iterations > max.it) 
  #   print(paste("Warning: Algorithm terminated at lambda =", 
  #               round(lambda.old, 1), ": Maximum iterations exceeded."))
  colnames(coefs) = intercept = NULL
  if (!is.null(grid)) {
    coefs = coefs/normx
    coefs = coefs[, order(grid)]
    b2index = b2index[, order(grid)]
    grid = sort(grid)
    intercept = mu - drop(t(coefs) %*% meanx)
  }
  list(coefs = coefs, lambda = grid, intercept = intercept, 
       error = (iterations > max.it), b2index = b2index)
}

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

#################################################
########## Vanilla Stacking Method  #############
#################################################
VanillaStacking <- function(X, y){
  ## Constrained Ridge: sum to one + positive + no intercept + prior weights 
  Z.mat <- as.matrix(X)
  y.vec <- y
  M <- ncol(Z.mat)
  weights <- rep(NA, M)
  tryCatch({
    model <- ConstrRidge.priorw(X = Z.mat, y = y.vec, l1.seq = 10^seq(-3, 6, 0.1), wprior = rep(1/M, M))
    weights <- model
  }, error = function(e){
    cat("No VanillaStacking!")
  })
  return(weights)
}





