setwd("")
#### load the libraries and functions ####
source("WeightingModels.R")
source("StackingAlgorithm.R")
#### simulation procedures: bootstrapping the judgment data ####
infdata <- read.csv("inf_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 9
X <- as.matrix(infdata[ , 4:12]) 
Y <- infdata[ , 3]
## estimate the bias and covariance matrix -> simulation seed
muy.emp <- mean(Y)
sdy.emp <- sd(Y)
mu.emp <- colMeans(X - Y)
sigma.emp <- cov(X)
## parameter settings for simulation 
n.vec.train <- 2^seq(4, 8, by = 1)
n.test <- 1000
draws <- 100
#rmse.test  <- array(NA, dim = c(length(n.vec.train), draws, 12+1)) # 12 component weighting methods + our stacking method
for(j in 1:length(n.vec.train)){
  n <- n.vec.train[j]
  for(k in 1:draws){
    set.seed(j*k*3927)
    
    # generate the true state
    y.true <- rnorm(100, muy.emp, sdy.emp)
    x.true <- rmvn(100, mu.emp+muy.emp, sigma.emp)
    muy.true <- mean(y.true)
    sdy.true <- sd(y.true)
    mu.true <- colMeans(x.true - y.true)
    sigma.true <- cov(x.true)
    
    # generate a set of judgments 
    Y.train <- rnorm(n, muy.true, sdy.true)
    X.train <- rmvn(n, mu.true+muy.true, sigma.true)
    Y.test <- rnorm(n.test, muy.true, sdy.true)
    X.test <- rmvn(n.test, mu.true+muy.true, sigma.true)
    
    # out-of-sample prediction 
    #Allweights <- WeightsCalculation(X.train, Y.train, weight_index = 1:11, w5.select = 5:8, w6.select = 9:10) # 9 is ridge regression; 10 is lasso regression 
    Allweights <- WeightsCalculation(X.train, Y.train, weight_index = 1:12, w5.select = 5:8, w6.select = 9) # 9 is ridge regression; 10 is lasso regression 
    # w3, w4, w5 need to use bias; w5 needs to be selected 
    ## select the best w5
    debias <- Allweights$bias
    debias.mat <- matrix(rep(debias, each = nrow(X.train)), ncol = ncol(X.train), nrow = nrow(X.train))
    w5.error <- colMeans(((X.train - debias.mat) %*% Allweights$weights[ , 5:8] - Y.train)^2)
    w5.select <- c(5:8)[which.min(w5.error)]
    w5.optw_1cor <- Allweights$weights[, w5.select]
    weights <- cbind(Allweights$weights[ , 1:4], w5.optw_1cor, Allweights$weights[ , 9:15])
    pred.ind <- MakePrediction(X.test, weights, debias)
    # correlation among different models 
    cormat <- cor(pred.ind)
    cor_vec <- cormat[upper.tri(cormat)] # upper half by column 
    # our stacking method 
    coef.stacking <- SuperLearner(X.train, Y.train, weights, meta_index = 4, w5.select, w6.select = 9)
    select.wa <- which(is.na(pred.ind[1, ]) == FALSE)
    pred.stacking <- cbind(rep(1, nrow(pred.ind)), pred.ind[ , select.wa]) %*% coef.stacking
    # combine all predictions 
    all.pred <- cbind(pred.ind, pred.stacking)
    
    # record all prediction errors 
    #rmse.test[j, k, ]  <- sqrt(colMeans((all.pred - Y.test)^2))
    print(paste(j, "-", k, sep = ""))
    # recording data line by line
    write.table(t(c(n.vec.train[j], k, cor_vec)), "sim_inf_corr.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(t(c(n.vec.train[j], k, sqrt(colMeans((all.pred - Y.test)^2)))), "sim_inf_rmse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    write.table(t(c(n.vec.train[j], k, coef.stacking)), "sim_inf_coef.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
}
