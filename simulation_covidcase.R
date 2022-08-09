setwd()
#### load the libraries and functions ####
source("WeightingModels.R")
source("StackingAlgorithm.R")
#### simulation procedures: bootstrapping the judgment data ####
## Prepare the data 
data_case <- read.csv("all_data_case_county.csv", header = T, stringsAsFactors = F)
data_case <- data_case[ , -1] #1:7, 66:70 are indicators
targets <- unique(data_case$target)
fips <- unique(data_case$fips)
na_threshold <- 0.1
n.vector <- 2^seq(4, 8, by = 1)
n.test <- 1000
for(i in 1:4){ # the first four targets
  for(j in 1:length(fips)){
    dat <- data_case[which(data_case$target == targets[i] & data_case$fips == fips[j]), ]
    na_ratio <- apply(dat, 2, function(x) length(which(is.na(x)))/nrow(dat))
    dat <- dat[ , -which(na_ratio >= na_threshold)]
    dat <- na.omit(dat)
    print(paste("#Observation: ", nrow(dat), sep = ""))
    print(paste("#Forecasters: ", ncol(dat)-11, sep = ""))
    if(nrow(dat) >= 32 & ncol(dat) >= 14){
      M <- ncol(dat)-11
      X <- as.matrix(dat[ , c(8:(ncol(dat)-4))]) # this is for inflation rate data
      Y <- dat$inc_case.y
      ## Estimate the bias and covariance matrix
      XY_all <- cbind(X, Y)
      mu_xy.emp <- colMeans(XY_all)
      Sigma_xy.emp <- cov(XY_all)
      for(r in 1:length(n.vector)){
        tryCatch({
          n <- n.vector[r]
          set.seed(i*j*r*361)
          # generate a set of judgments with sample size n 
          XY_all.train <- rmvn(n, mu_xy.emp, Sigma_xy.emp)
          XY_all.test  <- rmvn(n.test, mu_xy.emp, Sigma_xy.emp)
          X.train <- XY_all.train[ , 1:M]
          Y.train <- XY_all.train[ , (M+1)]
          X.test  <- XY_all.test[ , 1:M]
          Y.test  <- XY_all.test[ , (M+1)]
          
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
          print(paste(i, "-", j, "-", r, sep = ""))
          # recording data line by line
          write.table(t(c(i, j, n, cor_vec)), "sim_covidcase_corr.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
          write.table(t(c(i, j, n, sqrt(colMeans((all.pred - Y.test)^2)))), "sim_covidcase_rmse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
          write.table(t(c(i, j, n, coef.stacking)), "sim_covidcase_coef.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
          
        }, error = function(e){ cat("One Error Found!") })
      }
    }
  }
}
