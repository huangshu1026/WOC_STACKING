setwd()
#### load the libraries and functions ####
source("WeightingModels.R")
source("StackingAlgorithm.R")
#### Keck and Tang (2020)'s experimental data ####
library(readxl)
exp_dat <- read_excel("DataStudy1.xlsx")

# Start the anlaysis 
M <- 5 # 48 forecasters with 40 forecasts 
n.train.vec <- c(10, 15, 20, 25, 30)
draws <- 100
condition <- unique(exp_dat$con) # 3 conditions
for(i in 1:3){ # index of condition
  for(j in 1:length(n.train.vec)){
    n <- n.train.vec[j]
    for(k in 1:draws){
      tryCatch({
        set.seed(i*j*k*27)
        subdata <- exp_dat[which(exp_dat$con == condition[i]), ] 
        subdata <- na.omit(subdata)
        subdata <- subdata[sample(1:nrow(subdata), M, replace = FALSE), ]
        X.all <- t(subdata[ , 2:41])
        Y.all <- t(subdata[1, 42:81])
        index_train <- sample(1:nrow(X.all), n, replace = FALSE)
        index_test <- c(1:nrow(X.all))[-index_train]
        X.train <- as.matrix(X.all[index_train, ])
        Y.train <- as.vector(Y.all[index_train])
        X.test <- as.matrix(X.all[index_test, ])
        Y.test <- as.vector(Y.all[index_test])
        
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
        #cor_vec_record <- rbind(cor_vec_record, cor_vec)
        # our stacking method 
        coef.stacking <- SuperLearner(X.train, Y.train, weights, meta_index = 4, w5.select, w6.select = 9)
        #coef.stacking_record <- rbind(coef.stacking_record, coef.stacking)
        select.wa <- which(is.na(pred.ind[1, ]) == FALSE)
        pred.stacking <- cbind(rep(1, nrow(pred.ind)), pred.ind[ , select.wa]) %*% coef.stacking
        # combine all predictions 
        all.pred <- cbind(pred.ind, pred.stacking)
        
        # record all prediction errors 
        rmse.test  <- sqrt(colMeans((all.pred - Y.test)^2))
        print(paste(i, "-", j, "-", k, sep = ""))
        write.table(t(c(i, n.train.vec[j], k, rmse.test)), "emp_keck_rmse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
        write.table(t(c(i, n.train.vec[j], k, cor_vec)), "emp_keck_corr.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
        write.table(t(c(i, n.train.vec[j], k, coef.stacking)), "emp_keck_coef.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
      }, error = function(e){cat("Error Report!")})
    }
  }
}

