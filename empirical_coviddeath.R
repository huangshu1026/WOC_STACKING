####################################################################
#### Solve the empirical results for vanilla stacking algorithm ####
####################################################################
setwd()
#### load the libraries and functions ####
source("WeightingModels.R")
source("StackingAlgorithm.R")

##### COVID-19 death #####
data_death <- read.csv("all_data_death.csv", header = T, stringsAsFactors = F)
data_death <- data_death[ , -1] #1:4, 105:107 are indicators
# grouped data 
targets <- unique(data_death$target)
fips <- unique(data_death$fips)

na_threshold <- 0.2 # removing missing data
data_split_ratio <- 0.8 # setting up the training data
for(i in 1:4){ # for the first four targets 
  for(j in 1:length(fips)){
    tryCatch({
      subdata <- data_death[which(data_death$target == targets[i] & data_death$fips == fips[j]), ]
      na_ratio <- apply(subdata, 2, function(x) length(which(is.na(x)))/nrow(subdata))
      subdata <- subdata[ , -which(na_ratio >= na_threshold)]
      subdata <- na.omit(subdata)
      #print(paste("#Observation: ", nrow(subdata), sep = ""))
      #print(paste("#Forecasters: ", ncol(subdata)-7, sep = ""))
      sample_size_train <- ceiling(nrow(subdata) * data_split_ratio)
      data_train <- subdata[1:sample_size_train, ]
      data_test  <- subdata[(sample_size_train+1):nrow(subdata), ]
      X.train <- as.matrix(data_train[ , c(5:(ncol(data_train)-3))])
      Y.train <- as.vector(data_train$deaths)
      X.test <- as.matrix(data_test[ , c(5:(ncol(data_test)-3))])
      Y.test <- as.vector(data_test$deaths)
      
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
      
      rmse <- sqrt(colMeans((all.pred - Y.test)^2))
      print(paste(i, "-", j, sep = ""))
      write.table(t(c(i, j, cor_vec)), "emp_coviddeath_corr_new.csv", append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
      write.table(t(c(i, j, coef.stacking)), "emp_coviddeath_coef_new.csv", append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
      write.table(t(c(i, j, rmse)), "emp_coviddeath_rmse_new.csv", append = TRUE, col.names = FALSE, row.names = FALSE, sep = ",")
    }, error = function(e) {cat("Error Report!")})
  }
}

