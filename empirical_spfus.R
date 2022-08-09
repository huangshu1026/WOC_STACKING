####################################################################
#### Solve the empirical results for vanilla stacking algorithm ####
####################################################################
setwd()
#### load the libraries and functions ####
source("WeightingModels.R")
source("StackingAlgorithm.R")

#####  SPF FRB data #####
item.name <- c(rep("NGDP", 5), rep("PGDP", 4), rep("RCONSUM", 7), rep("RGDP", 7), rep("RNRESIN", 7), rep("RRESINV", 7))
time.name <- c("currentQ", "futureQ1", "futureQ2", "futureQ3", "futureQ4", 
               "futureQ1", "futureQ2", "futureQ3", "futureQ4",
               rep(c("currentQ", "currentY", "futureQ1", "futureQ2", "futureQ3", "futureQ4", "futureY1"), 4))
M.name <- c(rep(4, 5), rep(5, 4), rep(4, 28))

rmse <- c()
crit <- 0.8
# for each subset of data 
for(i in 14:length(item.name)){
  tryCatch({
    filename <- paste(item.name[i], "_", time.name[i], "_M", M.name[i], ".csv", sep = "")
    data <- read.csv(filename, stringsAsFactors = F)
    data <- data[ , -1]
    training.index <- 1:floor(nrow(data)*crit)
    testing.index <- (floor(nrow(data)*crit)+1):nrow(data)
    training <- data[training.index, ]
    testing <- data[testing.index, ]
    X.train <- as.matrix(training[ , 5:(ncol(training)-1)])
    Y.train <- as.vector(training$actual)
    X.test <- as.matrix(testing[ , 5:(ncol(testing)-1)])
    Y.test <- as.vector(testing$actual)
    
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
    
    #rmse <- rbind(rmse, sqrt(colMeans((all.pred - Y.test)^2)))
    print(i)
    write.table(t(c(i, cor_vec)), "emp_spfus_corr.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    write.table(t(c(i, sqrt(colMeans((all.pred - Y.test)^2)))), "emp_spfus_rmse.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    write.table(t(c(i, coef.stacking)), "emp_spfus_coef.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }, error = function(e){cat("error report!")})
} 
#write.csv(cbind(1:length(item.name), rmse), "emp_spfus_rmse.csv")


