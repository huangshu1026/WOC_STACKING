####################################################################
#### Solve the empirical results for vanilla stacking algorithm ####
####################################################################
setwd()
#### load the libraries and functions ####
source("WeightingModels.R")
source("StackingAlgorithm.R")

##### M4 data #####
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
training.index.end <- c(rep(c(11, 14, 6, 10, 5), 5), c(11, 38, 14, 6, 10, 5))
testing.index.end <- c(rep(c(14, 18, 8, 13, 6), 5), c(14, 48, 18, 8, 13, 6))
# RUNNING getM4Data for original data set
for(i in 16:20){
  # prepare the data 
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  filetag <- paste("_", item.name[i],"_",time.name[i],".csv", sep = "")
  data <- read_delim(filename, ",", col_names = TRUE)
  data <- data[ , 1:28]
  training.index <- 1:training.index.end[i]
  testing.index <- (training.index.end[i] + 1):testing.index.end[i]
  item <- unique(data$label)
  item.num <- length(unique(data$label))
  series.num <- nrow(data)/item.num # series.num - 2 > 25
  #cor_vec_record <- c()
  #coef.stacking_record <- c()
  #rmse <- c()
  for(l in 1:item.num){
    tryCatch({
      # organize the data to analyze
      subdat <- data[which(data$label ==  item[l]), ]
      training <- subdat[training.index, ]
      testing <- subdat[testing.index, ]
      X.train <- as.matrix(training[ , 4:ncol(training)])
      Y.train <- as.vector(training$target)
      X.test <- as.matrix(testing[ , 4:ncol(testing)])
      Y.test <- as.vector(testing$target)
      
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
      # cor_vec_record <- rbind(cor_vec_record, cor_vec)
      # our stacking method 
      coef.stacking <- SuperLearner(X.train, Y.train, weights, meta_index = 4, w5.select, w6.select = 9)
      #  coef.stacking_record <- rbind(coef.stacking, coef.stacking_record)
      select.wa <- which(is.na(pred.ind[1, ]) == FALSE)
      pred.stacking <- cbind(rep(1, nrow(pred.ind)), pred.ind[ , select.wa]) %*% coef.stacking
      # combine all predictions 
      all.pred <- cbind(pred.ind, pred.stacking)
      
      #rmse <- rbind(rmse, sqrt(colMeans((all.pred - Y.test)^2)))
      rmse <- sqrt(colMeans((all.pred - Y.test)^2))
      write.table(t(c(i, l, cor_vec)), file = paste("emp_m4_", item.name[i], time.name[i], "_corr.csv", sep = ""), sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(c(i, l, coef.stacking)), file = paste("emp_m4_", item.name[i], time.name[i], "_coef.csv", sep = ""), sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      write.table(t(c(i, l, rmse)), file = paste("emp_m4_", item.name[i], time.name[i], "_rmse.csv", sep = ""), sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      print(paste(i, "-", l, sep = ""))
    }, error = function(e) {cat("Error Report!")})
  }
}
