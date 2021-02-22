setwd("~/Document/weights inference/Stacking_2021/")
##### loading all the functions #####
source("WeightingModels.R")
source("UsefulFunctions.R")
source("StackingAlgorithm.R")
####################################
######## data analysis #############
####################################
item <- c("wILI1", "wILI2", "wILI3", "wILI4", "peak")
location <- c("loc1", "loc2", "loc3", "loc4", "loc5", "loc6", "loc7", "loc8", "loc9", "loc10", "loc11")
rmse.all <- c()
crit <- 0.5
for(i in 1:length(item)){
  for(j in 1:length(location)){
    # prepare the data 
    filename <- paste("data/flu/", item[i], "_", location[j], ".csv", sep = "")
    data <- read.csv(filename, header = T, stringsAsFactors = F)
    training.index <- 1:floor(nrow(data)*crit)
    testing.index <- (floor(nrow(data)*crit)+1):nrow(data)
    training <- data[training.index, ]
    testing <- data[testing.index, ]
    X <- as.matrix(training[ , c(9:23, 25:30)])
    y <- as.vector(training$obs_value)
    # getting individual weights and prediction
    weights <- WeightsCalculation(X, y)
    pred.ind <- as.matrix(testing[ , c(9:23, 25:30)]) %*% weights
    # getting the stacking coefficients and prediction
    stacking <- SuperLearner(X, y, weights) 
    coef.stacking <- stacking$coef.stacking
    select.model <- stacking$model.index
    pred.stacking <- cbind(rep(1, nrow(pred.ind)), pred.ind[, select.model]) %*% coef.stacking
    # combine all predictions and calculate rmse 
    all.pred <- cbind(pred.ind, pred.stacking)
    rmse <- sqrt(colMeans((all.pred - as.vector(testing$obs_value))^2))
    rmse.all <- rbind(rmse.all, rmse)
    print(j)
  }
  print(i)
}

write.table(rmse.all, "data/rmse_flu_5_5.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
summary(rmse.all)
###############################################
###### analyze results and visualization ######
###############################################
results <- read.csv("data/flu/rmse_flu.csv", header = T)
# visualization
modelnames <- c("SA", "OW", "CWM", "SSIN", "SSDE", "RP", "TOP5", "ApproxGopt", "CLAS0", "CLASE", "Stacking")
prmse <- (results[ , 1] - results)/results[ , 1]
# wILI
round(apply(prmse[1:11, ], 2, function(x) mean(na.omit(x))), 4)
round(apply(prmse[12:22, ], 2, function(x) mean(na.omit(x))), 4)
round(apply(prmse[23:33, ], 2, function(x) mean(na.omit(x))), 4)
round(apply(prmse[34:44, ], 2, function(x) mean(na.omit(x))), 4)
round(apply(prmse[1:44, ], 2, function(x) mean(na.omit(x))), 4)

summary(prmse)
boxplot(prmse)






boxplot(rmse.all, pch = 16, xaxt = "n", 
        ylab = "RMSE", xlab = "Models")
#abline(h = min(apply(rmse.all, 2, median)), col = "red", lty = 2)
axis(1, 1:ncol(rmse.all), labels = FALSE)
text(x = 1:ncol(rmse.all), y = par("usr")[3] - 0.4,
     labels = modelnames, xpd = NA, srt = 90, cex = 0.8)

# divided by variables 
for(i in 1:4){
  rmse.sub <- rmse.all[(11*(i-1)+1):(11*i), ]
  #rmse.tosa.sub <- rmse.tosa[(11*(i-1)+1):(11*i), ]
  #boxplot(rmse.tosa.sub, col = colors, pch = 20, xaxt = "n", 
  #        ylab = "Error Ratio of RMSE to SA", xlab = "Models",ylim = c(-1, 0.6), main = item[i])
  boxplot(rmse.sub, pch = 16, xaxt = "n", 
          ylab = "RMSE", xlab = "Models",ylim = c(0, 3))
  abline(v = 8.5, lty = 1, lwd = 2)
  #abline(h = 0, col = "red", lty = 2)
  axis(1, 1:ncol(rmse.sub), labels = FALSE)
  text(x = 1:ncol(rmse.sub), y = par("usr")[3] - 0.5,
       labels = modelnames, xpd = NA, srt = 90, cex = 0.8)
  text(1:ncol(rmse.sub), rep(2.9, ncol(rmse.sub)), round(apply(rmse.sub, 2, median), 2), cex = 0.8)
  #print(which.min(colMeans(rmse.sub)))
}




