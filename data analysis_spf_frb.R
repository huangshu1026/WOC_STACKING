setwd("~/Document/weights inference/Stacking_2021/")
##### loading all the functions #####
source("WeightingModels.R")
source("UsefulFunctions.R")
source("StackingAlgorithm.R")
####################################################################################
######### comparing all weighting models by using all subsets of SPF data ##########
####################################################################################
item.name <- c(rep("NGDP", 5), rep("PGDP", 4), rep("RCONSUM", 7), rep("RGDP", 7), rep("RNRESIN", 7), rep("RRESINV", 7))
time.name <- c("currentQ", "futureQ1", "futureQ2", "futureQ3", "futureQ4", 
               "futureQ1", "futureQ2", "futureQ3", "futureQ4",
               rep(c("currentQ", "currentY", "futureQ1", "futureQ2", "futureQ3", "futureQ4", "futureY1"), 4))
M.name <- c(rep(4, 5), rep(5, 4), rep(4, 28))

rmse.all <- c()
crit <- 0.8
# for each subset of data 
for(i in 1:length(item.name)){
  filename <- paste("data/spf_usa/", item.name[i], "_", time.name[i], "_M", M.name[i], ".csv", sep = "")
  data <- read.csv(filename, stringsAsFactors = F)
  data <- data[ , -1]
  training.index <- 1:floor(nrow(data)*crit)
  testing.index <- (floor(nrow(data)*crit)+1):nrow(data)
  training <- data[training.index, ]
  testing <- data[testing.index, ]
  X <- as.matrix(training[ , 5:(ncol(training)-1)])
  y <- as.vector(training$actual)
  # getting individual weights and prediction
  weights <- WeightsCalculation(X, y)
  pred.ind <- as.matrix(testing[ , 5:(ncol(testing)-1)]) %*% weights
  # getting the stacking coefficients and prediction
  stacking <- SuperLearner(X, y, weights) 
  coef.stacking <- stacking$coef.stacking
  select.model <- stacking$model.index
  pred.stacking <- cbind(rep(1, nrow(pred.ind)), pred.ind[, select.model]) %*% coef.stacking
  # combine all predictions and calculate rmse 
  all.pred <- cbind(pred.ind, pred.stacking)
  rmse <- sqrt(colMeans((all.pred - as.vector(testing$actual))^2))
  # recording all the rmse results 
  rmse.all <- rbind(rmse.all, rmse)
  print(i)
} # change i from 1 to 37
write.csv(rmse.all, "data/rmse_spf.csv")
#write.csv(rmse.all, "spf_usa/rmse_allsets_0.7.csv")
summary(rmse.all)

###############################################
###### analyze results and visualization ######
###############################################
results <- read.csv("data/spf_usa/rmse_spf.csv")
results <- results[ , -1]
prmse <- (results[ , 1] - results)/results[ , 1]

prmse$w8.ConLas[which(prmse$w8.ConLas <= -0.2)] <- NA
prmse$w9.OurConLas[which(prmse$w9.OurConLas <= -0.2)] <- NA
prmse$X.1[which(prmse$X.1 <= -0.2)] <- NA
modelnames <- c("SA", "OW", "CWM", "SSIN", "SSDE", "RP", "TOP5", "ApproxGopt", "CLAS0", "CLASE", "Stacking")

# NGDP
apply(prmse[1:3, ], 2, function(x) mean(na.omit(x)))
# RGDP
apply(prmse[17:20, ], 2, function(x) mean(na.omit(x)))
#RCONSUM
apply(prmse[10:13, ], 2, function(x) mean(na.omit(x)))
# RNRESIN
apply(prmse[24:27, ], 2, function(x) mean(na.omit(x)))
#RRESINV
apply(prmse[31:34, ], 2, function(x) mean(na.omit(x)))

apply(prmse[c(1:3, 17:20, 10:13, 24:27, 31:34), ], 2, function(x) mean(na.omit(x)))











boxplot(prmse, ylab = "RMSE error ratio to SA", xlab = "Models", main = "Situation 1", 
        xaxt = "n", pch = 16, ylim = c(min(prmse), 0.05))
abline(h = 0, lty = 2, col = "Red")
axis(1, at = 1:length(modelnames), labels =FALSE)
text(1:length(modelnames), rep(par("usr")[3] - 0.03, length(modelnames)),
     labels = modelnames, xpd = NA, srt = 90, cex = 0.8)
text(1:length(modelnames), rep(0.045, length(modelnames)), 
     round(apply(rmse1.tosa, 2, median), 2))
#text(1:length(modelnames), rep(-2.3, length(modelnames)), 
#     round(colMeans(rmse1.tosa), 2))

boxplot(rmse2.tosa, ylab = "RMSE error ratio to SA", xlab = "Models", main = "Situation 2", 
        xaxt = "n", pch = 16, ylim = c(-2.2, 0.5))
abline(h = 0, lty = 2, col = "Red")
axis(1, at = 1:length(modelnames), labels =FALSE)
text(1:length(modelnames), rep(par("usr")[3] - 0.4, length(modelnames)),
     labels = modelnames, xpd = NA, srt = 90, cex = 0.8)
text(1:length(modelnames), rep(0.45, length(modelnames)), 
     round(apply(rmse2.tosa, 2, median), 2))
#text(1:length(modelnames), rep(-2.3, length(modelnames)), 
#     round(colMeans(rmse2.tosa), 2))

boxplot(rmse3.tosa, ylab = "RMSE error ratio to SA", xlab = "Models", main = "Situation 3", 
        xaxt = "n", pch = 16, ylim = c(-0.1, 0.08))
abline(h = 0, lty = 2, col = "Red")
axis(1, at = 1:length(modelnames), labels =FALSE)
text(1:length(modelnames), rep(par("usr")[3] - 0.03, length(modelnames)),
     labels = modelnames, xpd = NA, srt = 90, cex = 0.8)
text(1:length(modelnames), rep(0.075, length(modelnames)), 
     round(apply(rmse3.tosa, 2, median), 2))
#text(1:length(modelnames), rep(-2.3, length(modelnames)), 
#     round(colMeans(rmse3.tosa), 2))

