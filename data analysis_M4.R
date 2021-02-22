setwd("~/Document/weights inference/Stacking_2021/")
##### loading all the functions #####
source("WeightingModels.R")
source("UsefulFunctions.R")
source("StackingAlgorithm.R")
########################################
####### for time series data ###########
########################################
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))


training.index.end <- c(rep(c(11, 14, 6, 10, 5), 5), c(11, 38, 14, 6, 10, 5))
testing.index.end <- c(rep(c(14, 18, 8, 13, 6), 5), c(14, 48, 18, 8, 13, 6))

for(i in 1:length(item.name)){
  # prepare the data 
  print(paste("data/m4/", item.name[i],"_",time.name[i], sep = ""))
  filename <- paste("~/Document/weights inference/RegularizedWeights/1-data/2M4comp2018/", item.name[i],"_",time.name[i],".csv", sep = "")
  data <- read_delim(filename, ",", col_names = TRUE)
  data <- data[ , 1:28]
  training.index <- 1:training.index.end[i]
  testing.index <- (training.index.end[i] + 1):testing.index.end[i]
  item <- unique(data$label)
  item.num <- length(unique(data$label))
  series.num <- nrow(data)/item.num # series.num - 2 > 25
  # train and test models 
  rmse.allitem <- c()
  for(l in 1:item.num){
    rmse.peritem <- rep(NA, 11)
    tryCatch({
      # organize the data to analyze
      subdat <- data[which(data$label ==  item[l]), ]
      training <- subdat[training.index, ]
      testing <- subdat[testing.index, ]
      X <- as.matrix(training[ , 4:ncol(training)])
      y <- as.vector(training$target)
      # getting individual weights and prediction
      weights <- WeightsCalculation(X, y)
      pred.ind <- as.matrix(testing[ , 4:ncol(testing)]) %*% weights
      # getting the stacking coefficients and prediction
      stacking <- SuperLearner(X, y, weights) 
      coef.stacking <- stacking$coef.stacking
      select.model <- stacking$model.index
      pred.stacking <- cbind(rep(1, nrow(pred.ind)), pred.ind[, select.model]) %*% coef.stacking
      # combine all predictions and calculate rmse 
      all.pred <- cbind(pred.ind, pred.stacking)
      rmse.peritem <- sqrt(colMeans((all.pred - as.vector(testing$target))^2))
    }, error = function(e){cat("No solution for this item.")})
    rmse.allitem <- rbind(rmse.allitem, rmse.peritem)
    print(l)
  }
  # recording
  print(i)
  filetag <- paste("_", item.name[i],"_",time.name[i],".csv", sep = "")
  write.table(rmse.allitem, file = paste("data/m4/rmse", filetag, sep = ""), sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
}

###############################################
###### analyze results and visualization ######
###############################################
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
modelnames <- c("SA", "OW", "CWM", "SSIN", "SSDE", "RP", "TOP5", "ApproxGopt", "CLAS0", "CLASE", "Stacking")

results.demo <- c()
for(i in 1:5){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  results <- read.csv(paste("data/m4/rmse_", filename, sep = ""), header = FALSE)
  results.demo <- rbind(results.demo, results)
}
results.finance <- c()
for(i in 6:10){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  results <- read.csv(paste("data/m4/rmse_", filename, sep = ""), header = FALSE)
  results.finance <- rbind(results.finance, results)
}
results.industry <- c()
for(i in 11:15){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  results <- read.csv(paste("data/m4/rmse_", filename, sep = ""), header = FALSE)
  results.industry <- rbind(results.industry, results)
}
results.macro <- c()
for(i in 16:20){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  results <- read.csv(paste("data/m4/rmse_", filename, sep = ""), header = FALSE)
  results.macro <- rbind(results.macro, results)
}
results.micro <- c()
for(i in 21:25){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  results <- read.csv(paste("data/m4/rmse_", filename, sep = ""), header = FALSE)
  results.micro <- rbind(results.micro, results)
}
results.other <- c()
for(i in 26:31){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  results <- read.csv(paste("data/m4/rmse_", filename, sep = ""), header = FALSE)
  results.other <- rbind(results.other, results)
}
results.all <- rbind(results.demo, results.finance, results.industry, results.macro, results.micro, results.other)

results <- results.demo[-which(is.na(results.demo[ , 1])), ]
results <- results.finance[-which(is.na(results.finance[ , 1])), ]
results <- results.industry[-which(is.na(results.industry[ , 1])), ]
results <- results.macro[-which(is.na(results.macro[ , 1])), ]
results <- results.micro[-which(is.na(results.micro[ , 1])), ]
results <- results.other[-which(is.na(results.other[ , 1])), ]
results <- results.all[-which(is.na(results.all[ , 1])), ]
prmse <- (results[ , 1] - results)/results[ , 1]
apply(prmse, 2, function(x) median(na.omit(x)))
summary(prmse)









setEPS()
postscript(paste("m4/figures/", item.name[i], "_", time.name[i], ".eps", sep = ""))
boxplot(rmse_scaled, pch = 16, xaxt = "n", 
        ylab = "Scaled RMSE", xlab = "Models", main = paste(item.name[i], "_",time.name[i], sep = ""), 
        ylim = c(0, 0.5))
axis(1, 1:length(modelnames), labels = FALSE)
text(x = 1:length(modelnames), y = par("usr")[3] - 2*abs(par("usr")[3]),
     labels = modelnames, xpd = NA, srt = 90, cex = 0.8)
text(x = 1:length(modelnames), y = apply(rmse_scaled, 2, median) + 0.02,
     round(apply(rmse_scaled, 2, median), 3), cex = 0.8)
dev.off()

# combining data according to the data frequency
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
modelnames <- c("SA", "OW", "CWM", "SSIN", "SSDE", "RP", "CLAS0", "CLASE")
for(j in c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly")){
  rmse_scaled.all <- c()
  for(i in which(time.name == j)){
    filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
    rmse <- read.csv(paste("m4/rmse_", filename, sep = ""), header = FALSE)
    mean.testing <- read.csv(paste("m4/mean of testing/testing_", filename, sep = ""), header = F)
    rmse_scaled <- matrix(unlist(apply(rmse, 2, function(x) x/mean.testing)), ncol = ncol(rmse), nrow = nrow(rmse), byrow = FALSE)
    rmse_scaled <- na.omit(rmse_scaled)
    rmse_scaled.all <- rbind(rmse_scaled.all, rmse_scaled)
  }
  print(round(colMeans(rmse_scaled.all), 4))
}
print(nrow(rmse_scaled.all))

setEPS()
postscript(paste("m4/figures/all_", time.name[i], ".eps", sep = ""))
boxplot(rmse_scaled.all, pch = 16, xaxt = "n", 
        ylab = "Scaled RMSE", xlab = "Models", main = "", 
        ylim = c(0, 0.5))
axis(1, 1:length(modelnames), labels = FALSE)
text(x = 1:length(modelnames), y = par("usr")[3] - 2*abs(par("usr")[3]),
     labels = modelnames, xpd = NA, srt = 90, cex = 0.8)
text(x = 1:length(modelnames), y = apply(rmse_scaled.all, 2, median) + 0.02,
     round(apply(rmse_scaled.all, 2, median), 3), cex = 0.8)
dev.off()
 


