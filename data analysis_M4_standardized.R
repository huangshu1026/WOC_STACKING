setwd("~/Document/weights inference/LassoEW/empirical/")
##### loading all the functions #####
source("WeightingModels.R")

########################################
####### for time series data ###########
########################################
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
draws <- 10
for(i in 19:length(item.name)){
  # prepare the data 
  print(paste(item.name[i],"_",time.name[i], sep = ""))
  filename <- paste("~/Document/weights inference/RegularizedWeights/1-data/2M4comp2018/", item.name[i],"_",time.name[i],".csv", sep = "")
  data <- read_delim(filename, ",", col_names = TRUE)
  data <- data[ , 1:28]
  item <- unique(data$label)
  item.num <- length(unique(data$label))
  series.num <- nrow(data)/item.num # series.num - 2 > 25
  for(j in 1:draws){
    train.item.index <- sample(1:item.num, ceiling(item.num/2), replace = F)
    test.item.index <- c(1:item.num)[-train.item.index]
    training <- data[which(data$label %in% item[train.item.index]), ]
    testing <- data[which(data$label %in% item[test.item.index]), ]
    training.mu <- rowMeans(training[ , 4:ncol(training)])
    testing.mu <- rowMeans(testing[ , 4:ncol(testing)])
    training.sigma <- apply(training[ , 4:ncol(training)], 1, sd)
    testing.sigma <- apply(testing[ , 4:ncol(testing)], 1, sd)
    training.std.x <- (training[ , 4:ncol(training)] - training.mu)/training.sigma
    training.std.y <- (training$target - training.mu)/training.sigma
    testing.std.x <- (testing[ , 4:ncol(testing)] - testing.mu)/testing.sigma
    # estimate the weights
    weights <- WeightsCalculation(as.matrix(training.std.x), training.std.y)
    pred.std <- as.matrix(testing.std.x) %*% weights
    pred <- pred.std * testing.sigma + testing.mu    
    actual <- testing$target
    rmse <- sqrt(colMeans((pred - actual)^2))
    mape <- colMeans(abs((actual - pred)/actual))
    rmse.all <- c(item.name[i], time.name[i], rmse)
    mape.all <- c(item.name[i], time.name[i], mape)
    # recording
    write.table(t(rmse.all), file = "m4/rmse_split var/rmse_all.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
    write.table(t(mape.all), file = "m4/rmse_split var/mape_all.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
    print(j)
  }
  print(i)
}

################################
#### RESULTS VISUALIZATION #####
################################
rmse.all.mark <- read.csv("m4/rmse_split var/mape_all.csv", header = F)
modelnames <- c("SA", "OW", "CWM", "SSIN", "SSDE", "RP", "CLAS0", "CLASE")
modelpch <- c(15, 16, 17, 4, 5, 6)
# for hourly
index <- which(rmse.all.mark[ , 2] == "weekly")
plotdata <- rmse.all.mark[index, ]
plot(colMeans(plotdata[, 3:10]), pch = modelpch[6], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Other"), pch = modelpch[6])




modelnames <- c("SA", "OW", "CWM", "SSIN", "SSDE", "RP", "CLAS0", "CLASE")
modelpch <- c(15, 16, 17, 4, 5, 6)
# for hourly
index <- which(rmse.all.mark[ , 2] == "hourly")
plotdata <- rmse.all.mark[index, ]
plot(unlist(plotdata[3:10]), pch = modelpch[6], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Other"), pch = modelpch[6])
# for daily 
index <- which(rmse.all.mark[ , 2] == "daily")
plotdata <- rmse.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE", ylim = c(0, 3000))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for weekly
index <- which(rmse.all.mark[ , 2] == "weekly")
plotdata <- rmse.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE", ylim = c(0, 1500))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for monthly
index <- which(rmse.all.mark[ , 2] == "monthly")
plotdata <- rmse.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE", ylim = c(0, 5000))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topright", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for quarterly
index <- which(rmse.all.mark[ , 2] == "quarterly")
plotdata <- rmse.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE", ylim = c(0, 2500))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for yearly
index <- which(rmse.all.mark[ , 2] == "yearly")
plotdata <- rmse.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "RMSE", ylim = c(0, 2500))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("bottomright", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)




# for hourly
index <- which(mape.all.mark[ , 2] == "hourly")
plotdata <- mape.all.mark[index, ]
plot(plotdata[3:10], pch = modelpch[6], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "MAPE", ylim = c(0, 0.2))
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Other"), pch = modelpch[6])
# for daily 
index <- which(mape.all.mark[ , 2] == "daily")
plotdata <- mape.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "MAPE", ylim = c(0, 0.2))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for weekly
index <- which(mape.all.mark[ , 2] == "weekly")
plotdata <- mape.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "MAPE", ylim = c(0, 0.2))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for monthly
index <- which(mape.all.mark[ , 2] == "monthly")
plotdata <- mape.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "MAPE", ylim = c(0, 0.2))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("bottomright", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for quarterly
index <- which(mape.all.mark[ , 2] == "quarterly")
plotdata <- mape.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "MAPE", ylim = c(0, 0.2))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("topleft", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
# for yearly
index <- which(mape.all.mark[ , 2] == "yearly")
plotdata <- mape.all.mark[index, ]
plot(plotdata[1, 3:10], pch = modelpch[1], type = "b", xaxt = "n", 
     xlab = "Models", ylab = "MAPE", ylim = c(0, 0.21))
lines(plotdata[2, 3:10], pch = modelpch[2], type = "b")
lines(plotdata[3, 3:10], pch = modelpch[3], type = "b")
lines(plotdata[4, 3:10], pch = modelpch[4], type = "b")
lines(plotdata[5, 3:10], pch = modelpch[5], type = "b")
lines(plotdata[6, 3:10], pch = modelpch[6], type = "b")
axis(1, 1:8, labels = modelnames)
legend("bottomright", legend = c("Demographic", "Finance", "Industry", "Macro", "Micro", "Other"), pch = modelpch)
