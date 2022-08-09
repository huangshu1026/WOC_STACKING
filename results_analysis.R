setwd()
##### simulation results: inflation rate base simulation #####
rmse_inf <- read.csv("sim_inf_rmse.csv", header = FALSE)
colnames(rmse_inf) <- c("SampleSize", "RandomDraw", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2", "Stacking")
mse_inf <- rmse_inf^2
# get the average individual rmse 
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
AvgInd <- c()
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
    
    AvgInd <- c(AvgInd, mean(colMeans((X.test - Y.test)^2)))
  }
}
mse_inf <- cbind(mse_inf, AvgInd)
colMeans(mse_inf)
samplesize <- unique(mse_inf$SampleSize)
avg_mse_ratio <- c()
med_mse_ratio <- c()
for(i in 1:length(samplesize)){
  sub_mse_inf <- mse_inf[which(mse_inf$SampleSize == samplesize[i]), ]
  sub_mse_ratio <- (sub_mse_inf$AvgInd - sub_mse_inf)/sub_mse_inf$AvgInd
  sub_mse_ratioavg <- colMeans(sub_mse_ratio)[3:ncol(sub_mse_ratio)]
  sub_mse_ratiomed <- unlist(apply(sub_mse_ratio, 2, function(x) median(x)))[3:ncol(sub_mse_ratio)]
  avg_mse_ratio <- rbind(avg_mse_ratio, sub_mse_ratioavg)
  med_mse_ratio <- rbind(med_mse_ratio, sub_mse_ratiomed)
}
print(round(med_mse_ratio, 4)*100)
print(round(avg_mse_ratio, 4)*100)
# coefficients
coef_inf <- read.csv("sim_inf_coef.csv", header = F)
colnames(coef_inf) <- c("SampleSize", "RandomDraw", "Cont", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
head(coef_inf)
samplesize <- unique(coef_inf$SampleSize)
avg_coef <- c()
for(i in 1:length(samplesize)){
  temp <- coef_inf[which(coef_inf$SampleSize == samplesize[i]), ]
  temp <- colMeans(temp, na.rm = TRUE)
  avg_coef <- rbind(avg_coef, temp[4:length(temp)])
}
print(round(avg_coef, 3))
# correlation 
corr_inf <- read.csv("sim_inf_corr.csv", header = F)
samplesize <- unique(corr_inf$V1)
avg_corr <- array(NA, dim = c(12, 12, 5))
for(i in 1:length(samplesize)){
  temp <- corr_inf[which(corr_inf$V1 == samplesize[i]), ]
  temp <- colMeans(temp)
  corr_mat <- matrix(NA, nrow = 12, ncol = 12)
  corr_mat[upper.tri(corr_mat)] <- temp[3:length(temp)]
  corr_mat[lower.tri(corr_mat)] <- t(corr_mat)[lower.tri(t(corr_mat))]
  diag(corr_mat) <- 1
  avg_corr[ , , i] <- corr_mat
}
corr_mat_plot <- avg_corr[ , , 5]
colnames(corr_mat_plot) <- rownames(corr_mat_plot) <- c("SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
corrplot::corrplot(corr_mat_plot, method = "color", type="upper", addCoef.col = "darkgray", tl.col="black")

##### simulation results: unemployment rate base simulation #####
rmse_unemp <- read.csv("sim_unemp_rmse.csv", header = FALSE)
colnames(rmse_unemp) <- c("SampleSize", "RandomDraw", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2", "Stacking")
mse_unemp <- rmse_unemp^2
# get the average individual rmse 
unempdata <- read.csv("unemp_filtered.csv", header = T, sep = ",")
unempdata <- unempdata[ , -1]
M <- 8
X <- as.matrix(unempdata[ , 4:11]) 
Y <- unempdata[ , 3]
## estimate the bias and covariance matrix -> simulation seed
muy.emp <- mean(Y)
sdy.emp <- sd(Y)
mu.emp <- colMeans(X - Y)
sigma.emp <- cov(X)
## parameter settings for simulation 
n.vec.train <- 2^seq(4, 8, by = 1)
n.test <- 1000
draws <- 100
AvgInd <- c()
for(j in 1:length(n.vec.train)){
  n <- n.vec.train[j]
  for(k in 1:draws){
    set.seed(j*k*4718)
    
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
    
    AvgInd <- c(AvgInd, mean(colMeans((X.test - Y.test)^2)))
  }
}
mse_unemp <- cbind(mse_unemp, AvgInd)
colMeans(mse_unemp)
samplesize <- unique(mse_unemp$SampleSize)
avg_mse_ratio <- c()
med_mse_ratio <- c()
for(i in 1:length(samplesize)){
  sub_mse_unemp <- mse_unemp[which(mse_unemp$SampleSize == samplesize[i]), ]
  sub_mse_ratio <- (sub_mse_unemp$AvgInd - sub_mse_unemp)/sub_mse_unemp$AvgInd
  sub_mse_ratioavg <- colMeans(sub_mse_ratio)[3:ncol(sub_mse_ratio)]
  sub_mse_ratiomed <- unlist(apply(sub_mse_ratio, 2, function(x) median(x)))[3:ncol(sub_mse_ratio)]
  avg_mse_ratio <- rbind(avg_mse_ratio, sub_mse_ratioavg)
  med_mse_ratio <- rbind(med_mse_ratio, sub_mse_ratiomed)
}
print(round(med_mse_ratio, 4)*100)
print(round(avg_mse_ratio, 4)*100)

# coefficients
coef_unemp <- read.csv("sim_unemp_coef.csv", header = F)
colnames(coef_unemp) <- c("SampleSize", "RandomDraw", "Cont", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
head(coef_unemp)
samplesize <- unique(coef_unemp$SampleSize)
avg_coef <- c()
for(i in 1:length(samplesize)){
  temp <- coef_unemp[which(coef_unemp$SampleSize == samplesize[i]), ]
  temp <- colMeans(temp, na.rm = TRUE)
  avg_coef <- rbind(avg_coef, temp[4:length(temp)])
}
print(round(avg_coef, 3))
# correlation 
corr_unemp <- read.csv("sim_unemp_corr.csv", header = F)
samplesize <- unique(corr_unemp$V1)
avg_corr <- array(NA, dim = c(12, 12, 5))
for(i in 1:length(samplesize)){
  temp <- corr_unemp[which(corr_unemp$V1 == samplesize[i]), ]
  temp <- colMeans(temp)
  corr_mat <- matrix(NA, nrow = 12, ncol = 12)
  corr_mat[upper.tri(corr_mat)] <- temp[3:length(temp)]
  corr_mat[lower.tri(corr_mat)] <- t(corr_mat)[lower.tri(t(corr_mat))]
  diag(corr_mat) <- 1
  avg_corr[ , , i] <- corr_mat
}
corr_mat_plot <- avg_corr[ , , 5]
colnames(corr_mat_plot) <- rownames(corr_mat_plot) <- c("SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
corrplot::corrplot(corr_mat_plot, method = "color", type="upper", addCoef.col = "darkgray", tl.col="black")

##### simulation results: covid19 cases based simulation #####
rmse_case1 <- read.csv("1sim_covidcase_rmse.csv", header = FALSE)
rmse_case2 <- read.csv("2sim_covidcase_rmse.csv", header = FALSE)
rmse_case3 <- read.csv("3sim_covidcase_rmse.csv", header = FALSE)
rmse_case4 <- read.csv("4sim_covidcase_rmse.csv", header = FALSE)
rmse_case <- rbind(rmse_case1, rmse_case2, rmse_case3, rmse_case4)
nrow(rmse_case)/5
head(rmse_case)
colnames(rmse_case) <- c("target", "location", "samplesize", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2", "Stacking")
mse_all <- rmse_case^2
# get the average individual mse 
data_case <- read.csv("all_data_case_county.csv", header = T, stringsAsFactors = F)
data_case <- data_case[ , -1] #1:7, 66:70 are indicators
targets <- unique(data_case$target)
fips <- unique(data_case$fips)
na_threshold <- 0.1
n.vector <- 2^seq(4, 8, by = 1)
n.test <- 1000
mse_indAvg <- c()
for(i in 1:4){ # the first four targets
  for(j in 1:length(fips)){
    dat <- data_case[which(data_case$target == targets[i] & data_case$fips == fips[j]), ]
    na_ratio <- apply(dat, 2, function(x) length(which(is.na(x)))/nrow(dat))
    dat <- dat[ , -which(na_ratio >= na_threshold)]
    dat <- na.omit(dat)
    #print(paste("#Observation: ", nrow(dat), sep = ""))
    #print(paste("#Forecasters: ", ncol(dat)-11, sep = ""))
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
          
          mse_indAvg <- rbind(mse_indAvg, c(i, j, n.vector[r], mean(colMeans((X.test - Y.test)^2))))
          
        }, error = function(e){ cat("One Error Found!") })
      }
    }
  }
}
mse_all$target <- sqrt(mse_all$target)
mse_all$location <- sqrt(mse_all$location)
mse_all$samplesize <- sqrt(mse_all$samplesize)
mse_indAvg <- data.frame(mse_indAvg)
colnames(mse_indAvg) <- c("target", "location", "samplesize", "IndAvg")
alloutput <- merge(mse_all, mse_indAvg, by = c("target", "location", "samplesize"), all.x = TRUE)
alloutput_ratio <- (alloutput$IndAvg - alloutput)/alloutput$IndAvg
samplesize <- unique(alloutput$samplesize)
avg_mse_ratio <- c()
med_mse_ratio <- c()
for(i in 1:length(samplesize)){
  sub_mse <- alloutput[which(alloutput$samplesize == samplesize[i]), ]
  sub_mse_ratio <- (sub_mse$IndAvg - sub_mse)/sub_mse$IndAvg
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(sub_mse_ratio, 2, function(x) mean(x, na.rm = TRUE)))[4:ncol(sub_mse_ratio)])
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(sub_mse_ratio, 2, function(x) median(x, na.rm = TRUE)))[4:ncol(sub_mse_ratio)])
}
print(round(med_mse_ratio, 4)*100)
print(round(avg_mse_ratio, 4)*100)

# coefficients
coef_case1 <- read.csv("1sim_covidcase_coef.csv", header = F)
coef_case2 <- read.csv("2sim_covidcase_coef.csv", header = F)
coef_case3 <- read.csv("3sim_covidcase_coef.csv", header = F)
coef_case4 <- read.csv("4sim_covidcase_coef.csv", header = F)
coef_case <- rbind(coef_case1, coef_case2, coef_case3, coef_case4)
head(coef_case)
colnames(coef_case) <- c("target", "location", "SampleSize", "Cont", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
samplesize <- unique(coef_case$SampleSize)
avg_coef <- c()
for(i in 1:length(samplesize)){
  temp <- coef_case[which(coef_case$SampleSize == samplesize[i]), ]
  temp <- colMeans(temp, na.rm = TRUE)
  avg_coef <- rbind(avg_coef, temp[4:length(temp)])
}
print(round(avg_coef, 3))
# correlation 
corr_case1 <- read.csv("1sim_covidcase_corr.csv", header = F)
corr_case2 <- read.csv("2sim_covidcase_corr.csv", header = F)
corr_case3 <- read.csv("3sim_covidcase_corr.csv", header = F)
corr_case4 <- read.csv("4sim_covidcase_corr.csv", header = F)
corr_case <- rbind(corr_case1, corr_case2, corr_case3, corr_case4)
head(corr_case)
samplesize <- unique(corr_case$V3)
avg_corr <- array(NA, dim = c(12, 12, 5))
for(i in 1:length(samplesize)){
  temp <- corr_case[which(corr_case$V3 == samplesize[i]), ]
  temp <- colMeans(temp, na.rm = TRUE)
  corr_mat <- matrix(NA, nrow = 12, ncol = 12)
  corr_mat[upper.tri(corr_mat)] <- temp[4:length(temp)]
  corr_mat[lower.tri(corr_mat)] <- t(corr_mat)[lower.tri(t(corr_mat))]
  diag(corr_mat) <- 1
  avg_corr[ , , i] <- corr_mat
}
corr_mat_plot <- avg_corr[ , , 5]
colnames(corr_mat_plot) <- rownames(corr_mat_plot) <- c("SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
corrplot::corrplot(corr_mat_plot, method = "color", type="upper", addCoef.col = "darkgray", tl.col="black")

##### empirical analysis: m4 competition #####
# get the average individual mse 
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
training.index.end <- c(rep(c(11, 14, 6, 10, 5), 5), c(11, 38, 14, 6, 10, 5))
testing.index.end <- c(rep(c(14, 18, 8, 13, 6), 5), c(14, 48, 18, 8, 13, 6))
## find the avarage out-of-sample individual MSEs
mse_indAvg <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  filename <- paste(item.name[i],"_",time.name[i],".csv", sep = "")
  filetag <- paste("_", item.name[i],"_",time.name[i],".csv", sep = "")
  data <- read_delim(filename, ",", col_names = TRUE)
  data <- data[ , 1:28]
  training.index <- 1:training.index.end[i]
  testing.index <- (training.index.end[i] + 1):testing.index.end[i]
  item <- unique(data$label)
  item.num <- length(unique(data$label))
  series.num <- nrow(data)/item.num # series.num - 2 > 25
  mse_ind <- c()
  for(l in 1:item.num){
    subdat <- data[which(data$label ==  item[l]), ]
    testing <- subdat[testing.index, ]
    testing.mat <- as.matrix(testing[ , c(3:ncol(testing))])
    mse_ind[l] <- mean(colMeans((testing.mat[ , c(2:ncol(testing.mat))] - testing.mat[ , 1])^2, na.rm = TRUE))
  }
  mse_indAvg <- rbind(mse_indAvg, cbind(rep(i, item.num), 1:item.num, mse_ind))
}
colnames(mse_indAvg) <- c("file", "item", "mse_ind")
dim(mse_indAvg)
# get the rmse of all the weighting methods 
rmse_all <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  temp <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_rmse.csv", sep = ""), header = F)
  rmse_all <- rbind(rmse_all, temp)
}
colnames(rmse_all) <- c("file", "item", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2", "Stacking")
mse_all <- rmse_all^2
mse_all$file <- rmse_all$file
mse_all$item <- rmse_all$item
mse_all <- merge(mse_all, mse_indAvg, by = c("file", "item"), all.x = TRUE)
# compare all weighting methods 
#target_time <- c("hourly", "daily", "weekly", "monthly", "quarterly")
target_time <- c("demo", "finance", "industry", "macro", "micro", "other")
avg_mse_ratio <- c()
med_mse_ratio <- c()
mse_ratio_all <- c()
avg_mse <- c()
med_mse <- c()
for(j in 1:length(target_time)){
  #mse_target <- mse_all[which(mse_all$file %in% which(time.name == target_time[j])), ]
  mse_target <- mse_all[which(mse_all$file %in% which(item.name == target_time[j])), ]
  mse_ratio <- c()
  #for(i in which(time.name == target_time[j])){
  for(i in which(item.name == target_time[j])){
    sub_mse <- mse_all[which(mse_all$file %in% i), ]
    sub_mse_ratio <- (sub_mse$mse_ind - sub_mse[ , 3:ncol(sub_mse)])/sub_mse$mse_ind
    mse_ratio <- rbind(mse_ratio, sub_mse_ratio)
    mse_ratio_all <- rbind(mse_ratio_all, sub_mse_ratio)
  }
  avg_mse <- rbind(avg_mse, unlist(apply(mse_target, 2, function(x) mean(x, na.rm = T))))
  med_mse <- rbind(med_mse, unlist(apply(mse_target, 2, function(x) median(x, na.rm = T))))
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(mse_ratio, 2, function(x) mean(x, na.rm = T))))
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(mse_ratio, 2, function(x) median(x, na.rm = T))))
}
print(round(med_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))
print(round(avg_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))

print(round(med_mse/100000, 3))
apply(mse_all, 2, function(x) median(x, na.rm = TRUE))/100000
print(round(avg_mse/100000, 3))
apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))/100000

# coefficients
coef_all <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  temp <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_coef.csv", sep = ""), header = F)
  coef_all <- rbind(coef_all, temp)
}
colnames(coef_all) <- c("file", "item","Cont", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
target_time <- c("hourly", "daily", "weekly", "monthly", "quarterly")
avg_coef <- c()
for(j in 1:length(target_time)){
  coef_temp <- coef_all[which(coef_all$file %in% which(time.name == target_time[j])), ]
  avg_coef <- rbind(avg_coef, colMeans(coef_temp, na.rm = TRUE))
}
print(round(avg_coef, 2))
# correlation 
corr_all <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  temp <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_corr.csv", sep = ""), header = F)
  corr_all <- rbind(corr_all, temp)
}
corr_avg <- colMeans(corr_all, na.rm = TRUE)
corr_mat_plot <- matrix(NA, nrow = 12, ncol = 12)
corr_mat_plot[upper.tri(corr_mat_plot)] <- corr_avg[3:length(corr_avg)]
corr_mat_plot[lower.tri(corr_mat_plot)] <- t(corr_mat_plot)[lower.tri(t(corr_mat_plot))]
diag(corr_mat_plot) <- 1
colnames(corr_mat_plot) <- rownames(corr_mat_plot) <- c("SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
corrplot::corrplot(corr_mat_plot, method = "color", type="upper", addCoef.col = "darkgray", tl.col="black")

##### empirical analysis: SPF USA data #####
rmse_all <- read.csv("emp_spfus_rmse.csv", header = F)
head(rmse_all)
colnames(rmse_all) <- c("item", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2", "Stacking")
mse_all <- rmse_all^2
mse_all$item <- rmse_all$item
# get the average individual mse 
item.name <- c(rep("NGDP", 5), rep("PGDP", 4), rep("RCONSUM", 7), rep("RGDP", 7), rep("RNRESIN", 7), rep("RRESINV", 7))
time.name <- c("currentQ", "futureQ1", "futureQ2", "futureQ3", "futureQ4", 
               "futureQ1", "futureQ2", "futureQ3", "futureQ4",
               rep(c("currentQ", "currentY", "futureQ1", "futureQ2", "futureQ3", "futureQ4", "futureY1"), 4))
M.name <- c(rep(4, 5), rep(5, 4), rep(4, 28))
crit <- 0.8
mse_indAvg <- c()
for(i in 1:length(item.name)){
  filename <- paste(item.name[i], "_", time.name[i], "_M", M.name[i], ".csv", sep = "")
  data <- read.csv(filename, stringsAsFactors = F)
  data <- data[ , -1]
  testing.index <- (floor(nrow(data)*crit)+1):nrow(data)
  testing <- data[testing.index, ]
  testing.mat <- as.matrix(testing[ , c(5:ncol(testing))])
  mse_indAvg <- rbind(mse_indAvg, c(i, mean(colMeans((testing.mat[ , c(1:(ncol(testing.mat)-1))] - testing.mat[ , ncol(testing.mat)])^2, na.rm = TRUE))))
}
colnames(mse_indAvg) <- c("item", "IndAvg")
mse_all <- merge(mse_all, mse_indAvg, by= "item", all.y = TRUE)
mse_all <- mse_all[ , -1]
#colnames(mse_all) <- c("SA", "OW", "CovW", "VarW", "CCor", "Reg", "CWM", "SSin", "SSde", "RP", "Top3", "Stacking", "VanillaStacking", "IndAvg")
mse_ratio_all <- (mse_all$IndAvg - mse_all)/mse_all$IndAvg
# get the aggregate results --> in tables 
#target_time <- c("currentQ", "futureQ1", "futureQ2", "futureQ3", "futureQ4", "currentY", "futureY1")
target_time <- c("RCONSUM", "RGDP", "RRESINV", "RNRESIN", "NGDP", "PGDP")
avg_mse_ratio <- c()
med_mse_ratio <- c()
avg_mse <- c()
med_mse <- c()
for(j in 1:length(target_time)){
  #rowindex <- which(time.name == target_time[j])
  rowindex <- which(item.name == target_time[j])
  avg_mse <- rbind(avg_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse <- rbind(med_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
}
print(round(med_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))*100
output <- rbind(round(med_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")

print(round(avg_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))*100

print(round(med_mse/10000, 3))
apply(mse_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse/10000, 3), apply(mse_all, 2, function(x) median(x, na.rm = TRUE))/10000)
write.csv(output, "tempfile.csv")

print(round(avg_mse/10000, 3))
apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse/10000, 3), apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))/10000)
write.csv(output, "tempfile.csv")

coef_spfus <- read.csv("emp_spfus_coef.csv", header = F)
colnames(coef_spfus) <- c("item", "Cont", "SA", "OW", "CovW", "VarW", "CCorW", "VanillaS1", "CWM", "SSIN", "SSDE", "RP", "TOP3", "VanillaS2")
target_time <- c("currentQ", "futureQ1", "futureQ2", "futureQ3", "futureQ4", "currentY", "futureY1")
avg_coef <- c()
for(i in 1:length(target_time)){
  rowindex <- which(coef_spfus$item %in% which(time.name == target_time[i]))
  avg_coef <- rbind(avg_coef, colMeans(coef_spfus[rowindex, ], na.rm = TRUE))
}
write.csv(avg_coef, "tempfile.csv")
##### empirical analysis: Flu forecasting #####
rmse_all <- read.csv("emp_flu_rmse.csv", header = F)
colnames(rmse_all) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "Reg", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2", "Stacking")
mse_all <- rmse_all^2
mse_all$Target <- rmse_all$Target
mse_all$Location <- rmse_all$Location
# get the average individual mse 
item <- c("wILI1", "wILI2", "wILI3", "wILI4", "peak")
location <- c("loc1", "loc2", "loc3", "loc4", "loc5", "loc6", "loc7", "loc8", "loc9", "loc10", "loc11")
crit <- 0.8
mse_indAvg <- c()
for(i in 1:length(item)){
  for(j in 1:length(location)){
    filename <- paste(item[i], "_", location[j], ".csv", sep = "")
    data <- read.csv(filename, header = T, stringsAsFactors = F)
    testing.index <- (floor(nrow(data)*crit)+1):nrow(data)
    testing <- data[testing.index, ]
    testing.mat <- as.matrix(testing[ , c(8, 9:23, 25:30)])
    mse_indAvg <- c(mse_indAvg, mean(colMeans((testing.mat[ , c(2:ncol(testing.mat))] - testing.mat[ , 1])^2, na.rm = TRUE)))
  }
}
# get MSE of all other weighting methods 
mse_all <- cbind(mse_all, mse_indAvg[1:44])
mse_all <- mse_all[ , -c(1:2)]
colnames(mse_all) <- c("SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2", "Stacking", "IndAvg")
mse_ratio_all <- (mse_all$IndAvg - mse_all)/mse_all$IndAvg
# get the aggregate level statistics 
#target_time <- c("wILI1", "wILI2", "wILI3", "wILI4")
target_time <- c("loc1", "loc2", "loc3", "loc4", "loc5", "loc6", "loc7", "loc8", "loc9", "loc10", "loc11")
avg_mse_ratio <- c()
med_mse_ratio <- c()
avg_mse <- c()
med_mse <- c()
for(j in 1:length(target_time)){
  #rowindex <- ((j-1)*11+1):(j*11)
  rowindex <- c(j, j+11, j+22, j+33)
  avg_mse <- rbind(avg_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse <- rbind(med_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
}
print(round(med_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")

print(round(avg_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")

print(round(med_mse, 3))
apply(mse_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse, 3), apply(mse_all, 2, function(x) median(x, na.rm = TRUE)))
write.csv(output, "tempfile.csv")

print(round(avg_mse, 3))
apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse, 3), apply(mse_all, 2, function(x) mean(x, na.rm = TRUE)))
write.csv(output, "tempfile.csv")

coef_flu <- read.csv("emp_flu_coef.csv", header = F)
coefs <- coef_flu$V3
coefs <- matrix(coefs, nrow = 44, ncol = 13, byrow = TRUE)
coef_flu <- cbind(rep(1:4, each = 11), rep(1:11, 4), coefs)
coef_flu <- as.data.frame(coef_flu)
colnames(coef_flu) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2")
avg_coef <- c()
for(i in 1:4){
  rowindex <- which(coef_flu$Target == i)
  avg_coef <- rbind(avg_coef, colMeans(coef_flu[rowindex, ], na.rm = TRUE))
}
write.csv(avg_coef, "tempfile.csv")
##### empirical analysis: COVID 19 death forecasting #####
rmse_all <- read.csv("emp_coviddeath_rmse_new.csv", header = F)
colnames(rmse_all) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2", "Stacking")
mse_all <- rmse_all^2
mse_all$Target <- rmse_all$Target
mse_all$Location <- rmse_all$Location
# get the average individual mse 
data_death <- read.csv("all_data_death.csv", header = T, stringsAsFactors = F)
data_death <- data_death[ , -1] #1:4, 105:107 are indicators
# grouped data 
targets <- unique(data_death$target)
fips <- unique(data_death$fips)
na_threshold <- 0.2 # removing missing data
data_split_ratio <- 0.8 # setting up the training data
mse_indAvg <- c()
for(i in 1:4){
  mse_ind <- c()
  for(j in 1:length(fips)){
    subdata <- data_death[which(data_death$target == targets[i] & data_death$fips == fips[j]), ]
    na_ratio <- apply(subdata, 2, function(x) length(which(is.na(x)))/nrow(subdata))
    subdata <- subdata[ , -which(na_ratio >= na_threshold)]
    subdata <- na.omit(subdata)
    sample_size_train <- ceiling(nrow(subdata) * data_split_ratio)
    data_test  <- subdata[(sample_size_train+1):nrow(subdata), ]
    testing.mat <- as.matrix(data_test[ , c(5:(ncol(data_test)-3), ncol(data_test))])
    if(ncol(testing.mat) > 2){
      mse_ind[j] <- mean(colMeans((testing.mat[ , c(1:(ncol(testing.mat)-1))] - testing.mat[ , ncol(testing.mat)])^2, na.rm = TRUE))
    }
    if(ncol(testing.mat) == 2){
      mse_ind[j] <- mean((testing.mat[ , c(1:(ncol(testing.mat)-1))] - testing.mat[ , ncol(testing.mat)])^2)
    }
  }
  mse_indAvg <- rbind(mse_indAvg, cbind(rep(i, length(fips)), 1:length(fips), mse_ind))
}
# get mse for all weighting methods 
colnames(mse_indAvg) <- c("Target", "Location", "IndAvg")
mse_all <- merge(mse_all, mse_indAvg, by = c("Target", "Location"), all.y = TRUE)
#mse_all <- mse_all[ , -c(1:2)]
mse_all <- mse_all[-which(mse_all$IndAvg == 0), ]
#colnames(mse_all) <- c("SA", "OW", "CovW", "VarW", "CCor", "Reg", "CWM", "SSin", "SSde", "RP", "Top3", "Stacking", "VanillaStacking", "IndAvg")
mse_ratio_all <- (mse_all$IndAvg - mse_all)/mse_all$IndAvg
# summarize the aggregate statistics 
avg_mse_ratio <- c()
med_mse_ratio <- c()
avg_mse <- c()
med_mse <- c()
for(j in 1:4){
  rowindex <- which(mse_all$Target == j)
  #rowindex <- ((j-1)*55+1):(j*55)
  avg_mse <- rbind(avg_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse <- rbind(med_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
}
output <- rbind(round(med_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")
print(round(avg_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")

print(round(med_mse/10000, 3))
apply(mse_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse/10000, 3), apply(mse_all, 2, function(x) median(x, na.rm = TRUE))/10000)
write.csv(output, "tempfile.csv")
print(round(avg_mse/10000, 3))
apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse/10000, 3), apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))/10000)
write.csv(output, "tempfile.csv")

# coefficients
coef_death <- read.csv("emp_coviddeath_coef_new.csv", header = F)
head(coef_death)
avg_coef <- c()
for(i in 1:4){
  temp <- coef_death[which(coef_death$V1 == i), ]
  avg_coef <- rbind(avg_coef, colMeans(temp, na.rm = TRUE))
}
colnames(avg_coef) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2")
write.csv(avg_coef, "tempfile.csv")

##### empirical analysis: COVID 19 case forecasting #####
rmse_case1 <- read.csv("emp_covidcase1_rmse.csv", header = F)
rmse_case2 <- read.csv("emp_covidcase2_rmse.csv", header = F)
rmse_case3 <- read.csv("emp_covidcase3_rmse.csv", header = F)
rmse_all <- rbind(rmse_case1, rmse_case2, rmse_case3)
head(rmse_all)
colnames(rmse_all) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2", "Stacking")
mse_all <- rmse_all^2
mse_all$Target <- rmse_all$Target
mse_all$Location <- rmse_all$Location
# get the average individual mse 
data_case <- read.csv("all_data_case_county.csv", header = T, stringsAsFactors = F)
data_case <- data_case[ , -1] #1:7, 66:70 are indicators
na_ratio <- apply(data_case, 2, function(x) length(which(is.na(x)))/nrow(data_case))
targets <- unique(data_case$target)
fips <- unique(data_case$fips)
length(fips)
na_threshold <- 0.2 # removing missing data
data_split_ratio <- 0.8 # setting up the training data
mse_indAvg <- c()
for(i in 1:4){
  mse_ind <- c()
  for(j in 1:length(fips)){
    subdata <- data_case[which(data_case$target == targets[i] & data_case$fips == fips[j]), ]
    na_ratio <- apply(subdata, 2, function(x) length(which(is.na(x)))/nrow(subdata))
    subdata <- subdata[ , -which(na_ratio >= na_threshold)]
    subdata <- na.omit(subdata)
    sample_size_train <- ceiling(nrow(subdata) * data_split_ratio)
    data_test  <- subdata[(sample_size_train+1):nrow(subdata), ]
    testing.mat <- as.matrix(data_test[ , c(8:(ncol(data_test)-5), ncol(data_test))])
    if(ncol(testing.mat) > 2){
      mse_ind[j] <- mean(colMeans((testing.mat[ , c(1:(ncol(testing.mat)-1))] - testing.mat[ , ncol(testing.mat)])^2, na.rm = TRUE))
    }
    if(ncol(testing.mat) == 2){
      mse_ind[j] <- mean((testing.mat[ , c(1:(ncol(testing.mat)-1))] - testing.mat[ , ncol(testing.mat)])^2)
    }
  }
  mse_indAvg <- rbind(mse_indAvg, cbind(rep(i, length(fips)), 1:length(fips), mse_ind))
}

# get mse for all weighting methods 
mse_indAvg <- as.data.frame(mse_indAvg)
colnames(mse_indAvg) <- c("Target", "Location", "IndAvg")
mse_all <- merge(mse_all, mse_indAvg, by = c("Target", "Location"), all.y = TRUE)
mse_all <- mse_all[-which(mse_all$IndAvg == 0), ]
#mse_all <- mse_all[ , -c(1:2)]
mse_ratio_all <- (mse_all$IndAvg - mse_all)/mse_all$IndAvg
# summarize the aggregate statistics 
avg_mse_ratio <- c()
med_mse_ratio <- c()
avg_mse <- c()
med_mse <- c()
for(j in 1:4){
  #rowindex <- ((j-1)*3132+1):(j*3132)
  rowindex <- which(mse_all$Target == j)
  avg_mse <- rbind(avg_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse <- rbind(med_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
}
print(round(med_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")
print(round(avg_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")

print(round(med_mse/10000, 3))
apply(mse_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse/10000, 3), apply(mse_all, 2, function(x) median(x, na.rm = TRUE))/10000)
write.csv(output, "tempfile.csv")
print(round(avg_mse/10000, 3))
apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse/10000, 3), apply(mse_all, 2, function(x) mean(x, na.rm = TRUE))/10000)
write.csv(output, "tempfile.csv")

# coefficients 
coef_case1 <- read.csv("emp_covidcase1_coef.csv", header = F)
coef_case2 <- read.csv("emp_covidcase3_coef.csv", header = F)
coef_all <- rbind(coef_case1, coef_case2)
head(coef_all)
avg_coef <- c()
for(i in 1:4){
  temp <- coef_all[which(coef_all$V1 == i), ]
  avg_coef <- rbind(avg_coef, colMeans(temp, na.rm = TRUE))
}
colnames(avg_coef) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2")
write.csv(avg_coef, "tempfile.csv")

##### Empirical analysis: Keck and Tang (2020) #####
rmse_all <- read.csv("emp_keck_rmse.csv", header = F)
colnames(rmse_all) <- c("condition", "samplesize", "randomdraw", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2", "Stacking")
mse_all <- rmse_all^2
mse_all$condition <- rmse_all$condition
mse_all$samplesize <- rmse_all$samplesize
mse_all$randomdraw <- rmse_all$randomdraw
# get average individual mse
library(readxl)
exp_dat <- read_excel("DataStudy1.xlsx")
# Start the anlaysis 
M <- 5 # 48 forecasters with 40 forecasts 
n.train.vec <- c(10, 15, 20, 25, 30)
draws <- 100
condition <- unique(exp_dat$con) # 3 conditions
mse_indAvg <- c()
for(i in 1:3){ # index of condition
  for(j in 1:length(n.train.vec)){
    n <- n.train.vec[j]
    for(k in 1:draws){
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
      
      mse_indAvg <- rbind(mse_indAvg, c(i, n, k, mean(colMeans((X.test - Y.test)^2))))
    }
  }
}
mse_indAvg <- as.data.frame(mse_indAvg)
colnames(mse_indAvg) <- c("condition", "samplesize", "randomdraw", "mse_indAvg")
mse_all <- merge(mse_all, mse_indAvg, by = c("condition", "samplesize", "randomdraw"), all.y = TRUE)
mse_ratio_all <- (mse_all$mse_indAvg - mse_all)/mse_all$mse_indAvg
# summarize the aggregate statistics 
#samplesize <- unique(mse_all$samplesize)
condition <- unique(mse_all$condition)
avg_mse_ratio <- c()
med_mse_ratio <- c()
avg_mse <- c()
med_mse <- c()
#for(j in 1:length(samplesize)){
for(j in 1:length(condition)){
  #rowindex <- which(mse_all$samplesize == samplesize[j])
  rowindex <- which(mse_all$condition == condition[j])
  avg_mse <- rbind(avg_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse <- rbind(med_mse, unlist(apply(mse_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
  avg_mse_ratio <- rbind(avg_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) mean(x, na.rm = T))))
  med_mse_ratio <- rbind(med_mse_ratio, unlist(apply(mse_ratio_all[rowindex, ], 2, function(x) median(x, na.rm = T))))
}
print(round(med_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))
output <- rbind(round(med_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) median(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")
print(round(avg_mse_ratio, 4) * 100)
apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))
output <- rbind(round(avg_mse_ratio, 4) * 100, apply(mse_ratio_all, 2, function(x) mean(x, na.rm = TRUE))*100)
write.csv(output, "tempfile.csv")

print(round(med_mse/10000, 3))
apply(mse_all, 2, function(x) median(x, na.rm = TRUE)/10000)
output <- rbind(round(med_mse/10000, 3), apply(mse_all, 2, function(x) median(x, na.rm = TRUE)/10000))
write.csv(output, "tempfile.csv")
print(round(avg_mse/10000, 3))
apply(mse_all, 2, function(x) mean(x, na.rm = TRUE)/10000)
output <- rbind(round(avg_mse/10000, 3), apply(mse_all, 2, function(x) mean(x, na.rm = TRUE)/10000))
write.csv(output, "tempfile.csv")

# coefficients
coef_all <- read.csv("emp_keck_coef.csv", header = F)
head(coef_all)
colnames(coef_all) <- c("condition", "samplesize", "randomdraw","Cont", "SA", "OW", "CovW", "VarW", "CCor", "VanillaS1", "CWM", "SSin", "SSde", "RP", "Top3", "VanillaS2")
samplesize <- unique(coef_all$samplesize)
avg_coef <- c()
for(i in 1:length(samplesize)){
  temp <- coef_all[which(coef_all$samplesize == samplesize[i]), ]
  avg_coef <- rbind(avg_coef, colMeans(temp, na.rm = T))
}
write.csv(avg_coef, "tempfile.csv")
