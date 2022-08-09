setwd()
################################################################################
######## get the correlation matrix among different weighting models ###########
################################################################################
##### M4 data #####
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
rmse_all <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  temp <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_rmse.csv", sep = ""), header = F)
  rmse_all <- rbind(rmse_all, temp)
}
colnames(rmse_all) <- c("file", "item", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_m4 <- rmse_all^2
cordata <- mse_m4[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW", "Stacking")]
cordata_m4 <- cordata
##### SPF US data #####
rmse_all <- read.csv("emp_spfus_rmse.csv", header = F)
colnames(rmse_all) <- c("item", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_spfus <- rmse_all^2
cordata <- mse_spfus[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW", "Stacking")]
cordata_spfus <- cordata
##### flu forecasting #####
rmse_all <- read.csv("emp_flu_rmse.csv", header = F)
colnames(rmse_all) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_flu <- rmse_all^2
cordata <- mse_flu[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW", "Stacking")]
cordata_flu <- cordata
##### covid 19 cumulative death forecasts #####
rmse_all <- read.csv("emp_coviddeath_rmse_new.csv", header = F)
colnames(rmse_all) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_coviddeath <- rmse_all^2
cordata <- mse_coviddeath[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW", "Stacking")]
cordata_coviddeath <- cordata
##### covid 19 new cases forecasts #####
rmse_case1 <- read.csv("emp_covidcase1_rmse.csv", header = F)
rmse_case2 <- read.csv("emp_covidcase2_rmse.csv", header = F)
rmse_case3 <- read.csv("emp_covidcase3_rmse.csv", header = F)
rmse_all <- rbind(rmse_case1, rmse_case2, rmse_case3)
colnames(rmse_all) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_covidcase <- rmse_all^2
cordata <- mse_covidcase[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW", "Stacking")]
cordata_covidcase <- cordata
##### keck and tang (2020) #####
rmse_all <- read.csv("emp_keck_rmse.csv", header = F)
colnames(rmse_all) <- c("condition", "samplesize", "randomdraw", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_keck <- rmse_all^2
cordata <- mse_keck[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW", "Stacking")]
cordata_keck <- cordata
##### combine all data #####
cordata_all <- rbind(cordata_m4, cordata_spfus, cordata_flu, cordata_coviddeath, cordata_covidcase, cordata_keck)


##### pearson correlation #####
cormat <- cor(cordata_all, use = "na.or.complete")
library(corrplot)
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/cor_all.eps", 
           width = 6, height = 5.5)
corrplot(cormat, method="color", col=col(200),  
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         number.cex = 0.8, 
         diag=FALSE 
)
dev.off()
################################################################################
######## get the Spearman correlation between MSE and meta coefficient #########
################################################################################
cor_coef <- c()
cor_pvalue <- c()
##### M4 data #####
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
msedata <- c()
coefdata <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  rmse_sub <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_rmse.csv", sep = ""), header = F)
  colnames(rmse_sub) <- c("file", "item", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
  mse_sub <- rmse_sub^2
  msedata <- rbind(msedata, mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")])
  coef_sub <- read.csv(paste("emp_m4_", item.name[i], time.name[i], "_coef.csv", sep = ""), header = F)
  colnames(coef_sub) <- c("file", "item", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
  coefdata <- rbind(coefdata, coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")])
}
cor_coef_each <- c()
cor_pvalue_each <- c()
for(j in 1:ncol(msedata)){
  testmodel <- cor.test(msedata[ , j], coefdata[ , j], method = "spearman")
  cor_coef_each[j] <- testmodel$estimate
  cor_pvalue_each[j] <- testmodel$p.value
}
cor_coef <- rbind(cor_coef, cor_coef_each)
cor_pvalue <- rbind(cor_pvalue, cor_pvalue_each)
##### SPF US data #####
rmse_sub <- read.csv("emp_spfus_rmse.csv", header = F)
colnames(rmse_sub) <- c("item", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_spfus_coef.csv", header = F)
colnames(coef_sub) <- c("item", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_each <- c()
cor_pvalue_each <- c()
for(j in 1:ncol(msedata)){
  testmodel <- cor.test(msedata[ , j], coefdata[ , j], method = "spearman")
  cor_coef_each[j] <- testmodel$estimate
  cor_pvalue_each[j] <- testmodel$p.value
}
cor_coef <- rbind(cor_coef, cor_coef_each)
cor_pvalue <- rbind(cor_pvalue, cor_pvalue_each)
##### Flu data #####
rmse_sub <- read.csv("emp_flu_rmse.csv", header = F)
colnames(rmse_sub) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_flu_coef.csv", header = F)
coefs <- coef_sub$V3
coefs <- matrix(coefs, nrow = 44, ncol = 13, byrow = TRUE)
coef_sub <- cbind(rep(1:4, each = 11), rep(1:11, 4), coefs)
coef_sub <- as.data.frame(coef_sub)
colnames(coef_sub) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_each <- c()
cor_pvalue_each <- c()
for(j in 1:ncol(msedata)){
  testmodel <- cor.test(msedata[ , j], coefdata[ , j], method = "spearman")
  cor_coef_each[j] <- testmodel$estimate
  cor_pvalue_each[j] <- testmodel$p.value
}
cor_coef <- rbind(cor_coef, cor_coef_each)
cor_pvalue <- rbind(cor_pvalue, cor_pvalue_each)
##### covid death data #####
rmse_sub <- read.csv("emp_coviddeath_rmse_new.csv", header = F)
colnames(rmse_sub) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_coviddeath_coef_new.csv", header = F)
head(coef_sub)
colnames(coef_sub) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_each <- c()
cor_pvalue_each <- c()
for(j in 1:ncol(msedata)){
  testmodel <- cor.test(msedata[ , j], coefdata[ , j], method = "spearman")
  cor_coef_each[j] <- testmodel$estimate
  cor_pvalue_each[j] <- testmodel$p.value
}
cor_coef <- rbind(cor_coef, cor_coef_each)
cor_pvalue <- rbind(cor_pvalue, cor_pvalue_each)
##### covid case data #####
rmse_case1 <- read.csv("emp_covidcase1_rmse.csv", header = F)
rmse_case2 <- read.csv("emp_covidcase2_rmse.csv", header = F)
rmse_case3 <- read.csv("emp_covidcase3_rmse.csv", header = F)
rmse_sub <- rbind(rmse_case1, rmse_case2, rmse_case3)
colnames(rmse_sub) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_case1 <- read.csv("emp_covidcase1_coef.csv", header = F)
coef_case2 <- read.csv("emp_covidcase2_coef.csv", header = F)
coef_case3 <- read.csv("emp_covidcase3_coef.csv", header = F)
coef_sub <- rbind(coef_case1, coef_case2, coef_case3)
head(coef_sub)
colnames(coef_sub) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_each <- c()
cor_pvalue_each <- c()
for(j in 1:ncol(msedata)){
  testmodel <- cor.test(msedata[ , j], coefdata[ , j], method = "spearman")
  cor_coef_each[j] <- testmodel$estimate
  cor_pvalue_each[j] <- testmodel$p.value
}
cor_coef <- rbind(cor_coef, cor_coef_each)
cor_pvalue <- rbind(cor_pvalue, cor_pvalue_each)
##### keck and tang (2020) #####
rmse_sub <- read.csv("emp_keck_rmse.csv", header = F)
colnames(rmse_sub) <- c("condition", "samplesize", "randomdraw", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_keck_coef.csv", header = F)
head(coef_sub)
colnames(coef_sub) <- c("condition", "samplesize", "randomdraw", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_each <- c()
cor_pvalue_each <- c()
for(j in 1:ncol(msedata)){
  testmodel <- cor.test(msedata[ , j], coefdata[ , j], method = "spearman")
  cor_coef_each[j] <- testmodel$estimate
  cor_pvalue_each[j] <- testmodel$p.value
}
cor_coef <- rbind(cor_coef, cor_coef_each)
cor_pvalue <- rbind(cor_pvalue, cor_pvalue_each)
##### visualization: bar plot #####
setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/barplot_m4.eps", 
           width = 8, height = 5)
bp <- barplot(cor_coef[1, ], ylim = c(min(cor_coef)-0.2, max(cor_coef)+0.1), 
              ylab = "Correlation Coefficient", xlab = "Weighting Methods")
text(bp, rep(min(cor_coef), length(bp)), srt=45,
     c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW"))
text(bp, cor_coef[1, ]+c(rep(-0.05, 6), rep(0.05, 2), rep(-0.05, 4)), 
     round(cor_coef[1, ], 2))
round(cor_pvalue[1, ], 3)
text(bp, cor_coef[1, ]+c(rep(-0.1, 6), rep(0.1, 2), rep(-0.1, 4)), 
     rep("***", length(bp)))
dev.off()

setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/barplot_spfus.eps", 
           width = 8, height = 5)
bp <- barplot(cor_coef[2, ], ylim = c(min(cor_coef)-0.2, max(cor_coef)+0.1), 
              ylab = "Correlation Coefficient", xlab = "Weighting Methods")
text(bp, rep(min(cor_coef), length(bp)), srt=45,
     c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW"))
text(bp, cor_coef[2, ]+c(rep(-0.05, 2), rep(0.05, 3), -0.05, rep(0.05, 3), rep(-0.05, 3)), 
     round(cor_coef[2, ], 2))
round(cor_pvalue[2, ], 3)
#text(bp, cor_coef[1, ]+c(rep(-0.1, 6), rep(0.1, 2), rep(-0.1, 4)), 
#     rep("***", length(bp)))
dev.off()

setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/barplot_flu.eps", 
           width = 8, height = 5)
bp <- barplot(cor_coef[3, ], ylim = c(min(cor_coef)-0.2, max(cor_coef)+0.1), 
              ylab = "Correlation Coefficient", xlab = "Weighting Methods")
text(bp, rep(min(cor_coef), length(bp)), srt=45,
     c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW"))
text(bp, cor_coef[3, ]+c(rep(0.05, 3), rep(-0.05, 3), 0.05, -0.05, 0.05, rep(-0.05, 3)), 
     round(cor_coef[3, ], 2))
round(cor_pvalue[3, ], 3)
text(bp[1], cor_coef[3, 1]+0.1, 
     rep("**", length(bp)))
dev.off()

setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/barplot_coviddeath.eps", 
           width = 8, height = 5)
bp <- barplot(cor_coef[4, ], ylim = c(min(cor_coef)-0.3, max(cor_coef)+0.1), 
              ylab = "Correlation Coefficient", xlab = "Weighting Methods")
text(bp, rep(min(cor_coef)-0.2, length(bp)), srt=45,
     c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW"))
text(bp, cor_coef[4, ]+c(-0.05, rep(0.05, 8), -0.05, 0.05, -0.05), 
     round(cor_coef[4, ], 2))
round(cor_pvalue[4, ], 3)
text(bp[c(1, 3, 6, 7, 9, 11)], cor_coef[4, c(1, 3, 6, 7, 9, 11)]+c(-0.1, 0.1, 0.1, 0.1, 0.1, 0.1), 
     c("***", "**", "*", "***", "*", "***"))
dev.off()

setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/barplot_covidcase.eps", 
           width = 8, height = 5)
bp <- barplot(cor_coef[5, ], ylim = c(min(cor_coef)-0.2, max(cor_coef)+0.1), 
              ylab = "Correlation Coefficient", xlab = "Weighting Methods")
text(bp, rep(min(cor_coef), length(bp)), srt=45,
     c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW"))
text(bp, cor_coef[5, ]+c(rep(-0.05, 8), rep(0.05, 4)), 
     round(cor_coef[5, ], 2))
round(cor_pvalue[5, ], 3)
text(bp[c(1,2,7,9,10,11,12)], cor_coef[5, c(1,2,7,9,10,11,12)]+c(-0.1, -0.1, -0.1, 0.1, 0.1, 0.1,0.1), 
     c("***", "***", "***", "***", "***", "**", "***"))
dev.off()


setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/barplot_keck.eps", 
           width = 8, height = 5)
bp <- barplot(cor_coef[6, ], ylim = c(min(cor_coef)-0.2, max(cor_coef)+0.1), 
              ylab = "Correlation Coefficient", xlab = "Weighting Methods")
text(bp, rep(min(cor_coef), length(bp)), srt=45,
     c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW"))
text(bp, cor_coef[6, ]+c(-0.05,rep(0.05, 4),-0.05,0.05,-0.05,rep(0.05, 4)), 
     round(cor_coef[6, ], 2))
round(cor_pvalue[6, ], 4)
text(bp[c(2,3,6,9,11,12)], cor_coef[6, c(2,3,6,9,11,12)]+c(0.1,0.1,-0.1,0.1,0.1,0.1), 
     c("**", "**", "**", "***", "***", "**"))
dev.off()

######################################################################################################
######## get the Spearman correlation between MSE and meta coefficient for each target value #########
######################################################################################################
cor_coef_avg <- c()
cor_coef_se <- c()
##### M4 data #####
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
msedata <- c()
coefdata <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  rmse_sub <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_rmse.csv", sep = ""), header = F)
  colnames(rmse_sub) <- c("file", "item", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
  mse_sub <- rmse_sub^2
  msedata <- rbind(msedata, mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")])
  coef_sub <- read.csv(paste("emp_m4_", item.name[i], time.name[i], "_coef.csv", sep = ""), header = F)
  colnames(coef_sub) <- c("file", "item", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
  coefdata <- rbind(coefdata, coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")])
}
cor_coef_m4 <- c()
for(j in 1:nrow(msedata)){
  if(is.na(coefdata[j, 1]) == FALSE){
    #testmodel <- cor.test(unlist(msedata[j, ]), unlist(coefdata[j, ]), method = "spearman", exact=FALSE)
    testmodel <- cor.test(unlist(msedata[j, 2:12]), unlist(coefdata[j, 2:12]), method = "spearman", exact=FALSE)
    cor_coef_m4[j] <- testmodel$estimate
  }
}
cor_coef_avg[1] <- mean(cor_coef_m4, na.rm = TRUE)
cor_coef_se[1] <- t.test(cor_coef_m4)$stderr
##### SPF US data #####
rmse_sub <- read.csv("emp_spfus_rmse.csv", header = F)
colnames(rmse_sub) <- c("item", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_spfus_coef.csv", header = F)
colnames(coef_sub) <- c("item", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_spfus <- c()
for(j in 1:nrow(msedata)){
  if(is.na(coefdata[j, 1]) == FALSE){
    #testmodel <- cor.test(unlist(msedata[j, ]), unlist(coefdata[j, ]), method = "spearman", exact=FALSE)
    testmodel <- cor.test(unlist(msedata[j, 2:12]), unlist(coefdata[j, 2:12]), method = "spearman", exact=FALSE)
    cor_coef_spfus[j] <- testmodel$estimate
  }
}
cor_coef_avg[2] <- mean(cor_coef_spfus, na.rm = TRUE)
cor_coef_se[2] <- t.test(cor_coef_spfus)$stderr

##### Flu data #####
rmse_sub <- read.csv("emp_flu_rmse.csv", header = F)
colnames(rmse_sub) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_flu_coef.csv", header = F)
coefs <- coef_sub$V3
coefs <- matrix(coefs, nrow = 44, ncol = 13, byrow = TRUE)
coef_sub <- cbind(rep(1:4, each = 11), rep(1:11, 4), coefs)
coef_sub <- as.data.frame(coef_sub)
colnames(coef_sub) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_flu <- c()
for(j in 1:nrow(msedata)){
  if(is.na(coefdata[j, 1]) == FALSE){
    #testmodel <- cor.test(unlist(msedata[j, ]), unlist(coefdata[j, ]), method = "spearman", exact=FALSE)
    testmodel <- cor.test(unlist(msedata[j, 2:12]), unlist(coefdata[j, 2:12]), method = "spearman", exact=FALSE)
    cor_coef_flu[j] <- testmodel$estimate
  }
}
cor_coef_avg[3] <- mean(cor_coef_flu, na.rm = TRUE)
cor_coef_se[3] <- t.test(cor_coef_flu)$stderr

##### covid death data #####
rmse_sub <- read.csv("emp_coviddeath_rmse_new.csv", header = F)
colnames(rmse_sub) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_coviddeath_coef_new.csv", header = F)
head(coef_sub)
colnames(coef_sub) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_death <- c()
for(j in 1:nrow(msedata)){
  if(is.na(coefdata[j, 1]) == FALSE){
    #testmodel <- cor.test(unlist(msedata[j, ]), unlist(coefdata[j, ]), method = "spearman", exact=FALSE)
    testmodel <- cor.test(unlist(msedata[j, 2:12]), unlist(coefdata[j, 2:12]), method = "spearman", exact=FALSE)
    cor_coef_death[j] <- testmodel$estimate
  }
}
cor_coef_avg[4] <- mean(cor_coef_death, na.rm = TRUE)
cor_coef_se[4] <- t.test(cor_coef_death)$stderr
##### covid case data #####
rmse_case1 <- read.csv("emp_covidcase1_rmse.csv", header = F)
rmse_case2 <- read.csv("emp_covidcase2_rmse.csv", header = F)
rmse_case3 <- read.csv("emp_covidcase3_rmse.csv", header = F)
rmse_sub <- rbind(rmse_case1, rmse_case2, rmse_case3)
colnames(rmse_sub) <- c("Target", "Location", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_case1 <- read.csv("emp_covidcase1_coef.csv", header = F)
coef_case2 <- read.csv("emp_covidcase2_coef.csv", header = F)
coef_case3 <- read.csv("emp_covidcase3_coef.csv", header = F)
coef_sub <- rbind(coef_case1, coef_case2, coef_case3)
head(coef_sub)
colnames(coef_sub) <- c("Target", "Location", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_case <- c()
for(j in 1:nrow(msedata)){
  if(is.na(coefdata[j, 1]) == FALSE){
    #testmodel <- cor.test(unlist(msedata[j, ]), unlist(coefdata[j, ]), method = "spearman", exact=FALSE)
    testmodel <- cor.test(unlist(msedata[j, 2:12]), unlist(coefdata[j, 2:12]), method = "spearman", exact=FALSE)
    cor_coef_case[j] <- testmodel$estimate
  }
}
cor_coef_avg[5] <- mean(cor_coef_case, na.rm = TRUE)
cor_coef_se[5] <- t.test(cor_coef_case)$stderr
##### keck and tang (2020) #####
rmse_sub <- read.csv("emp_keck_rmse.csv", header = F)
colnames(rmse_sub) <- c("condition", "samplesize", "randomdraw", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos", "Stacking")
mse_sub <- rmse_sub^2
msedata <- mse_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
coef_sub <- read.csv("emp_keck_coef.csv", header = F)
head(coef_sub)
colnames(coef_sub) <- c("condition", "samplesize", "randomdraw", "Cont", "SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
coefdata <- coef_sub[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
cor_coef_keck <- c()
for(j in 1:nrow(msedata)){
  if(is.na(coefdata[j, 1]) == FALSE){
    #testmodel <- cor.test(unlist(msedata[j, ]), unlist(coefdata[j, ]), method = "spearman", exact=FALSE)
    testmodel <- cor.test(unlist(msedata[j, 2:12]), unlist(coefdata[j, 2:12]), method = "spearman", exact=FALSE)
    cor_coef_keck[j] <- testmodel$estimate
  }
}
cor_coef_avg[6] <- mean(cor_coef_keck, na.rm = TRUE)
cor_coef_se[6] <- t.test(cor_coef_keck)$stderr
##### visualization: bar plot #####
setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/Spreaman_Corr_Avg.eps", 
           width = 6, height = 5)
bp <- barplot(cor_coef_avg, ylim = c(-0.5, 0.5),
              ylab = "Average Correlation Coefficient", xlab = "Data")
text(bp, rep(-0.35, length(bp)), srt = 30, 
     c("M4", "SPF", "Flu", "COVID-19 Death", "COVID-19 Case", "Keck&Tang"))
points(bp, cor_coef_avg, pch = 16, cex = 0.5)
points(bp, cor_coef_avg+1.96*cor_coef_se, pch = "-", cex = 2)
points(bp, cor_coef_avg-1.96*cor_coef_se, pch = "-", cex = 2)
segments(bp, cor_coef_avg+1.96*cor_coef_se, bp, cor_coef_avg-1.96*cor_coef_se, lwd = 2)
text(bp, rep(0.35, length(bp)), round(cor_coef_avg, 2))
abline(h = 0)
dev.off()

##################################################################################################
######## get the correlation matrix among different weighting models using predictions ###########
##################################################################################################
##### M4 data #####
item.name <- c(rep("demo", 5), rep("finance", 5), rep("industry", 5), 
               rep("macro", 5), rep("micro", 5), rep("other", 6))
time.name <- c(rep(c("daily", "monthly", "quarterly", "weekly", "yearly"), 5), 
               c("daily", "hourly", "monthly", "quarterly", "weekly", "yearly"))
corr_all <- c()
for(i in c(1:4, 6:9, 11:14, 16:19, 21:24, 26:30)){
  temp <- read.csv(paste("emp_m4_", item.name[i], time.name[i],"_corr.csv", sep = ""), header = F)
  corr_all <- rbind(corr_all, temp)
}
corr_array <- array(NA, dim = c(12, 12, nrow(corr_all)))
for(i in 1:nrow(corr_all)){
  temp <- unlist(corr_all[i, 3:ncol(corr_all)])
  temp_mat <- matrix(NA, nrow = 12, ncol = 12)
  temp_mat[upper.tri(temp_mat)] <- temp
  temp_mat[lower.tri(temp_mat)] <- t(temp_mat)[lower.tri(t(temp_mat))]
  diag(temp_mat) <- 1
  corr_array[ , , i] <- temp_mat
}
corr_avg_m4 <- apply(corr_array, c(1, 2), function(x) mean(x, na.rm = TRUE))
colnames(corr_avg_m4) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
rownames(corr_avg_m4) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")

##### SPF US data #####
corr_all <- read.csv("emp_spfus_corr.csv", header = F)
head(corr_all)
corr_array <- array(NA, dim = c(12, 12, nrow(corr_all)))
for(i in 1:nrow(corr_all)){
  temp <- unlist(corr_all[i, 2:ncol(corr_all)])
  temp_mat <- matrix(NA, nrow = 12, ncol = 12)
  temp_mat[upper.tri(temp_mat)] <- temp
  temp_mat[lower.tri(temp_mat)] <- t(temp_mat)[lower.tri(t(temp_mat))]
  diag(temp_mat) <- 1
  corr_array[ , , i] <- temp_mat
}
corr_avg_spf <- apply(corr_array, c(1, 2), function(x) mean(x, na.rm = TRUE))
colnames(corr_avg_spf) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
rownames(corr_avg_spf) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
##### flu forecasting #####
corr_all <- read.csv("emp_flu_corr.csv", header = F)
head(corr_all)
corr_array <- array(NA, dim = c(12, 12, nrow(corr_all)))
for(i in 1:nrow(corr_all)){
  temp <- unlist(corr_all[i, 3:ncol(corr_all)])
  temp_mat <- matrix(NA, nrow = 12, ncol = 12)
  temp_mat[upper.tri(temp_mat)] <- temp
  temp_mat[lower.tri(temp_mat)] <- t(temp_mat)[lower.tri(t(temp_mat))]
  diag(temp_mat) <- 1
  corr_array[ , , i] <- temp_mat
}
corr_avg_flu <- apply(corr_array, c(1, 2), function(x) mean(x, na.rm = TRUE))
colnames(corr_avg_flu) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
rownames(corr_avg_flu) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
##### covid 19 cumulative death forecasts #####
corr_all <- read.csv("emp_coviddeath_corr_new.csv", header = F)
head(corr_all)
corr_array <- array(NA, dim = c(12, 12, nrow(corr_all)))
for(i in 1:nrow(corr_all)){
  temp <- unlist(corr_all[i, 3:ncol(corr_all)])
  temp_mat <- matrix(NA, nrow = 12, ncol = 12)
  temp_mat[upper.tri(temp_mat)] <- temp
  temp_mat[lower.tri(temp_mat)] <- t(temp_mat)[lower.tri(t(temp_mat))]
  diag(temp_mat) <- 1
  corr_array[ , , i] <- temp_mat
}
corr_avg_death <- apply(corr_array, c(1, 2), function(x) mean(x, na.rm = TRUE))
colnames(corr_avg_death) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
rownames(corr_avg_death) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
##### covid 19 new cases forecasts #####
corr_case1 <- read.csv("emp_covidcase1_corr.csv", header = F)
corr_case2 <- read.csv("emp_covidcase2_corr.csv", header = F)
corr_case3 <- read.csv("emp_covidcase3_corr.csv", header = F)
corr_all <- rbind(corr_case1, corr_case2, corr_case3)
head(corr_all)
corr_array <- array(NA, dim = c(12, 12, nrow(corr_all)))
for(i in 1:nrow(corr_all)){
  temp <- unlist(corr_all[i, 3:ncol(corr_all)])
  temp_mat <- matrix(NA, nrow = 12, ncol = 12)
  temp_mat[upper.tri(temp_mat)] <- temp
  temp_mat[lower.tri(temp_mat)] <- t(temp_mat)[lower.tri(t(temp_mat))]
  diag(temp_mat) <- 1
  corr_array[ , , i] <- temp_mat
}
corr_avg_case <- apply(corr_array, c(1, 2), function(x) mean(x, na.rm = TRUE))
colnames(corr_avg_case) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
rownames(corr_avg_case) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
##### keck and tang (2020) #####
corr_all <- read.csv("emp_keck_corr.csv", header = F)
head(corr_all)
corr_array <- array(NA, dim = c(12, 12, nrow(corr_all)))
for(i in 1:nrow(corr_all)){
  temp <- unlist(corr_all[i, 4:ncol(corr_all)])
  temp_mat <- matrix(NA, nrow = 12, ncol = 12)
  temp_mat[upper.tri(temp_mat)] <- temp
  temp_mat[lower.tri(temp_mat)] <- t(temp_mat)[lower.tri(t(temp_mat))]
  diag(temp_mat) <- 1
  corr_array[ , , i] <- temp_mat
}
corr_avg_keck <- apply(corr_array, c(1, 2), function(x) mean(x, na.rm = TRUE))
colnames(corr_avg_keck) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")
rownames(corr_avg_keck) <- c("SA", "OW", "CovW", "VarW", "CCor", "RegSA", "CWM", "SSIN", "SSDE", "RP", "Top3", "RegSAPos")

##### pearson correlation #####
#plotdata <- corr_avg_m4[ , c("SA", "CWM", "SSIN", "SSDE", "RP", "Top3", "VarW", "CCor", "RegSA", "RegSAPos", "OW", "CovW")]
plotdata <- corr_avg_m4
plotdata <- corr_avg_spf
plotdata <- corr_avg_flu
plotdata <- corr_avg_death
plotdata <- corr_avg_case
plotdata <- corr_avg_keck
library(corrplot)
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/correlation matrix/pred_cor_case.eps", 
           width = 6, height = 5.5)
corrplot(plotdata, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         number.cex = 0.8, 
         diag=FALSE 
)
dev.off()
###################### 
##### barplot ######
overMSE <- read.csv("overallMSE.csv", header = T)
avgMetaCoef <- read.csv("avgMetaCoef.csv", header = T)

pearson_corr <- c()
pearson_corr.pvalue <- c()
for(i in 1:nrow(overMSE)){
  x <- as.numeric(unlist(overMSE[i, ]))
  y <- unlist(avgMetaCoef[i, ])
  temp_test <- cor.test(x, y)
  pearson_corr[i] <- temp_test$estimate
  pearson_corr.pvalue[i] <- temp_test$p.value
}
pearson_corr
pearson_corr.pvalue


setEPS()
postscript("~/Documents/WOC_STACKING/science_advances/R1/barplot_corr.eps", 
           width = 6, height = 4.5)
bp <- barplot(pearson_corr, xlab = "Data", ylab = "Pearson Correlation", ylim = c(-0.25, 0.3))
text(bp, rep(-0.15, length(bp)), c("M4", "SPF", "Flu", "COVID19", "COVID19", "Keck and"), cex = 0.8)
text(bp[4:6], rep(-0.18, 3), c("Deaths", "New", "Tang (2020)"), cex = 0.8)
text(bp[5:6], rep(-0.21, 2), c("Cases", "Experiment"), cex = 0.8)
text(bp, pearson_corr+c(-0.03, rep(0.03, 5)), round(pearson_corr, 2))
dev.off()
