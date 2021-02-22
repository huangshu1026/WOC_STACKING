# Applying Stacking Techniques on the Wisdom of Crowds Problem

This repository contains R code for implementing algorithms, replicating simulations, tables and figures 
described in the working paper "Combining the Aggregated Judgments: An Efficient Method for Improving Accuracy byStacking Multiple Weighting Models" 
by Shu Huang, Stephen B. Broomell, and Russell Golman. 

## How to use
1. Clone this repository to your local computer
2. For each script to run, be sure to change the working directory in each script (searching "setwd()")

## Algorithm
1. The file "WeightingModels.R" contains the replication of previous common approaches to combine judgments including the theoretically optimal weighting method, several crowd selection method, our new parametric crowd selection method and the regularized weighting method with equal-weight prior. 
2. The file "StackingAlgorithm.R" contains the stacking algorithm to combine all individual weighting models. We adopt the LASSO regression as the meta algorithm to select base learners. 
3. The file "UsefulFunctions" contains other useful functions, such as the MSE calculation and plotting setting functions. 
4. The file "data analysis_spf_frb.R" contains the sample code of analyzing the forecasting data from US SPF. 
5. The file "data analysis_M4.R" contains the sample code of analyzing the M4 point forecasting data. 
6. The file "data analysis_epidemic.R" contains the sample code of analyzing the epidemic forecasting data.

## Data
The raw data can be obtained by email authors. 

