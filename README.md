# Applying Stacking Techniques on the Wisdom of Crowds Problem

This repository contains R code for implementing algorithms, replicating simulations, tables and figures 
described in the working paper "Combining the Aggregated Judgments: An Efficient Method for Improving Accuracy byStacking Multiple Weighting Models" 
by Shu Huang, Russell Golman, and Stephen B. Broomell.

## How to use
1. Clone this repository to your local computer
2. For each script to run, be sure to change the working directory in each script (searching "setwd()")

## Algorithm
1. The file "WeightingModels.R" contains the replication of previous common approaches to combine judgments including the theoretically optimal weighting method, several crowd selection method, our new parametric crowd selection method and the regularized weighting method with equal-weight prior. 
2. The file "StackingAlgorithm.R" contains the stacking algorithm to combine all individual weighting models. We adopt the LASSO regression as the meta algorithm to select base learners. 
3. The file "simulation_XXX.R" contains the sample code of generating simulated data and running the simulation. 
4. The file "empirical_XXX.R" contains the sample code of analyzing real-world data. 
5. The file "results_analysis.R" contains the sample code of result analysis after obtaining all simulation and empirical analysis results. 
6. The file "correlation analysis.R" contains the sample code of calculating correlation involved in the working paper. 

## Data
The raw data can be obtained by email authors. 

