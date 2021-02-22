##############################################################
#### useful functions to calculate the mean squared error ####
##############################################################
###### Calculate MSE ######
CalMSE <- function(weights, X, y){
  # X is a n * M matrix
  # weights is a M * K matrix where K indicates the number of models
  # y is a n*1 vector 
  mse <- colMeans((X %*% weights - y)^2)
  return(mse)
}

CalRMSE <- function(weights, X, y){
  # X is a n * M matrix
  # weights is a M * K matrix where K indicates the number of models
  # y is a n*1 vector 
  rmse <- sqrt(colMeans((X %*% weights - y)^2))
  return(rmse)
}

###### Data Visualization ######
makeTransparent <- function(..., alpha) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}


