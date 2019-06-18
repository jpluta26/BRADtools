# function to generate ROC curve statistics
# estimates out of sample error with leave-one-out cross validation

getROCstats <- function( x, y, covars = NULL )
# input: x (vector), the independent variable
#        y (binary vector), the dependent variable
#        covars (vector, data.frame, or matrix), other independent variable(s)
# output: a list containing TPR and FPR values; AUC; and number of cases and controls
{
  library(pROC)
  library(ROCR)
  
  # dependent variable must be an integer of values 0 and 1
  if( class(y) != "integer" )
  {
    print("attempting to coerce y to integer...")
    y <- as.integer(y) - 1
   
  }
  
  if(any(!(y %in% c(0,1))))
  {
    stop("y must be an integer with values 0 and 1")
  }
  
  if( !is.null(covars))
  {
    if( !( class(covars) %in% c("data.frame", "matrix")))
    {
      covars <- as.data.frame(covars)
    }
  }
  
  n <- length(x)
  pprob.vec <- rep(0, n)
  
  # LOOCV
  for( i in 1:n )
  {
    # multivariate case
    if( !is.null(covars))
    {
      out <- glm(  y[-i] ~ x[-i] + covars[-i,], family = "binomial")
      B <- out$coefficients
      pprob.vec[i] <- unlist(c(1, x[i], covars[i,])) %*% B
      
    } else
    
    # univariate
    {
      out <- glm( y[-i] ~ x[-i], family = "binomial")
      B <- out$coefficients
      pprob.vec[i] <- c(1, x[i]) %*% B
    }
   
    
  }
  
  pred <- prediction(pprob.vec, y)
  perf <- performance(pred, 'tpr', 'fpr')
  AUC <- performance(pred, 'auc')@y.values
  roc.out <- roc(response = y, predictor = pprob.vec)
  
  # perf: x and y values used in plotting the ROC curve
  # AUC: the AUC value in a single slot
  # roc.out: contains number of cases and controls
  return( list(perf, AUC, roc.out) )
}
