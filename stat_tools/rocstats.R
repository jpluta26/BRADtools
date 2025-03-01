
# pluta 5/3/21

# set of functions for ROC curve analysis
# functions to compute ROC curves/AUC for a given data set with various methods of cross-validation
# statistically compare ROC curves
# plotting functionality

library(ggplot2)
library(ROCR)

# make this a basic plot that can be added to
# ------------------------------------------------------------------ #
makeROCPlot <- function( res, plot.title, grp.names )
  # maybe turn this into a base function, and seperate function can add on features
  # like annotation of p-values and AUC?
  #
  # function to plot AUC curve from ROC analysis
  # input:
  #   res (list), results from getROCStats
  #   plot.title (string), title on the plot
  # output:
  #   p1 (ggplot), returned plot of ROC
{
  if( all(sapply(res, class) == "list"))
  {
    dat <- c()
    
    auc <- rep("", length(res))
    
    for(i  in 1:length(res))
    {
      n <- length(unlist(res[[i]][[1]]@x.values))
      tmp <- data.frame(x = rep(0,n),  y= rep(0,n), grp = rep("na",n))
      tmp$x <- unlist(res[[i]][[1]]@x.values)
      tmp$y <- unlist(res[[i]][[1]]@y.values)
      tmp$grp <- rep(i, dim(tmp)[1])
      auc[i] <- paste0(grp.names[i], ": ", round(res[[i]][[2]]$auc,3))
      dat <- rbind(tmp,dat)
    }
    
    colnames(dat) <- c("x", "y", "grp")
    
    p1 <- ggplot(data = dat, aes(x=x,y=y)) +
      geom_line(aes(color=as.factor(grp))) +
      geom_abline(slope=1, intercept =0, color = "red", linetype = "dashed") +
      theme_minimal() 
    
    for(i in 1:length(auc))
    {
      p1 <- p1 + geom_text(x=0.7, y = 0.45 - (0.1 * (i - 1)), label = auc[i])
    }
    
    p1 <- p1 + labs(title = plot.title)
    
    return(p1)
  } 
  
  dat  <- data.frame(x = res[[1]]@x.values, y =  res[[1]]@y.values)
  colnames(df) <- c("x","y")
  
  p1 <- ggplot(data  = df, aes(x = x, y= y)) +
    geom_line(color  = "blue") +
    geom_abline(slope=1, intercept =0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = plot.title,
         subtitle = paste0(" AUC = ", round(res[[2]]$auc,2)))
  return(p1)
}
# ------------------------------------------------------------------ #


# ------------------- makeBaselineROC ------------------------------ #
makeBaselineROC <- function( n )
# create simulated ROC data with AUC of exactly 0.50, currently only 
# works for data with even n. used to test real data against the null of 0.5
# input:
#   n (integer), number of subjects in the real data
# output:
#   roc.out (roc object), simulated ROC data with AUC of 0.5
{
  x <- c(rep(0, n / 2), rep(1, n / 2))
  y <- rep(0,n)
  
  # sample 1/2 of the subjects and match the values
  ind <- sample.int( n, n / 2)
  y[ind] <- x[ind]
  
  # for the other half, make sure the value DONT match
  for( i in seq(1,n,by=1)[-ind] )
  {
    if( y[i] == x[i] )
    {
      if( y[i] == 1 )
      { y[i] <- 0 } else
      { y[i] <- 1 }
    }
  }
  
  # generate the ROC object
  fit <- glm(y ~ x, family = "binomial")
  B <- fit$coefficients
  pprob.vec <- cbind(1,x)  %*% B
  roc.out <- roc(response = y, predictor = as.numeric(pprob.vec))
  return(roc.out)
}
# ------------------------------------------------------------------ #


# ------------------------------------------------------------------ #
# recode a categorical variable to dim(x)[1] - 1 dummies
# this is usually done automatically in glm, but it has to be done explicitly here
# to keep all covariate data in a single general matrix
recodeDummy <- function( x )
{
  dummies <- model.matrix(~as.factor(x))
  
  # drop the intercept
  dummies <- dummies[,2:dim(dummies)[2]]
  colnames(dummies) <- sub("as.factor\\(x\\)", "", colnames(dummies))
  return(dummies)
}
# ------------------------------------------------------------------ #


# ------------------------------------------------------------------ #
# convert factors to dummy columns and insert into data
recodeDummies <- function( dat )
{
  ind <- which(sapply(dat, class) == "factor")
  for( j in ind )
  {
    x <- recodeDummy( dat[,j])
    dat <- dat[,-j]
    dat <- cbind(dat, x)
  }
  
  return(dat)
  
}
# ------------------------------------------------------------------ #



# --------------------------- LOOCV ----------------------- #
# implement leave-one-out cross-validation
LOOCV <- function(X,y)
# input: X, matrix of independent variables
#        y, binary vector of outcome values
# 
# output: pprob.vec (numeric), a vector of the calculated probabilities
{
  n <- dim(X)[1]
  pprob.vec <- rep(0, n)
  
  for( i in 1:n )
  {
    fit <- glm( y[-i] ~ X[-i,], family = "binomial")
    B <- fit$coefficients
    
    pprob.vec[i] <- c(1,X[i,]) %*% B
  }
  
  return(pprob.vec)
}
# ------------------------------------------------------------------ #

# -------------------- getROC --------------------------------------- #
# function to getROC curve statistics and data necessary for plotting
getROC <- function( X, y, cross.validation = NULL )
  # input: x (vector), the independent variable
  #        y (binary vector), the dependent variable
  #        cross.validation (string), one of "LOOCV", "split", or NULL
  #           LOOCV implements leave-one-out cross validation
  #           split divides the data evenly and uses a training and testing set
  # 
  # output: list of perf (roc object containing data points for the ROC curve)
  #   and roc.out object (contains sample size and AUC)
{
  
  
  # dependent variable must be an integer of values 0 and 1
  if( class(y) == "numeric" )
  {
    y <- as.integer(y)
  } else 
    if( class(y) != "integer" )
    {
      print("attempting to coerce y to integer...")
      y <- as.integer(y) - 1
    }
  
  if(any(is.na(y)))
  {
    
    X <- X[!is.na(y),]
    y <- y[!is.na(y)]
    print("warning! NA detected in y. dropping these subjects")
  }
  
  if(any(!(y %in% c(0,1))))
  {
    stop("y must be an integer with values 0 and 1")
  }
  
  # if there are any factors in the data, convert them to dummy variables
  # need the fully specified design matrix to calculate predicted probabilites
  if(any(sapply(X, class) == "factor"))
  {
    X <- recodeDummies(X)
  }
  
  X <- as.matrix(X)
  
  # select validation method
  if( is.null(cross.validation))
  # no validation
  {
    fit <- glm( y ~ X, family = "binomial" )
    B <- fit$coefficients
    pprob.vec <- cbind(1,X) %*% B
  } else if(cross.validation == "LOOCV")
    # leave-one-out cross validation
    {
     
      pprob.vec <- LOOCV(X,y)
      
    } else if(cross.validation == "split")
    # split 
      {   
        # split data evenly into training and testing sets
        n <- dim(X)[1]
        ind <- sample.int( n, n / 2)
        y.train <- y[ind]
        X.train <- X[ind,]
        y.test <- y[-ind]
        X.test <- X[-ind,]
        
        # fit training model and apply to design matrix of test data
        fit <- glm( y.train ~ X.train, family = "binomial")
        B <- fit$coefficients
        pprob.vec <- cbind(1, X.test) %*% B
        y <- y.test
        
      } else
        {
          stop(paste0(cross.validation, " is not a valid choice for cross validation"))
        }
  
  pred <- prediction(pprob.vec, y)
  perf <- performance(pred, 'tpr', 'fpr')
  roc.out <- roc(response = y, predictor = as.numeric(pprob.vec))
  
  # perf: x and y values used in plotting the ROC curve
  # roc.out: contains number of cases and controls, AUC
  return( list(perf, roc.out) )
}
# ------------------------------------------------------------------ #