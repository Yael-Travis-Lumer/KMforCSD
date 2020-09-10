#' Simulate Current Status Data from Exponential Failure and Censoring Times
#'
#' \code{exp_data} returns a list that includes a data frame of size [n,3]. The first column corresponds to the univarite covarite Z, the second column is the current status indicator, and the third column is the censoring times. Note that when test=TRUE we get an extra column with the true uncensored failure time T.
#'
#' @param n A positive integer (sample size)
#' @param tau A non-negative value (defining the support of T and C). Default value is 3.
#' @param test Logical indicating whether the function should return in addition also the true (uncensored) failure time, together with Bayes risk, for testing. Default value is False and then only current status data is returned (for training).
#' @return A list containing: (1) the simulated data in the format of a data frame with the columns Z, delta, C, and possibly also T (for testing only), and (2) the Bayes risk which can only be evaluated when test=TRUE (otherwise it is NA).
#' @examples
#' d <- exp_data(n=100)
#' @export
exp_data <- function(n,tau=3,test=FALSE){
  Bayes_risk <- NA
  Z <- stats::runif(n,0,1)
  T <- stats::rexp(n,exp(-0.5*Z))
  C <- stats::runif(n,0,tau)
  T[T>tau] <- tau
  C[C>tau] <- tau
  delta <- (T<=C)
  delta <- as.numeric(delta)
  d <- data.frame(cbind(Z,delta,C))
  covar_names <- paste("Z",1, sep = "")
  colnames(d) <- c(covar_names,"delta","C")
  if (isTRUE(test)){
    expect <- exp(0.5*Z)*(1-exp(-tau*exp(-0.5*Z)))
    d <- data.frame(cbind(Z,delta,C,T))
    colnames(d) <- c(covar_names,"delta","C","T")
    Bayes_risk <- mean((T-expect)^2)
  }
  return(list(data=d,Bayes_risk=Bayes_risk))
}

#' Simulate Current Status Data from Triangle Shaped Failure Time and Uniform Censoring Times
#'
#' \code{triangle_data} returns a list that includes a data frame of size [n,3]. The first column corresponds to the univarite covarite Z, the second column is the current status indicator, and the third column is the censoring times. Note that when test=TRUE we get an extra column with the true uncensored failure time T.
#'
#' @param n A positive integer (sample size)
#' @param tau A non-negative value (defining the support of T and C). Default value is 8.
#' @param test Logical indicating whether the function should return in addition also the true (uncensored) failure time, together with Bayes risk, for testing. Default value is False and then only current status data is returned (for training).
#' @return A list containing: (1) the simulated data in the format of a data frame with the columns Z, delta, C, and possibly also T (for testing only), and (2) the Bayes risk which can only be evaluated when test=TRUE (otherwise it is NA).
#' @examples
#' d <- triangle_data(n=100)
#' @export
triangle_data <- function(n,tau=8,test=FALSE){
  Bayes_risk <- NA
  Z <- stats::runif(n,0,1)
  C <- stats::runif(n,0,tau)
  epsilon <- stats::rnorm(n)
  T <- (4+6*Z)*(Z<0.5)+(10-6*Z)*(Z>0.5)+epsilon
  T[T>tau] <- tau
  C[C>tau] <- tau
  delta <- (T<=C)
  delta <- as.numeric(delta)
  d <- data.frame(cbind(Z,delta,C))
  covar_names <- paste("Z",1, sep = "")
  colnames(d) <- c(covar_names,"delta","C")
  if (isTRUE(test)){
    expect <- ((4+6*Z-tau)*stats::pnorm(tau-4-6*Z)-((2*pi)^-0.5)*exp(-0.5*(tau-4-6*Z)^2)+tau)*(Z<0.5)+((10-6*Z-tau)*stats::pnorm(tau-10+6*Z)-((2*pi)^-0.5)*exp(-0.5*(tau-10+6*Z)^2)+tau)*(Z>0.5) #should this be returned? maybe as a function for plotting?
    d <- data.frame(cbind(Z,delta,C,T))
    colnames(d) <- c(covar_names,"delta","C","T")
    Bayes_risk <- mean((T-expect)^2)
  }
  return(list(data=d,Bayes_risk=Bayes_risk))
}

#' Simulate Current Status Data from Weibull Failure Time and Uniform Censoring Times
#'
#' \code{weibull_data} returns a list that includes a data frame of size [n,3]. The first column corresponds to the univarite covarite Z, the second column is the current status indicator, and the third column is the censoring times. Note that when test=TRUE we get an extra column with the true uncensored failure time T.
#'
#' @param n A positive integer (sample size)
#' @param tau A non-negative value (defining the support of T and C). Default value is 1.
#' @param test Logical indicating whether the function should return in addition also the true (uncensored) failure time, together with Bayes risk, for testing. Default value is False and then only current status data is returned (for training).
#' @return A list containing: (1) the simulated data in the format of a data frame with the columns Z, delta, C, and possibly also T (for testing only), and (2) the Bayes risk which can only be evaluated when test=TRUE (otherwise it is NA).
#' @examples
#' d <- weibull_data(n=100)
#' @export
weibull_data <- function(n,tau=1,test=FALSE){
  Bayes_risk <- NA
  Z <- stats::runif(n,0,1)
  scale <- exp(-0.5*Z)
  shape <- 2
  T <- stats::rweibull(n,shape,scale)
  C <- stats::runif(n,0,tau)
  T[T>tau] <- tau
  C[C>tau] <- tau
  delta <- (T<=C)
  delta <- as.numeric(delta)
  d <- data.frame(cbind(Z,delta,C))
  covar_names <- paste("Z",1, sep = "")
  colnames(d) <- c(covar_names,"delta","C")
  if (isTRUE(test)){
    expect <- scale* stats::pgamma((tau/scale)^shape,1+1/shape)*gamma(1+1/shape)+tau*exp(-(tau/scale)^shape)
    d <- data.frame(cbind(Z,delta,C,T))
    colnames(d) <- c(covar_names,"delta","C","T")
    Bayes_risk <- mean((T-expect)^2)
  }
  return(list(data=d,Bayes_risk=Bayes_risk))
}


#' Simulate Current Status Data from Weibull Failure Time Conditional on Multivariate Covariates
#'
#' \code{multi_weibull_data} returns a list that includes a data frame of size [n,p+2]. The first p columns corresponds to the multivarite covarites Z, the p+1 column is the current status indicator, and the p+2 column is the censoring times. Note that when test=TRUE we get an extra column with the true uncensored failure time T.
#'
#' @param n A positive integer (sample size)
#' @param p a positive integer indictating the dimension of the covariates. Default value is 10.
#' @param tau A non-negative value (defining the support of T and C). Default value is 2.
#' @param test Logical indicating whether the function should return in addition also the true (uncensored) failure time, together with Bayes risk, for testing. Default value is False and then only current status data is returned (for training).
#' @return A list containing: (1) the simulated data in the format of a data frame with the columns of  the multivariate covariates Z, delta, C, and possibly also T (for testing only), and (2) the Bayes risk which can only be evaluated when test=TRUE (otherwise it is NA).
#' @examples
#' d <- multi_weibull_data(n=100)
#' @export
multi_weibull_data <- function(n,p=10,tau=2,test=FALSE){
  Bayes_risk <- NA
  beta <- matrix(c(-0.5,2,-1,matrix(0,p-3,1)),p,1)
  Z <- matrix(stats::runif(n*p,0,1),n,p)
  scale <- exp(Z%*%beta)
  shape <- 2
  T <- stats::rweibull(n,shape,scale)
  C <- stats::runif(n,0,tau)
  T[T>tau] <- tau
  C[C>tau] <- tau
  delta <- (T<=C)
  delta <- as.numeric(delta)
  d <- data.frame(cbind(Z,delta,C))
  covar_names <- paste("Z",1:p, sep = "")
  colnames(d) <- c(covar_names,"delta","C")
  if (isTRUE(test)){
    expect <- exp(Z%*%beta)* stats::pgamma((tau/exp(Z%*%beta))^shape,(1+1/shape)*matrix(1,n,1))*gamma((1+1/shape)*matrix(1,n,1))+tau*exp(-(tau/exp(Z%*%beta))^shape)
    d <- data.frame(cbind(Z,delta,C,T))
    colnames(d) <- c(covar_names,"delta","C","T")
    Bayes_risk <- mean((T-expect)^2)
  }
  return(list(data=d,Bayes_risk=Bayes_risk))
}


#' Simulate Current Status Data from Log-Normal Failure Time Conditional on Multivariate Covariates
#'
#' \code{multi_LN_data} returns a list that includes a data frame of size [n,p+2]. The first p columns corresponds to the multivarite covarites Z, the p+1 column is the current status indicator, and the p+2 column is the censoring times. Note that when test=TRUE we get an extra column with the true uncensored failure time T.
#' @param n A positive integer (sample size)
#' @param p a positive integer indictating the dimension of the covariates. Default value is 10.
#' @param tau A non-negative value (defining the support of T and C). Default value is 7.
#' @param test Logical indicating whether the function should return in addition also the true (uncensored) failure time, together with Bayes risk, for testing. Default value is False and then only current status data is returned (for training).
#' @param sigma A numeric value indicating the standard deviation of the log-normal distribution.
#' @return A list containing: (1) the simulated data in the format of a data frame with the columns of  the multivariate covariates Z, delta, C, and possibly also T (for testing only), and (2) the Bayes risk which can only be evaluated when test=TRUE (otherwise it is NA).
#' @examples
#' d <- multi_LN_data(n=100)
#' @export
multi_LN_data <- function(n,p=10,tau=7,test=FALSE,sigma=1){
  Bayes_risk <- NA
  beta <- matrix(c(0.3,0.5,0.2,matrix(0,p-3,1)),p,1)
  Z <- matrix(stats::runif(n*p,0,1),n,p)
  x <- (Z%*%beta)
  T <- stats::rlnorm(n,meanlog=0.5*x,sdlog=sigma)
  C <- stats::runif(n,0,tau)
  T[T>tau] <- tau
  C[C>tau] <- tau
  delta <- (T<=C)
  delta <- as.numeric(delta)
  d <- data.frame(cbind(Z,delta,C))
  covar_names <- paste("Z",1:p, sep = "")
  colnames(d) <- c(covar_names,"delta","C")
  if (isTRUE(test)){
    expect <- stats::pnorm(q=log(tau),mean=0.5*x+sigma^2,sd=sigma)*exp(0.5*x+0.5*sigma^2)+tau*(1-stats::pnorm(q=(log(tau)-0.5*x)/sigma))
    d <- data.frame(cbind(Z,delta,C,T))
    colnames(d) <- c(covar_names,"delta","C","T")
    Bayes_risk <- mean((T-expect)^2)
  }
  return(list(data=d,Bayes_risk=Bayes_risk))
}


#' Simulate Current Status Data from Half Circle Shaped Failure Time and Uniform Censoring Times
#'
#' \code{half_circle_MVdata} returns a list that includes a data frame of size [n,p+2]. The first p columns corresponds to the multivarite covarites Z, the p+1 column is the current status indicator, and the p+2 column is the censoring times. Note that when test=TRUE we get an extra column with the true uncensored failure time T.
#' @param n A positive integer (sample size)
#' @param p a positive integer indicating the dimension of the covariates. Default value is 2.
#' @param tau A non-negative value (defining the support of T and C). Default value is 9.
#' @param test Logical indicating whether the function should return in addition also the true (uncensored) failure time, together with Bayes risk, for testing. Default value is False and then only current status data is returned (for training).
#' @return A list containing: (1) the simulated data in the format of a data frame with the columns Z, delta, C, and possibly also T (for testing only), and (2) the Bayes risk which can only be evaluated when test=TRUE (otherwise it is NA).
#' @examples
#' d <- half_circle_MVdata(n=100)
#' @export
half_circle_MVdata <- function(n,p=2,tau=9,test=FALSE){
  Bayes_risk <- NA
  Z <- matrix(stats::runif(n*p,0,1),n,p)
  beta <- matrix(1,p,1)
  C <- stats::runif(n,0,tau)
  epsilon <- stats::rnorm(n)
  T <- 1+8*(Z%*%beta-1)^2+epsilon
  T[T>tau] <- tau
  C[C>tau] <- tau
  delta <- (T<=C)
  delta <- as.numeric(delta)
  d <- data.frame(cbind(Z,delta,C))
  covar_names <- paste("Z",1:p, sep = "")
  colnames(d) <- c(covar_names,"delta","C")
  if (isTRUE(test)){
    expect <-  1+8*(Z%*%beta-1)^2
    d <- data.frame(cbind(Z,delta,C,T))
    colnames(d) <- c(covar_names,"delta","C","T")
    Bayes_risk <- mean((T-expect)^2)
  }
  return(list(data=d,Bayes_risk=Bayes_risk))
}

