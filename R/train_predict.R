#' Kernel Machine for Current Status Data - Training
#'
#' \code{KMforCSD} returns the KM vector of coefficients, together with the kernel function and the support vectors.
#'
#' @param data A data frame consisting of the p-dimensional covariates Z, the current status indicator vector delta, and the censoring times C.
#' @param cost A numeric non-negative hyper-parameter. The cost parameter defines a trade-off between model fit and  regularization, and should be fine-tuned for best results. Default value is cost=1. Note that cost=1/(n*lambda), where n is the sample size and lambda is the KM regularization parameter.
#' @param kernel A string indicating what type of kernel function should be used. Possible values are: "rbfdot", "polydot", "tanhdot", "vanilladot", "laplacedot", "besseldot", "anovadot", "splinedot". These correspond to \code{\link[kernlab]{dots}}. Default value is kernel="rbfdot" which corresponds to a Gaussian RBF kernel.
#' @param gamma A numeric non-negative kernel hyper-parameter; see \code{\link[kernlab]{dots}}. Note that in \code{\link[kernlab]{dots}} gamma is called sigma. We prefer gamma since it is the inverse of twice the RBF kernel width (and sigma is usually reserved for the kernel width).
#' @param scale A numeric non-negative kernel hyper-parameter; see \code{\link[kernlab]{dots}}.
#' @param offset A numeric kernel hyper-parameter; see \code{\link[kernlab]{dots}}.
#' @param degree A positive integer kernel hyper-parameter; see \code{\link[kernlab]{dots}}.
#' @param order A numeric kernel hyper-parameter; see \code{\link[kernlab]{dots}}.
#' @param alpha_cutoff A small positive numeric value indicating the cutoff for the support vectors. Default value is 0.00001. The support vectors are the covariates that correspond to the non-zero coefficients. Due to numerical errors, we assume that any coefficient with absolute value smaller than the cutoff is essentially zero.
#' @param g_unknown Logical - indicating whether the censoring distribution is unknown and needs to be estimated. Default is TRUE and then a kernel density estimate is used to estimate the censoring density g. If FALSE, then the censoring distribution is assumed to be U(0, tau), which only corresponds to the simulated examples in this package.
#' @param misspecification Logical - indicating whether the censoring distribution is misspecified. Default is FALSE. If TRUE, then the censoring distribution is estimated via tau\*Beta(0.9,0.9) instead of tau\*Beta(1,1)=U(0, tau), which corresponds to the distribution of the simulated examples in this package.
#' @return The function returns a list with the kernel machine fitted model. The list contains 4 items: (1) the n-dimensional vector of coefficients alpha, (2) the intercept b, (3) the kernel function used during training, and (4) the training covariates.
#' @examples
#' d <- exp_data(n=100)
#' data <- d[[1]]
#' sol <- KMforCSD(data)
#' @export
KMforCSD <- function(data,cost=1,kernel="rbfdot",gamma=1,scale=1,offset=1,degree=1,order=1, alpha_cutoff=0.00001,g_unknown=TRUE,misspecification=FALSE){
  n <- nrow(data)
  p <- ncol(data)-2
  lambda <- 1/(n*cost)
  delta <- data$delta
  C <- data$C
  cols_remove <- c("C","delta")
  Z <- data[, !(colnames(data) %in% cols_remove)]
  Z <- as.matrix(Z)
  tau <- max(C)
  ## define dot product function
  dotproduct <- switch(kernel,
                       "rbfdot" = kernlab::rbfdot(sigma=gamma),
                       "polydot" = kernlab::polydot(degree,scale,offset),
                       "tanhdot" = kernlab::tanhdot(scale,offset),
                       "vanilladot" = kernlab::vanilladot(),
                       "laplacedot" = kernlab::laplacedot(sigma = gamma),
                       "besseldot" = kernlab::besseldot(sigma = gamma, order, degree),
                       "anovadot" = kernlab::anovadot(sigma = gamma, degree),
                       "splinedot" = kernlab::splinedot()
  )
  ## calculate kernel
  K <- kernlab::kernelMatrix(dotproduct, Z)@.Data
  Q <- K+n*lambda*diag(n)
  vec1 <- matrix(1,n,1)
  vec2 <- rbind(vec1,0)
  A <- cbind(Q,vec1)
  A <- rbind(A,t(vec2))
  if (!isTRUE(misspecification)){
    if (isTRUE(g_unknown)){
      g_hat <- ks::kde(C,eval.points = C)$estimate
    } else {
      g_hat <- 1/tau
    }
  } else {
    g_hat <- (1/tau)*stats::dbeta(C/tau,0.9,0.9)
  }
  V <- as.matrix((1-delta)/g_hat,n,1)
  V <- rbind(V,0)
  sol <- solve(A,V)
  alpha <- sol[1:n]
  b <- sol[n+1]
  SVs_ind=which(abs(alpha)>alpha_cutoff,arr.ind = TRUE)
  alpha=alpha[SVs_ind]
  SVs=Z[SVs_ind,]
  model = structure(list(alpha=alpha,b=b,used_kernel=dotproduct,support_vectors=SVs), class = "KMforCSD")
  return(model)
}


#' Kernel Machine for Current Status Data - Prediction
#'
#' \code{predict.KMforCSD} returns the KM decision function, together with the test data (which can be equivalent to the training data).
#'
#' @param object A KMforCSD fitted model. This is the output of \code{\link{KMforCSD}}.
#' @param newdata A data frame consisting of the p-dimensional covariates Z, the current status indicator vector delta, and the censoring times C.
#' @return The function returns a data frame consisting of newdata, and the kernel machine predictions.
#' @examples
#' d <- exp_data(n=100) #training set
#' sol <- KMforCSD(data=d[[1]]) #training
#' d_test <- exp_data(n=50) #test set
#' new_data <- predict(sol,d_test[[1]]) #prediction
#' decision_function <- new_data$prediction
#' @export
predict.KMforCSD <- function(object,newdata,...){
  n_new <- nrow(newdata)
  p <- ncol(newdata)-2
  cols_remove <- c("C","delta","T")
  Z_new <- newdata[, !(colnames(newdata) %in% cols_remove)]
  Z_new <- as.matrix(Z_new)
  alpha <- object$alpha
  b <- object$b
  dotproduct <- object$used_kernel
  SVs <- object$support_vectors
  SVs <- as.matrix(SVs,length(SVs),p)
  K_new <- kernlab::kernelMatrix(dotproduct,SVs ,Z_new)@.Data
  prediction <- t(K_new)%*%alpha+b
  return_data <- newdata
  return_data$response <- prediction
  return(return_data)
}
