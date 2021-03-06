% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_predict.R
\name{KMforCSD}
\alias{KMforCSD}
\title{Kernel Machine for Current Status Data - Training}
\usage{
KMforCSD(
  data,
  cost = 1,
  kernel = "rbfdot",
  gamma = 1,
  scale = 1,
  offset = 1,
  degree = 1,
  order = 1,
  alpha_cutoff = 1e-05,
  g_unknown = TRUE,
  misspecification = FALSE
)
}
\arguments{
\item{data}{A data frame consisting of the p-dimensional covariates Z, the current status indicator vector delta, and the censoring times C.}

\item{cost}{A numeric non-negative hyper-parameter. The cost parameter defines a trade-off between model fit and  regularization, and should be fine-tuned for best results. Default value is cost=1. Note that cost=1/(n*lambda), where n is the sample size and lambda is the KM regularization parameter.}

\item{kernel}{A string indicating what type of kernel function should be used. Possible values are: "rbfdot", "polydot", "tanhdot", "vanilladot", "laplacedot", "besseldot", "anovadot", "splinedot". These correspond to \code{\link[kernlab]{dots}}. Default value is kernel="rbfdot" which corresponds to a Gaussian RBF kernel.}

\item{gamma}{A numeric non-negative kernel hyper-parameter; see \code{\link[kernlab]{dots}}. Note that in \code{\link[kernlab]{dots}} gamma is called sigma. We prefer gamma since it is the inverse of twice the RBF kernel width (and sigma is usually reserved for the kernel width).}

\item{scale}{A numeric non-negative kernel hyper-parameter; see \code{\link[kernlab]{dots}}.}

\item{offset}{A numeric kernel hyper-parameter; see \code{\link[kernlab]{dots}}.}

\item{degree}{A positive integer kernel hyper-parameter; see \code{\link[kernlab]{dots}}.}

\item{order}{A numeric kernel hyper-parameter; see \code{\link[kernlab]{dots}}.}

\item{alpha_cutoff}{A small positive numeric value indicating the cutoff for the support vectors. Default value is 0.00001. The support vectors are the covariates that correspond to the non-zero coefficients. Due to numerical errors, we assume that any coefficient with absolute value smaller than the cutoff is essentially zero.}

\item{g_unknown}{Logical - indicating whether the censoring distribution is unknown and needs to be estimated. Default is TRUE and then a kernel density estimate is used to estimate the censoring density g. If FALSE, then the censoring distribution is assumed to be U(0, tau), which only corresponds to the simulated examples in this package.}

\item{misspecification}{Logical - indicating whether the censoring distribution is misspecified. Default is FALSE. If TRUE, then the censoring distribution is estimated via tau\*Beta(0.9,0.9) instead of tau\*Beta(1,1)=U(0, tau), which corresponds to the distribution of the simulated examples in this package.}
}
\value{
The function returns a list with the kernel machine fitted model. The list contains 4 items: (1) the n-dimensional vector of coefficients alpha, (2) the intercept b, (3) the kernel function used during training, and (4) the training covariates.
}
\description{
\code{KMforCSD} returns the KM vector of coefficients, together with the kernel function and the support vectors.
}
\examples{
d <- exp_data(n=100)
data <- d[[1]]
sol <- KMforCSD(data)
}
