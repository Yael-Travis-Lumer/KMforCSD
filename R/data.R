#' Online News Popularity Data Set
#'
#' Online news popularity data set after dimensionality reduction and standardization of the covariates.
#'
#' @format A data frame with 39644 rows and 7 variables:
#' \describe{
#'   \item{Z1}{first covariate}
#'   \item{Z2}{second covariate}
#'   \item{Z3}{third covariate}
#'   \item{Z4}{fourth covariate}
#'   \item{Z5}{fifth covariate}
#'   \item{Z6}{sixth covariate}
#'   \item{T}{number of shares}
#' }
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Online+News+Popularity}
"news_pop_data_lasso"


#' Online News Popularity Data Set
#'
#' Censored news popularity data set after dimensionality reduction and standardization of the covariates.
#'
#' @format A data frame with 39644 rows and 8 variables:
#' \describe{
#'   \item{Z1}{first covariate}
#'   \item{Z2}{second covariate}
#'   \item{Z3}{third covariate}
#'   \item{Z4}{fourth covariate}
#'   \item{Z5}{fifth covariate}
#'   \item{Z6}{sixth covariate}
#'   \item{delta}{current status indicator}
#'   \item{C}{censoring time}
#' }
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Online+News+Popularity}
"news_pop_censored"
