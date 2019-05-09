#' Sample Posterior Predictive Responses
#'
#' This function samples posterior predictive responses (x=1,...,K) given a
#' matrix of probabilities. It is implemented in C++ to increase speed.
#'
#' @param prob Numeric matrix with K columns
#' @name sample_pp
#' @return Numeric vector of length \code{nrow(prob)} containing integerish
#'   responses between 1 and K
#' @examples
#' prob <- matrix(runif(5*20), 20)
#' mpt2irt:::sample_pp(prob)
NULL
