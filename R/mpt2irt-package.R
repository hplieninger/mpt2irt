#' mpt2irt: Bringing Multinomial Processing Tree Models To Item Response Theory To Investigate Response Styles
#'
#' Multinomial processing tree (MPT) models are used to differentiate
#' between different cognitive processes in categorical data. Those can be
#' lifted to hierarchical MPT models to account for person- and/or item effects
#' (Matzke et al., 2015): This is done by reparameterizing each MPT parameter
#' by means of an item response (IRT) model. Boeckenholt (2012) developed a
#' model to disentengle two reponse styles, namely midpoint (MRS) and extreme
#' responding (ERS), and the target trait. This package allows to fit an
#' extension of this Boeckenholt Model to acquiescence (ARS). It makes use of
#' Bayesian hierarchical models requiring JAGS or Stan.
#'
#' The package requires installation of JAGS
#' (\url{http://mcmc-jags.sourceforge.net}) and/or Stan
#' (\url{http://mc-stan.org}) to sample from the posterior.
#'
#' @docType package
#' @name mpt2irt-package
#' @aliases mpt2irt
#' @author Hansjoerg Plieninger & Daniel W. Heck
#' @references Plieninger, H., & Heck, D. W. (2018). A new model for
#'   acquiescence at the interface of psychometrics and cognitive psychology.
#'   Multivariate Behavioral Research. Advance online publication.
#'   doi:10.1080/00273171.2018.1469966
#' @useDynLib mpt2irt, .registration = TRUE 
#' 
#' @import methods
#' @importFrom rstan optimizing sampling vb constrain_pars extract
#'   extract_sparse_parts get_posterior_mean stanc
#' @import Rcpp
#' @import rstantools
#' @importFrom utils menu object.size sessionInfo setTxtProgressBar txtProgressBar tail combn
#' @import stats
#' @importFrom grDevices heat.colors
#' @importFrom graphics abline barplot hist legend lines par plot
#' @importFrom magrittr %>%
#' 
#' @examples 
#' \dontrun{
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' 
#' # fit model
#' res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200)
#' res2 <- summarize_irtree_fit(res1)
#' res3 <- tidyup_irtree_fit(res2, N = N, J = J, revItem = dat$revItem,
#'                           traitItem = dat$traitItem, fitModel = res$fitModel)
#' res3$plot
#' }
#' 
NULL
