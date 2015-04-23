#' Using MPT models to Get IRT Estimates
#'
#' Combination of Multinomial Processing Tree (MPT) Models with Item Response Theory (IRT) to obtain indivdual estimates for different processes. At the moment, standard models for response styles (choosing middle category, extreme responding, acquiesence; see Boeckenholt, 2012) are implemented in JAGS and Stan. 
#'
#' Requires JAGS (\url{http://mcmc-jags.sourceforge.net}) and/or Stan (\url{http://mc-stan.org}) to sample from the posterior using Markov Chain Monte Carlo (MCMC).
#'
#' @docType package
#' @name mpt2irt
#' @author Daniel W. Heck & Hansjoerg Plieninger
#' @references Boeckenholt (2012). Modeling Multiple Response Processes in Judgment and Choice. Psychological Methods, 17, 665-678.
NULL