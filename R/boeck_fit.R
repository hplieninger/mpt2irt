#' Analyze with Boeckenholt Model for Response Styles
#' 
#' Using either the standard 2012 model (only response styles towards middle + extreme responding) or the extended version (additionally, aquiescence).
#' 
#' The estimated parameters are ordered as follows:
#' \itemize{
#' \item Model "2012": S=2+number of traits
#' \itemize{
#'  \item theta[i,1:S] = c(middle, extreme, trait(s))
#'  \item beta[j,1:3] = c(middle, extreme, trait) (which trait depends on traitItem!)
#'  }
#' \item Model "ext": S=3+number of traits
#' \itemize{
#'  \item theta[i,1:S] = c(middle, extreme, acquiesence, trait(s))
#'  \item beta[j,1:4] = c(middle, extreme, acquiesence, trait) (which trait depends on traitItem!)
#'  }
#' }
#' If more than a single trait is measured, beta and theta have more columns accordingly (e.g., theta[i,1:5]=c(mid, extr, acq, trait1,..., trait5))
#' 
#' @param X an N x J matrix of observed responses for categories 1...5 (use \link{mult2cat} to transform a multinomial frequency matrix with 1s/0s to responses from 1...5)
#' @param revItem vector of length J specifying reversed items (1=reversed, 0=not reversed)
#' @param traitItem vector of length J specifying the underlying traits (e.g., indexed from 1...5). Standard: only a single trait is measured by all items. If the Big5 are measures, might be something like c(1,1,1,2,2,2,...,5,5,5,5)
#' @param df degrees of freedom for wishart prior on covariance of traits (standard/minimum: number of processes + 1)
#' @param items either "fixed" or "random" (with hierarchical normal-wishart structure estimated from the data)
#' @param V prior for wishart distribution (standard: diagonal matrix)
#' @param fitModel either "2012" (without acquiescence) or "ext"
#' @param fitMethod whether to use JAGS or Stan
#' @param format either "mcmc.list" (can be analyzed with coda package) or "stan" 
#' @param startSmall  only for JAGS: whether to use random starting values for beta sampled from "wide" (startSmall=F, generated within JAGS) or "narrow" priors (startSmall=T; beta and theta closer to 0; might solve problems with slow convergence of some chains for extreme starting values).
#' @param M number of MCMC samples (after warmup)
#' @param warmup number of samples for warmup (in JAGS: 3/4 for adaption, 1/4 for burnin)
#' @param n.chains number of MCMC chains (and number of CPUs used)
#' @param thin thinning of MCMC samples
#' @param ... further arguments passed to stan() or coda.samples() for JAGS
#' @details  Note that the progress of Stan is shown in a text file in the working directory ("_Stanprogress.txt")
#' @examples 
#' \dontrun{
#' # generate data
#' N <- 20; J <- 10
#' beta <- cbind(runif(J,0,1), seq(-2,2, length=J), runif(J,0,1))
#' gen <- boeck.gen.2012(N=N, J=J, beta=beta, cat=T)
#' 
#' # fit model
#' fit <- boeck.fit(gen$X, gen$revItem, M= 200, n.chains=2, 
#'                  fitMethod="jags", format="mcmc.list")
#' summary(fit[,paste0("beta[",1:J,",1]")])$stat
#' plot(fit[,paste0("beta[",1:3,",1]")])
#' }
#' @import runjags
#' @import rstan
#' @import parallel
#' @export
boeck.fit <- function(X, revItem, traitItem=rep(1, ncol(X)), df, V, items="fixed",
                      fitModel="ext", fitMethod="jags", format="stan", startSmall=T,
                      M=1000, warmup=1000, n.chains=2, thin=1, ...){
  fitModel <- as.character(fitModel)
  N <- nrow(X)
  J <- ncol(X)
  if(!is.vector(revItem) | length(revItem) != J)
    warning("Check definition of vector revItem!")
  n.trait <- length(table(traitItem))
  if(min(traitItem) != 1 | max(traitItem) != n.trait){
    warning("Check definition of traitItem")
  }
  # number of response processes/dimensions
  S <- switch(fitModel, "ext" = 3,  "2012"=2)+n.trait
  if(missing(df)){
    df <- S+1
  }
  if(missing(V)){
    V <- diag(S)
  }
  
  # adjust starting values
  inits <- list()
  if(startSmall){
    if(items == "fixed"){
      for(c in 1:n.chains){
        inits[[c]] <- list(beta = matrix(rnorm( J*(S-n.trait+1),0, .1), J) )
      }
    }else if (items == "random"){
      for(c in 1:n.chains){
        inits[[c]] <- list(beta.raw = matrix(rnorm( J*S,0, .1), J), 
                           xi.beta = rnorm(S, 1, .05))
      }
    }
  }
  
  if(items == "fixed"){
    datalist <- list(X=X, S=S, df=df, V=V, J=J, N=N, theta_mu=rep(0,S), 
                     revItem=revItem, traitItem=traitItem)
    varlist <-  c("theta", "beta", "Sigma", "T_obs","T_pred","post_p")
    
  }else if (items == "random"){
    datalist <- list(X=X, S=S, df=df, V=V, J=J, N=N, theta_mu=rep(0,S), V0.beta=diag(S),
                     revItem=revItem, traitItem=traitItem) #, selItem=cbind(matrix(1:3, J, 3, byrow=T), traitItem+3)
    varlist <-  c("theta", "beta", "Sigma.beta", "Sigma.theta", "T_obs","T_pred","post_p")
  }
  
  if(fitMethod == "jags"){
    ##################### fit JAGS ###############################################
    boeck.jags <- run.jags(model=paste0(.libPaths()[1],"/mpt2irt/models/jags_boeck_",
                                        fitModel,"_",ifelse(items=="fixed", "5d", "v2"), ".txt"), 
                           monitor = varlist, 
                           data=datalist, inits=inits,
                           n.chains=n.chains, burnin = ceiling(warmup/4), sample = M, adapt=ceiling(warmup*3/4),
                           thin = thin, method='parallel', summarise=FALSE, modules=c("glm"))
    
    boeck.samp <- as.mcmc.list(boeck.jags)
    if(format == "stan")
      boeck.samp <- mcmc.list2stan(boeck.samp)
    
  }else{  
    ##################### fit Stan ###############################################
    
    boeck.samp <- stan_parallel(file=paste0(.libPaths()[1],"/mpt2irt/models/stan_boeck_",
                                            fitModel,"_", ifelse(items=="fixed", "5d", "v2"), ".stan"), 
                                model_name = paste0("Boeckenholt_",  
                                                    fitModel,"_",n.trait,"d_items=",items),
                                data = datalist, pars = varlist, 
                                chains = n.chains, iter = warmup+M*thin, 
                                warmup=warmup, thin = thin,
                                path=paste0(getwd(),"\\"),
                                init = 'random')
    try( unlink(paste0(getwd(), "\\_StanProgress.txt")), silent=T)
    if(format == "mcmc.list"){
      boeck.samp <- stan2mcmc.list(boeck.samp)
    }
  }
  
  return(boeck.samp)
}