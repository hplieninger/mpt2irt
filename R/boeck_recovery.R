################################################
#### Recovery Simulation
################################################

#' Recovery Simulation for Response Styles
#' 
#' Repeatedly generates and fits data and returns a list of the true and fitted estimates. At the moment, only unidimensional traits are supported (i.e., no response scales like the Big 5)
#' 
#' @param N sample size
#' @param J number of items. Can be a vector for multiple traits (e.g., J=c(10,10,10))
#' @param prop.rev number of reversed items. Can be a vector for multiple traits(e.g., prop.rev=c(.5,.3,.5))
#' @param R number of replications
#' @param genModel data generating model (either "2012" or "ext")
#' @param fitModel model for data analysis ("2012", "ext", or both as vector c("2012", "ext"))
#' @param fitMethod whether to use "stan" or "jags"
#' @param trait.range range of normally distrroundeibuted item properties for trait (positive!)
#' @param RSmin minimum of range for uniformly distributed response styles (can be a vector with separate bounds for middle/extreme/(acquiesence))
#' @param RSmax maximum of range for response styles (can also be a vector)
#' @param theta.vcov true covariance matrix of response processes (order: middle, extreme, (acquiescence), trait). standard is diag(3) / diag(4)
#' @param df degrees of freedom for wishart hyperprior on trait covariance matrix. standard ist the number of trait parameters plus one. should be a vector of length two if two models are fitted.
#' @param V prior matrix for theta covariance. standard is diag(S). if two models are fitted, should be a list with two matrices of the correct size
#' @param M number of MCMC samples (adaptation + burnin = M/2)
#' @param n.chains number of chains
#' @param thin thinning of MCMC samples
#' @param nCPU number of CPUs for parallel run of MCMC sampler
# @param printProgress whether to print the progress in a text file (either at the location defined by \code{path} or at the working directory)
#' @param path where to save progress report (e.g., path="C:/"). disabled by using path=""
#' @param saveTemp save temporary results to same location as progress report
# @param rng_seed random number seed
# @param ... further arguments passed to stan() or coda.samples() for JAGS
#' @return a list of class "boeck.recov" containing true and estimated parameters
#' @import rjags
#' @import coda
#' @import rstan
#' @import doParallel
#' @import runjags
#' @import parallel
#' @import foreach
#' @examples
#' \dontrun{
#'  sim <- boeck.recovery(N=30, J=6, R=2, genModel="2012", fitModel=c("2012","ext"), 
#'                    fitMethod="jags", nCPU=4, RSmin=0, RSmax=1, prop.rev=.5)
#'  summary(sim, whichFit=1, plot=T)
#'  sim$dic
#'}
#' @export
boeck.recovery <- function(N, J, R, genModel, fitModel, fitMethod, 
                           prop.rev=.5, trait.range=1.5, RSmin=0, RSmax=1, 
                           theta.vcov, df, V, M=1000, n.chains=4, thin=1, 
                           nCPU=4, path=getwd(), saveTemp=F){
  
  genModel <- as.character(genModel)
  fitModel <- as.character(fitModel)
  S.gen <- ifelse(genModel == "ext", 4, 3)
  S.fit <- sapply(fitModel, function(x) ifelse(x == "ext", 4, 3))
  
  if(length(RSmin) == 1)
    RSmin <- rep(RSmin, S.gen-1)
  if(length(RSmax) == 1)
    RSmax <- rep(RSmax, S.gen-1)
  if(length(RSmin)!=(S.gen-1) | length(RSmax)!=(S.gen-1))
    warning("Check definition of RSmin and RSmax")
  
  genNames <- c(paste0("theta[",apply(expand.grid(1:N, 1:S.gen),1, 
                                      paste0, collapse=","),"]"),
                paste0("beta[",apply(expand.grid(1:J, 1:S.gen),1, 
                                     paste0, collapse=","),"]"), 
                paste0("Sigma[",apply(expand.grid(1:S.gen, 1:S.gen),1, 
                                      paste0, collapse=","),"]"))
  ttt <- c("T_obs","T_pred","post_p")
  
  # separate lists for fitting fitModel
  fitNames <- theta_mu <- V.ls <- list()
  for(sss in 1:length(S.fit)){
    fitNames[[sss]] <- c(paste0("theta[",apply(expand.grid(1:N, 1:S.fit[sss]),1, 
                                               paste0, collapse=","),"]"),
                         paste0("beta[",apply(expand.grid(1:J, 1:S.fit[sss]),1, 
                                              paste0, collapse=","),"]"), 
                         paste0("Sigma[",apply(expand.grid(1:S.fit[sss], 1:S.fit[sss]),1, 
                                               paste0, collapse=","),"]"),
                         ttt)
    theta_mu[[sss]] <- rep(0,S.fit[sss])
    
    ### set up parameters
    if(missing(V)){
      V.ls[[sss]] <- diag(S.fit[sss])
    }else if(length(S.fit) == 1){
      V.ls[[sss]] <- unlist(V)
    }else{
      V.ls[[sss]] <- V[[sss]]
    }
  }
  if(missing(df)){
    df <- S.fit+1
  }
  
  if(missing(theta.vcov)){
    theta.vcov <- diag(S.gen)  # middle  / extremity / acq / trait
  }else if(is.vector(theta.vcov)){
    theta.vcov <- diag(S.gen) * theta.vcov   # middle /  extremity / trait
  }
  if(any(dim(theta.vcov)!= S.gen)){
    warning("Check definition of theta.vcov!")
  }
  
  
  #### function for a single CPU ####################
  
  boeckrep <- function(rr){
    
    #### generate data
    if(genModel == "2012"){
      beta <- cbind(runif(J,RSmin[1],RSmax[1]),
                    runif(J,RSmin[2],RSmax[2]),
                    qnorm(seq(.05,.95, length=J),0,trait.range ))
      gen <- boeck.gen.2012(N=N, J=J, beta=beta, theta.vcov=theta.vcov, prop.rev=prop.rev, cat=T)
    }else if (genModel == "ext"){
      beta <- cbind(runif(J,RSmin[1],RSmax[1]),
                    runif(J,RSmin[2],RSmax[2]), 
                    runif(J,RSmin[3],RSmax[3]),
                    qnorm(seq(.05,.95, length=J),0,trait.range ))
      gen <- boeck.gen.ext(N=N, J=J, beta=beta, theta.vcov=theta.vcov, prop.rev=prop.rev, cat=T)
    }
    
    
    
    ##################### fit JAGS / Stan
    fitpar <- list() ; dic <- list()
    for(sss in 1:length(S.fit)){
      datalist <- list(S=S.fit[sss], df=df[sss], V=V.ls[[sss]], N=N, J=J, revItem=gen$revItem,
                       X=gen$X, theta_mu=theta_mu[[sss]])
      if(fitMethod == "jags"){
        #         boeck.jags <- run.jags(model=paste0( modelPath,"/mpt2irt/models/jags_boeck_",fitModel[sss],"_1d.txt"), 
        #                               monitor = c("theta", "beta", "Sigma", "T_obs","T_pred", "post_p",
        #                                           'deviance', 'pd', 'pd.i', 'popt', 'dic'), 
        #                               data=datalist, n.chains=n.chains,  
        #                               burnin = M/4, sample = M, adapt=M/4, 
        #                               summarise = T, thin = thin, method="simple")
        boeck.jags <- jags.model(file=paste0( modelPath,"/mpt2irt/models/jags_boeck_",
                                              fitModel[sss],"_1d.txt"), 
                                 data=datalist,  n.chains=n.chains, n.adapt=M/4)
        adapt <- F
        while(!adapt){
          adapt <- adapt(boeck.jags,M/4,end.adaptation = FALSE)
        }
        update(boeck.jags, M/4)
        boeck.samp <- coda.samples.dic(boeck.jags, variable.names = c("theta", "beta", "Sigma",
                                                                      "T_obs","T_pred", "post_p"),
                                       n.iter=M*thin, thin = thin)
        dic[[sss]] <- boeck.samp$dic
        boeck.stan <- mcmc.list2stan(boeck.samp$samples)
        rm(boeck.samp)
      }else{
        ### Stan
        boeck.stan <- stan(file=paste0( modelPath,"/mpt2irt/models/stan_boeck_",
                                        fitModel[sss],"_1d.stan"), 
                           model_name=paste0("Boeckenholt_", fitModel[sss]),
                           pars=c("theta","beta", "Sigma","T_obs","T_pred","post_p"),
                           data=datalist, chains = n.chains, thin=thin,
                           iter=M*thin)
        dic[[sss]] <- NULL
      }
      
      fitpar[[sss]] <- monitor(boeck.stan, print=FALSE)[fitNames[[sss]],] 
    }
    genpar <- c(gen$theta, gen$beta, gen$theta.vcov)
    names(genpar) <- genNames
    
    if(path != ""){
      try(write(paste0("\n|", paste0(rep("#",floor(rr/R*30)), collapse=""),
                       paste0(rep(" ",30-floor(rr/R*30)), collapse=""),
                       "| rr=",rr,"/",R," ; ",Sys.time()), append=T,
                file=paste0(path,"/boeck_recovery_progress.txt")))
      if(saveTemp)
        try(save(genpar, fitpar, dic, file=paste0(path,"/boeck_recov_",fitMethod,"_gen-",genModel,
                     "_fit-",paste0(fitModel, collapse="-") ,"_N=",N,"_J=",J,"_rep",rr,".RData")))
    }
    
    # clean up to get RAM
    gc(T, verbose=F)
    res <- list(gen=genpar, fit=fitpar, dic=dic)
    return(res)
  }
  
  # paths to save temporary results
  modelPath <- .libPaths()[1]
  if(path != "")
    try(write(paste0("Boeckenholt Recovery Simulation: \ngenModel = ", genModel, "; fitModel = ", 
                     paste0(fitModel, collapse="-"), " ; Start: ",Sys.time()),
              file=paste0(path,"/boeck_recovery_progress.txt")))
  
  ######## parallel computing
  
  cl <- makeCluster(nCPU)
  registerDoParallel(cl)
  
  runtime <- system.time({
    res <- foreach(rr=1:R, .packages = c("mpt2irt","rjags","rstan","coda"), .combine="c") %dopar% {boeckrep(rr)}
  })
  cat("\nTotal Runtime: ",round(runtime[3]/3600, 2), "Hours\n")
  
  stopCluster(cl)
  
  ###### collect results
  
  estimate1 <-  array(NA, c(R, N*S.fit[1] + J*S.fit[1] + S.fit[1]*S.fit[1] + 3 ,10), 
                      dimnames= list(NULL, fitNames[[1]], colnames(res[[2]][[1]])))
  if(length(S.fit) == 2){
    estimate2 <-  array(NA, c(R, N*S.fit[2] + J*S.fit[2] + S.fit[2]*S.fit[2] + 3 ,10), 
                        dimnames= list(NULL, fitNames[[2]], colnames(res[[2]][[2]])))
  }else{
    estimate2 <- NULL
  }
  true <- matrix(NA, R, length(res[[1]]),  dimnames=list(NULL, names(res[[1]])))
  dic <- data.frame(dev1=NA, penalty1=NA, dic1=NA, dev2=NA, penalty2=NA, dic2=NA)
  
  for(r in 1:R){
    true[r,] <- res[[(r-1)*3+1]]
    estimate1[r,,] <- res[[(r-1)*3+2]][[1]]
    ddd <- res[[(r-1)*3+3]]
    try(dic[r,1:3] <- c(ddd[[1]]$dev, ddd[[1]]$pen, ddd[[1]]$dev+ ddd[[1]]$pen), silent=T)
    if(length(S.fit) == 2){
      estimate2[r,,] <- res[[(r-1)*3+2]][[2]]
      try(dic[r,4:6] <- c(ddd[[2]]$dev, ddd[[2]]$pen, ddd[[2]]$dev+ ddd[[2]]$pen), silent=T)
    }
  }
  
  if(path != ""){
    try(unlink(paste0(path,"/boeck_recovery_progress.txt")), silent=T)
    for(rr in 1:R){
      try(unlink(paste0(path,"/boeck_recov_",fitMethod,"_gen-",genModel,
        "_fit-",paste0(fitModel, collapse="-") ,"_N=",N,"_J=",J,"_rep",rr,".RData")), silent=T)
    }
  }
  
  res <- list(estimate1=estimate1, estimate2=estimate2, true=true, dic=dic, N=N, J=J, S.gen=S.gen, 
              S.fit=S.fit, R=R, genModel=genModel, fitModel=fitModel, fitMethod=fitMethod, 
              RSmin=RSmin, RSmax=RSmax, theta.vcov=theta.vcov,
              prop.rev=prop.rev, trait.range=trait.range, thin=thin, M=M)
  class(res) <- "boecksim"
  return(res)
}

#' Summarize Results of Recovery Simulation
#' 
#' Compute correlations and mean difference of true and estimated values in recovery simulation.
#' 
#' @param sim simulation generated using \link{boeck.recovery}
#' @param whichFit whether to summarize for first or second fitted model
#' @param plot whether to plot the estimated and true values
#' @param rep which replication to plot
#' @param whether to plot estimated and true parameters
#' @param save if specified, plot will be saved in working directoy
#' @return a list containing matrices with correlations and mean differences for theta and beta of each replication
#' @export
summary.boecksim <- function(sim, whichFit=1, plot=T, rep=1, save=F){
  J <- sim$J
  N <- sim$N
  R <- sim$R
  S.fit <- sim$S.fit[whichFit]
  S.gen <- sim$S.gen
  S <- min(S.fit, S.gen)
  
  if(S == 3)
    parnames <- c("middle", "extreme", "trait")
  else if (S == 4)
    parnames <- c("middle", "extreme", "acquiesence","trait")
  
  cor.beta <- cor.theta <- meandiff.beta <- meandiff.theta <- matrix(NA, R, S, dimnames=list(NULL, parnames))
  mfrow <- par()$mfrow
  
  if(whichFit == 1){
    estimate <- sim$estimate1
  }else{
    estimate <- sim$estimate2
  }
  
  if(save)
    png(paste0("gen-",sim$genModel,"_fit-",sim$fitModel[whichFit],"_N=",N,"_J=",J,"_rep=",rep,".png"), 800,800)
  par(mfrow=c(S,4))
  for(rr in 1:R){
    for(s in 1:S){
      # correct indices for cases : generate "2012" / fit "ext" ....
      selb.est <- paste0("beta[",1:J, ",",s+ifelse(S==3 & S.fit==4 & s==3,1,0),"]")
      selt.est <- paste0("theta[",1:N, ",",s+ifelse(S==3 & S.fit==4 & s==3,1,0),"]")
      selb.true <- paste0("beta[",1:J, ",",s+ifelse(S==3 & S.gen==4 & s==3,1,0),"]")
      selt.true <- paste0("theta[",1:N, ",",s+ifelse(S==3 & S.gen==4 & s==3,1,0),"]")
      
      
      meandiff.theta[,s] <- rowMeans(estimate[,selt.est,1] - sim$true[,selt.true])
      meandiff.beta[,s] <- rowMeans(estimate[,selb.est,1] - sim$true[,selb.true])
      
      cor.theta[rr,s] <- cor(estimate[rr,selt.est,1], sim$true[rr,selt.true])
      cor.beta[rr,s] <- cor(estimate[rr,selb.est,1], sim$true[rr,selb.true])
      
      if(rr == rep & plot){
        hist(estimate[rep,selt.est,1] - sim$true[rep,selt.true],40, 
             main=paste("person mean difference=",round(meandiff.theta[rep,s],3)),
             xlab=paste(parnames[s],"(est-true)"))
        plot(sim$true[rep,selt.true],estimate[rep,selt.est,1], 
             main=paste("person cor =",round(cor.theta[rep,s],3)),
             xlab="true",ylab="estimated")
        hist(estimate[rep,selb.est,1] - sim$true[rep,selb.true],40,
             main=paste("item mean difference =",round(meandiff.beta[rep,s],3)),
             xlab=paste(parnames[s],"(est-true)"))
        plot(sim$true[rep,selb.true], estimate[rep,selb.est,1], 
             main=paste("item cor =",round(cor.beta[rep,s],3)),
             xlab="true",ylab="estimated")
      }
    }
  }
  if(save)
    dev.off()
  par(mfrow=mfrow)
  list(cor.theta=cor.theta, cor.beta=cor.beta,
       meandiff.theta=meandiff.theta, meandiff.beta=meandiff.beta)
}
