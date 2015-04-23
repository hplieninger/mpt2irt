#' Plot True and Estimated Parameters
#'
#' Plot estimates from parameter recovery simulation. Only two object jags.samp/stan.samp/true.par are allowed
#' 
#' @param N number of participants
#' @param J number of items
#' @param number of theta parameters to plot (either 3 or 4)
#' @param jags.samp a fitted mcmc.list from JAGS
#' @param stan.samp a fitted Stan array
#' @param true.par a list with true values for beta and theta
#' @export
plot_recovery <- function(N, J, S, jags.samp, stan.samp, true.par){
  
  if(S == 4)
    parnames <- c("middle","agree","extremity", "acquies")
  else
    parnames <- c("middle","agree","extremity")
  
  ###################### correlation to expected values
  theta.jags <- theta.stan <- matrix(NA, N, S); 
  beta.jags <- beta.stan <- betaSE.jags <- betaSE.stan <- matrix(NA, J, S)
  plotsel <- rep(F, 3)
  
  if(!missing(true.par)){
    plotsel[1] <- T
    theta.true <- true.par$theta.true
    beta.true <- true.par$beta.true  
  }
  
  if(!missing(jags.samp)){
    plotsel[2] <- T
    for(i in 1:S){
      theta.jags[,i] <- summary(boeck.samp[,paste0("theta[",1:N,",",i,"]")])$stat[,"Mean"]
      beta.jags[,i] <- summary(boeck.samp[,paste0("beta[",1:J,",",i,"]")])$stat[,"Mean"]
    }
  }
  
  if(!missing(stan.samp)){
    plotsel[3] <- T
    for(i in 1:S){
      theta.stan[,i] <- summary(boeck.stansamp)$summ[paste0("theta[",1:N,",",i,"]"),1]
      beta.stan[,i] <- summary(boeck.stansamp)$summ[paste0("beta[",1:J,",",i,"]"),1]
    }
  }
  if(sum(plotsel) != 2)
    warning("please provide only two objects for plotting")
  
  if(plotsel[1]){
    x1 <- theta.true; x2 <- beta.true
    xl <- "true"
    if(plotsel[2]){
      y1 <- theta.jags; y2 <- beta.jags;
      yl <- "JAGS"
    }else{
      y1 <- theta.stan; y2 <- beta.stan;
      yl <- "Stan"
    }}else{
      x1 <- theta.jags; x2 <- beta.jags;
      y1 <- theta.stan; y2 <- beta.stan;
      xl <- "JAGS"; yl <- "Stan"
    }
  
  for(i in 1:S){
    plot(x1[,i],y1[1:N,i], main=paste("Theta",parnames[i]), xlab=xl)
  }
  for(s in 1:S){
    plot(x2[,i],y2[1:J,i], col=revItem+10, main=paste("Beta",parnames[i]), ylab=yl)
  }
  
  par(mfrow=c(1,1))
  
}
