#' Plot Observed Frequencies
#' 
#' Plot histograms of observed response frequencies separately for each item. This is especially helpful to judge whether simulated data are reasonable (i.e., unimodal etc.).
#' @param X matrix with observed responses for all participnats (by row) and items (by column)
#' @param revItem use different colors, if specified (as a vector with 1=reversed, 0=not reversed)
#' @param traitItem use as labels for separate traits
#' @param points how many resposne categories in Likert scale
#' @export
plot_responses <- function(X, revItem, traitItem, points=5){
    J <- ncol(X)
    
    if(missing(revItem))
        revItem <- rep(1, J)
    if(missing(traitItem))
        traitItem <- rep(1,J)
    n1 <- floor(sqrt(J) )
    n2 <- ceiling(J/n1)
    
    mfrow <- par()$mfrow
    mar <- par()$mar
    par(mfrow=c(n1,n2), mar=c(4, 3, 3, 1))
    for(j in 1:J){
        barplot(table(factor(X[,j], levels=1:points)), main=paste("Trait", traitItem[j], "(rev:",revItem[j],")"), col=revItem[j]+1)
    }
    
    par(mfrow=mfrow, mar=mar)
    
    return()
}

#' Plot Observed Frequencies and Model Predictions
#' 
#' plot expected distribution of frequencies
#' @param fitModel either "ext" or "2012"
#' @export
boeck_predict <- function(fit_sum = NULL,
                          theta = matrix(0, N, S2),
                          betas = NULL, beta_ARS_extreme = NULL,
                          S = NULL, N = 1, J = NULL, revItem = NULL, traitItem = NULL,
                          fitModel = NULL, item_order = 1:J){
    if (is.null(betas)) betas <- fit_sum$beta
    if (is.null(beta_ARS_extreme)) beta_ARS_extreme <- fit_sum$beta_ARS_extreme
    if (is.null(J)) J <- nrow(fit_sum$beta)
    if (is.null(revItem)) revItem <- fit_sum$revItem
    if (is.null(traitItem)) traitItem <- fit_sum$traitItem
    if (is.null(fitModel)) fitModel <- fit_sum$fitModel
    
    if (is.null(S)) {
        if (fitModel == "ext") {
            S <- 4
        } else if (fitModel == "2012") {
            S <- 3
        }
    }
    S2 <- S - 1 + length(unique(traitItem))
    
    p <- array(NA, c(N, J, 5))
    m <- y <- e <- a <- matrix(NA, N, J)
    if (fitModel == "2012") {
        a <- matrix(1, N, J)
    }

    for(i in 1:N){
        for(j in 1:J){      
            m[i, j] <- pnorm(theta[i, 1] - betas[j, 1])
            e[i, j] <- pnorm(theta[i, 2] - betas[j, 2])
            if (fitModel != "2012") {
                a[i, j] <- pnorm(theta[i, 3] - betas[j, 3])
            }
            y[i, j] <- pnorm(theta[i, (S - 1) + traitItem[j]] - betas[j, S])
            
            # reversed items: only relevant for trait-dimension/response process
            p.trait <- ifelse(revItem[j] == 1, 1 - y[i, j], y[i, j])
            
            # type of ERS in case of ARS:
            # extr_ars <- switch(as.character(fitModel),
            #                    "ext" = e[i,j], 
            #                    "ext2"=extreme_ARS, 
            #                    "ext3"=pnorm(theta[i,S]+qnorm(extreme_ARS)),
            #                    "2012"=.5)
            extr_ars <- switch(as.character(fitModel),
                               "ext" = pnorm(theta[i, 2] - beta_ARS_extreme),
                               # "ext" = e[i,j], 
                               # "ext2" = beta_ARS_extreme, 
                               "ext3" = pnorm(theta[i, S] - beta_ARS_extreme),
                               "2012" = .5)
            
            # response probabilities: MPT model from -2, -1, 0, 1, 2 // 0,1,2,3,4
            p[i,j,1] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*e[i,j]
            p[i,j,2] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*(1-e[i,j])
            p[i,j,3] <- (1-a[i,j])*m[i,j]
            p[i,j,4] <- (1-a[i,j])*(1-m[i,j])*p.trait*(1-e[i,j]) +a[i,j]*(1-extr_ars)
            p[i,j,5] <- (1-a[i,j])*(1-m[i,j])*p.trait*e[i,j]     +a[i,j]*extr_ars
        }
    }
    colnames(p) <- item_order
    px <- reshape2::melt(p, varnames = c("Person", "Item", "Categ"), value.name = "Pred")
    px <- px[order(px$Item, px$Categ), ]
    px$Categ <- factor(px$Categ)
    rownames(px) <- NULL
    if (N == 1) {
        px <- subset(px, select = -Person)
    }
    return(px)
}


#' Plot Expected Frequencies
#' 
#' plot expected distribution of frequencies
#' @param type either "ext" or "2012"
#' @export
plot_expected <- function(model, type, X, revItem, traitItem=rep(1, ncol(X))){
    browser()
  J <- ncol(X)
  S <- ifelse(as.character(type) == "ext", 4, 3)
  
  est <- plot_irtree(fit = model, S = S, N = N, J = J, revItem = revItem,
                 traitItem = traitItem, trait= NULL, return_data = T)
  means <- matrix(est$Mean,J, S)
  
  predict <- function(par, type){
    p <- 
    if(as.character(type) == "ext"){
      
      
      p[i,j,1] <- (1-acq[i,j])*(1-mid[i,j])*(1-trait[i,j])*extr[i,j]
      p[i,j,2] <- (1-acq[i,j])*(1-mid[i,j])*(1-trait[i,j])*(1-extr[i,j])
      p[i,j,3] <- (1-acq[i,j])*mid[i,j]
      p[i,j,4] <- (1-acq[i,j])*(1-mid[i,j])*trait[i,j]*(1-extr[i,j]) +acq[i,j]*(1-extr[i,j])
      p[i,j,5] <- (1-acq[i,j])*(1-mid[i,j])*trait[i,j]*extr[i,j]     +acq[i,j]*extr[i,j]
    }else{
      
    }
  }
  
  n1 <- floor(sqrt(J) )
  n2 <- ceiling(J/n1)
  
  mfrow <- par()$mfrow
  mar <- par()$mar
  par(mfrow=c(n1,n2), mar=c(4, 3, 3, 1))
  for(j in 1:J){
    barplot(table(factor(X[,j], levels=1:5)), main=paste("Trait", traitItem[j], "(rev:",revItem[j],")"), col=revItem[j]+1)
  }
  
  par(mfrow=mfrow, mar=mar)
}