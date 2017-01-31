#' Generate data for Acquiescence Model.
#'
#' Function generates categorical data 1...5 for \code{N} persons and \code{J} items given the item parameters \code{betas}.
#' 
#' @param N number of persons
#' @param J number of items
#' @param betas Jx4 matrix with item parameters on four response dimensions (middle, extreme, acquiescence, relevant trait defined by \code{traitItem}).
#' @param theta.vcov 4x4 covariance matrix for middle, extremity, acquiescence, trait(s) (can be a vector of length 4 with variances for uncorrelated processes).
#' @param prop.rev proportion of reversed items (rounded to next integer). can be a vector if multiple traits are specified by \code{traitItem}.
#' @param genModel Character. Either \code{"2012"} (Boeckenholt Model without
#'   acquiescence) or \code{"ext"} (Acquiescence Model)
#' @param beta_ARS_extreme only for \code{genModel="ext"}: probability (on probit scale) of choosing category 5 (vs.4) in case of ARS
#' @param cat whether to return categorical data (response categories 1...5) or multinomial data (frequencies of 0 and 1)
#' @inheritParams fit_irtree
#' @return a list containing the generated matrix of responses X, a vector revItem indicating reversed items and true, latent values of the parameters
#' @examples
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' @export
# @import MASS
generate_irtree_ext <- function(N, J, betas, traitItem = rep(1,J), theta.vcov = NULL, 
                          prop.rev = .5, genModel = "ext", beta_ARS_extreme = NULL,
                          cat = TRUE, theta = NULL) {
    
    checkmate::assert_number(beta_ARS_extreme, finite = TRUE)
    
    # multiple traits
    n.trait <- length(unique(traitItem))
    if(n.trait != 1 & (min(traitItem)!=1 | max(traitItem) != n.trait))
        warning("Check definition of traitItem!")
    S_style <- ifelse(as.character(genModel) != "2012", 3, 2)
    S <- S_style + n.trait + ifelse(genModel=="ext3", 1, 0)
    if(length(prop.rev) == 1) {
        prop.rev <- rep(prop.rev, n.trait)
    }
    if (is.null(theta)) {
        if(missing(theta.vcov) | is.null(theta.vcov)) {
            theta.vcov <- diag(S)
        } else if (is.vector(theta.vcov)) {
            theta.vcov <- theta.vcov * diag(S)
        } else if (any(dim(theta.vcov) != S)) {
            warning(paste0("check definition of theta.vcov: wrong dimension (required: ", S,")"))
        }
        
        theta.mu <- rep(0, S)
        theta <- MASS::mvrnorm(N, theta.mu, theta.vcov)
    }
    
    
    # decompose to person ability and item difficulty
    p <- X <- array(NA, c(N, J, 5))
    m <- y <- e <- a <- matrix(NA, N, J)
    # reversed items
    revItem <- rep(0, J)
    for(tt in 1:n.trait){
        Jtmp <-  ceiling(prop.rev[tt] * sum(traitItem == tt))
        tmp <- rep(0:1, c(sum(traitItem == tt)-Jtmp,Jtmp))
        revItem[traitItem == tt] <- sample(  tmp )
    }
    #   if(prop.rev != 0)
    #     revItem[1:ceiling(J*prop.rev)] <- 1
    #   revItem <- sample(revItem, J)
    # generate data
    for(i in 1:N){
        for(j in 1:J){      
            m[i, j] <- pnorm(theta[i,1] - betas[j,1])
            e[i, j] <- pnorm(theta[i,2] - betas[j,2])
            if(as.character(genModel) != "2012"){
                a[i, j] <- pnorm(theta[i,3] - betas[j,3])
            }
            y[i, j] <- pnorm(theta[i, S_style + traitItem[j]] - betas[j, S_style + 1])
            
            # reversed items: only relevant for trait-dimension/response process
            p.trait <- ifelse(revItem[j] == 1, 1-y[i, j], y[i, j])
            
            # type of ERS in case of ARS:
            extr_ars <- switch(as.character(genModel),
                               "ext" = pnorm(theta[i, 2] - beta_ARS_extreme),
                               # "ext" = e[i,j], 
                               # "ext2" = beta_ARS_extreme, 
                               "ext3" = pnorm(theta[i, S] - beta_ARS_extreme),
                               "2012" = .5)
            
            # response probabilities: MPT model from -2, -1, 0, 1, 2 // 0,1,2,3,4
            p[i, j, 1] <- (1-a[i, j])*(1-m[i, j])*(1-p.trait)*e[i, j]
            p[i, j, 2] <- (1-a[i, j])*(1-m[i, j])*(1-p.trait)*(1-e[i, j])
            p[i, j, 3] <- (1-a[i, j])*m[i, j]
            p[i, j, 4] <- (1-a[i, j])*(1-m[i, j])*p.trait*(1-e[i, j])+a[i, j]*(1-extr_ars)
            p[i, j, 5] <- (1-a[i, j])*(1-m[i, j])*p.trait*e[i, j]    +a[i, j]*extr_ars
            
            
            X[i, j, ] <- rmultinom(1, 1, p[i,j,1:5])
        }
    }
    if(cat){
        X <- mult_to_cat(X)
    }

    res <- list(X = X, revItem = revItem, traitItem = traitItem, theta = theta,
                betas = betas, theta.vcov = theta.vcov,
                p = p, middle = m, trait = y, extreme = e, genModel = genModel)
    if(genModel != "2012") res$beta_ARS_extreme <- beta_ARS_extreme
    return(res)
}


#' Generate data for Boeckenholt Model.
#'
#' Function generates categorical data 1...5 for \code{N} persons and \code{J} items given the item parameters \code{betas}.
#' 
#' @param betas Jx3 matrix with item parameters on three response dimensions (middle, extreme, target trait defined by \code{traitItem}).
#' @param theta.vcov 3x3 covariance matrix for middle, extremity, trait(s) (can be a vector of length 3 with variances for uncorrelated processes).
#' @return a list containing the generated matrix of responses X, a vector revItem indicating reversed items and true, latent values of the parameters.
#' @inheritParams fit_irtree
#' @inheritParams generate_irtree_ext
#' @examples
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 0))
#' dat <- generate_irtree_2012(N = N, J = J, betas = betas)
#' @export
generate_irtree_2012 <- function(N, J, betas, traitItem = rep(1, J),
                                 theta.vcov = NULL, prop.rev = .5, cat = TRUE){
    
    
    n.trait <- length(unique(traitItem))
    if(n.trait != 1 & (min(traitItem)!=1 | max(traitItem) != n.trait))
        warning("Check definition of traitItem!")
    S <- 2 + n.trait
    
    if(missing(theta.vcov) | is.null(theta.vcov)){
        theta.vcov <- diag(S)
    }else if (is.vector(theta.vcov)){
        theta.vcov <- theta.vcov * diag(S)
    }
    if(length(prop.rev) == 1)
        prop.rev <- rep(prop.rev, n.trait)
    
    theta.mu <- rep(0, S)
    theta <- MASS::mvrnorm(N, theta.mu, theta.vcov)
    
    # reversed items
    revItem <- rep(0,J)
    for(tt in 1:n.trait){
        Jtmp <-  ceiling(prop.rev[tt] * sum(traitItem == tt))
        tmp <- rep(0:1, c(sum(traitItem == tt)-Jtmp,Jtmp))
        revItem[traitItem == tt] <- sample(  tmp )
    }
    #   if(prop.rev != 0)
    #     revItem[1:ceiling(J*prop.rev)] <- 1
    #   revItem <- sample(revItem, J)
    
    # decompose to person ability and item difficulty
    p <- X <- array(NA, c(N, J, 5))
    m <- y <- e <- matrix(NA, N, J)
    for(i in 1:N){
        for(j in 1:J){      
            m[i,j] <- pnorm(theta[i,1] - betas[j,1])
            e[i,j] <- pnorm(theta[i,2] - betas[j,2])
            y[i,j] <- pnorm(theta[i,2+traitItem[j]] - betas[j,3])
            
            # response probabilities: MPT model from -2, -1, 0, 1, 2
            p[i,j,1] <- (1-m[i,j])*(1-y[i,j])*e[i,j]
            p[i,j,2] <- (1-m[i,j])*(1-y[i,j])*(1-e[i,j])
            p[i,j,3] <- m[i,j]
            p[i,j,4] <- (1-m[i,j])*y[i,j]*(1-e[i,j])
            p[i,j,5] <- (1-m[i,j])*y[i,j]*e[i,j]
            
            if(revItem[j] == 1) 
                p[i,j,1:5] <- rev(p[i,j,1:5])
            X[i,j,] <- rmultinom(1, 1, p[i,j,1:5])
        }
    }
    if(cat){
        X <- mult_to_cat(X)
    }
    
    return(list(X=X, revItem=revItem, traitItem=traitItem, theta=theta, betas=betas, theta.vcov=theta.vcov,
                p=p, middle=m, trait=y, extreme=e))
}

#' Generate item parameters.
#' 
#' Function used internally in \code{\link{recovery_irtree}} to generate item parameters.
#' 
#' @inheritParams recovery_irtree
#' @return A matrix with \code{J} rows containing the item parameters for MRS, ERS,---if \code{genModel != "2012"}---ARS, and trait.
#' @seealso \code{\link[truncnorm]{rtruncnorm}}
# @export
gen_betas <- function(genModel = genModel, J = J, betas = betas) {
    J <- sum(J)
    if (any(!names(betas) %in% c("beta.mrs", "beta.ers", "beta.trait", "beta.ars"))) {
        stop("Argument 'betas' is an optional list with possible entries 'beta.mrs', 'beta.ers', 'beta.trait', and 'beta.ars'. Check definition and spelling of 'betas'.")
    }
    for (iii in seq_len(length(betas))) {
        if (any(!names(betas[[iii]]) %in% c("mean", "sd", "a", "b"))) {
            stop("Possible entries for lists in argument 'betas' are 'mean', 'sd', 'a', and 'b'. Check definition and spelling of 'betas'.")
        }
    }
    res <- matrix(NA, J, ifelse(genModel == "2012", 3, 4))
    # beta.mrs
    if (is.null(betas$beta.mrs$mean)) {
        betas$beta.mrs$mean <- qnorm(.7)
    }
    if (is.null(betas$beta.mrs$sd)) {
        betas$beta.mrs$sd <- sqrt(.1)
    }
    if (is.null(betas$beta.mrs$a)) {
        betas$beta.mrs$a <- qnorm(.5)
    }
    if (is.null(betas$beta.mrs$b)) {
        betas$beta.mrs$b <- qnorm(.9)
    }
    res[, 1] <- do.call(truncnorm::rtruncnorm, c(n = J, betas$beta.mrs))
    # beta.ers
    if (is.null(betas$beta.ers$mean)) {
        betas$beta.ers$mean <- qnorm(.7)
    }
    if (is.null(betas$beta.ers$sd)) {
        betas$beta.ers$sd <- sqrt(.1)
    }
    if (is.null(betas$beta.ers$a)) {
        betas$beta.ers$a <- qnorm(.5)
    }
    if (is.null(betas$beta.ers$b)) {
        betas$beta.ers$b <- qnorm(.9)
    }
    res[, 2] <- do.call(truncnorm::rtruncnorm, c(n = J, betas$beta.ers))
    # beta.trait
    if (is.null(betas$beta.trait$mean)) {
        betas$beta.trait$mean <- 0
    }
    if (is.null(betas$beta.trait$sd)) {
        betas$beta.trait$sd <- sqrt(.5)
    }
    if (is.null(betas$beta.trait$a)) {
        betas$beta.trait$a <- qnorm(.3)
    }
    if (is.null(betas$beta.trait$b)) {
        betas$beta.trait$b <- qnorm(.7)
    }
    res[, ncol(res)] <- do.call(truncnorm::rtruncnorm, c(n = J, betas$beta.trait))
    # beta.ars
    if (genModel != "2012") {
        if (is.null(betas$beta.ars$mean)) {
            betas$beta.ars$mean <- qnorm(.95)
        }
        if (is.null(betas$beta.ars$sd)) {
            betas$beta.ars$sd <- sqrt(.1)
        }
        if (is.null(betas$beta.ars$a)) {
            betas$beta.ars$a <- qnorm(.8)
        }
        if (is.null(betas$beta.ars$b)) {
            betas$beta.ars$b <- qnorm(.999)
        }
        res[, 3] <- do.call(truncnorm::rtruncnorm, c(n = J, betas$beta.ars))
    }
    return(res)
}

#' Convert to Categorical Data
#' 
#' multinomial (1/0 frequencies) to categorical data (responses 1...5)
#' 
#' @param X a N x J x C matrix with response frequencies of 0 and 1
#' @return N x J   matrix with responses 1...5
#' @export
mult_to_cat <- function(X){
    stopifnot(length(dim(X))==3, min(X)==0, max(X)==1, length(table(X))==2)
    dims <- dim(X)
    N <- dims[1]
    J <- dims[2]
    X.cat <- matrix(NA, N, J)
    for(i in 1:N){
        for(j in 1:J){
            X.cat[i,j] <- (1:5)[X[i,j,] ==1]
        }
    }
    return(X.cat)
}

#' Convert to Multinomial Data
#' 
#' categorical  (responses 1...5) to multinomial (1/0 frequencies) data
#' 
#' @param X a N x J   matrix with responses 1...5
#' @return  a N x J x C  matrix with response frequencies of 0 and 1
#' @export
cat_to_mult <- function(X, C=5){
    N <- nrow(X)
    J <- ncol(X)
    X.mult <- array(0, dim=c(N,J,C))
    for(j in 1:J){
        for(i in 1:N){
            X.mult[i,j,X[i,j]] <- 1
        }
    }
    return(X.mult)
}

