#' Generate Data for Extended Boeckenholt (2012)
#'
#' Generate data
#' 
#' @param N number of participants
#' @param J number of items
#' @param beta Jx4 matrix with item properties on four response dimensions (middle, extreme, acquiescence, relevant trait defined by traitItem)
#' @param traitItem vector of length J specifying latent traits, e.g., for Big5 from 1...5 . standard: unidimensional test, all items measure the same trait.
#' @param theta.vcov 4x4 covariance matrix for middle, extremity, acquiesence, trait(s) (can be a vector of length 4 with variances for uncorrelated processes)
#' @param prop.rev proportion of reversed items (rounded to next integer)
#' @param cat whether to return categorical data (response categories 1...5) or multinomial data (frequencies of 0 and 1)
#' @return a list containing the generated matrix of responses X, a vector revItem indicating reversed items and true, latent values of the parameters
#' @examples
#' # generate data for unidimensional scale
#' N <- 20; J <- 10
#' beta <- cbind(runif(J), seq(-2,2, length=J), runif(J), runif(J))
#' gendat <- boeck.gen.ext(N=N, J=J, beta=beta)
#' @export
#' @import MASS
boeck.gen.ext <- function(N, J, beta, traitItem=rep(1,J), theta.vcov=NULL, prop.rev=.5, cat=T){
  
  n.trait <- length(unique(traitItem))
  if(n.trait != 1 & (min(traitItem)!=1 | max(traitItem) != n.trait))
    warning("Check definition of traitItem!")
  S <- 3 + n.trait
  
  if(missing(theta.vcov) | is.null(theta.vcov)){
    theta.vcov <- diag(S)
  }else if (is.vector(theta.vcov)){
    theta.vcov <- theta.vcov * diag(S)
  }
  library(MASS)
  theta.mu <- rep(0,S)
  theta <- mvrnorm(N, theta.mu, theta.vcov)
  
  # decompose to person ability and item difficulty
  p <- X <- array(NA, c(N, J, 5))
  m <- y <- e <- a <- matrix(NA, N, J)
  # reversed items
  revItem <- rep(0,J)
  if(prop.rev != 0)
    revItem[1:ceiling(J*prop.rev)] <- 1
  revItem <- sample(revItem, J)
  # generate data
  for(i in 1:N){
    for(j in 1:J){      
      m[i,j] <- pnorm(theta[i,1] - beta[j,1])
      e[i,j] <- pnorm(theta[i,2] - beta[j,2])
      a[i,j] <- pnorm(theta[i,3] - beta[j,3])
      y[i,j] <- pnorm(theta[i,3+traitItem[j]] - beta[j,4])
      
      # reversed items: only relevant for trait-dimension/response process
      p.trait <- ifelse(revItem[j] == 1, 1-y[i,j], y[i,j])
      
      # response probabilities: MPT model from -2, -1, 0, 1, 2 // 0,1,2,3,4
      p[i,j,1] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*e[i,j]
      p[i,j,2] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*(1-e[i,j])
      p[i,j,3] <- (1-a[i,j])*m[i,j]
      p[i,j,4] <- (1-a[i,j])*(1-m[i,j])*p.trait*(1-e[i,j]) +a[i,j]*(1-e[i,j])
      p[i,j,5] <- (1-a[i,j])*(1-m[i,j])*p.trait*e[i,j]     +a[i,j]*e[i,j]
      
      
      X[i,j,] <- rmultinom(1, 1, p[i,j,1:5])
    }
  }
  if(cat){
    X <- mult.to.cat(X)
  }
  
  return(list(X=X, revItem=revItem, traitItem=traitItem, theta=theta, beta=beta, theta.vcov=theta.vcov,
              p=p, middle=m, trait=y, extreme=e))
}


#' Generate Data for Original Boeckenholt (2012)
#'
#' More general function: \link{boeck.gen.ext}
#' 
#' @param N number of participants
#' @param J number of items
#' @param theta.vcov covariance matrix for middle, extremity, trait(s) (can be a vector of length 3 for uncorrelated processes)
#' @param beta Jx3 matrix with item properties on three response dimensions (middle, extreme, relevant trait defined by traitItem)
#' @param traitItem vector of length J specifying latent traits, e.g., for Big5 from 1...5. standard: unidimensional test, all items measure the same trait.
#' @param prop.rev proportion of reversed items (rounded to next integer)
#' @param cat whether to return categorical data (response categories 1...5) or multinomial data (frequencies of 0 and 1)
#' @return a list containing the generated matrix of responses X, a vector revItem indicating reversed items and true, latent values of the parameters
#' @examples
#' # generate data
#' N <- 20; J <- 10
#' beta <- cbind(runif(J,0,1), seq(-2,2, length=J), runif(J,0,1))
#' gen <- boeck.gen.2012(N=N, J=J, beta=beta)
#' @export
boeck.gen.2012 <- function(N, J, beta, traitItem=rep(1,J), theta.vcov=NULL, prop.rev=.5, cat=T){
  
  
  n.trait <- length(unique(traitItem))
  if(n.trait != 1 & (min(traitItem)!=1 | max(traitItem) != n.trait))
    warning("Check definition of traitItem!")
  S <- 2 + n.trait
  
  if(missing(theta.vcov) | is.null(theta.vcov)){
    theta.vcov <- diag(S)
  }else if (is.vector(theta.vcov)){
    theta.vcov <- theta.vcov * diag(S)
  }
  
  library(MASS)
  theta.mu <- rep(0, S)
  theta <- mvrnorm(N, theta.mu, theta.vcov)
  
  # reversed items
  revItem <- rep(0,J)
  if(prop.rev != 0)
    revItem[1:ceiling(J*prop.rev)] <- 1
  revItem <- sample(revItem, J)
  
  # decompose to person ability and item difficulty
  p <- X <- array(NA, c(N, J, 5))
  m <- y <- e <- matrix(NA, N, J)
  for(i in 1:N){
    for(j in 1:J){      
      m[i,j] <- pnorm(theta[i,1] - beta[j,1])
      e[i,j] <- pnorm(theta[i,2] - beta[j,2])
      y[i,j] <- pnorm(theta[i,2+traitItem[j]] - beta[j,3])
      
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
    X <- mult.to.cat(X)
  }
  
  return(list(X=X, revItem=revItem, traitItem=traitItem, theta=theta, beta=beta, theta.vcov=theta.vcov,
              p=p, middle=m, trait=y, extreme=e))
}


#' Convert to Categorical Data
#' 
#' multinomial (1/0 frequencies) to categorical data (responses 1...5)
#' 
#' @param X a N x J x C matrix with response frequencies of 0 and 1
#' @return N x J   matrix with responses 1...5
#' @export
mult.to.cat <- function(X){
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
cat.to.mult <- function(X, C=5){
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

