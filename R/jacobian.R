#' Check identifiability of models.
#' 
#' Generate random parameter values and test, whether the Jacobian has full rank.
#' 
#' @param type either "2012", "ext", "ext2", or "ext3"
#' @param fixed.theta if TRUE, assumes that mean of theta parameters is assumed to be fixed at zero
#' @param rep number of points in parameter space to check identification
#' @param betas optional matrix with beta parameters
#' @param theta optional matrix with theta parameters (note that the parameters of the first person will be set to to minus the theta-column means if \code{fixed.theta=TRUE})
#' @inheritParams fit_irtree
#' @inheritParams generate_irtree_ext
# @importFrom numDeriv jacobian
#' @examples 
#' # Standard model identified, even without reversed items:
#' jacobian_irtree(N=6, J=8, revItem=rep(1,8), type=2012, rep=5)
#' 
#' @export
jacobian_irtree <- function(N, J, revItem=rep(1,J), traitItem=rep(1,J), 
                            type, fixed.theta=TRUE, rep=100, betas, theta){
    
    n.trait <- length(unique(traitItem))
    S.type <- ifelse(as.character(type) != "2012", 3, 2)
    cols <- S.type+n.trait 
    
    if(!missing(betas) & ! missing(theta))
        rep <- 1
    
    min.rank <- Inf; max.rank <- -Inf
    for(rr in 1:rep){
        if(missing(betas) | rr>1)
            betas <- matrix(rnorm(J*(S.type+1),0,1), J, S.type+1)
        if(missing(theta) | rr>1){
            cols2 <- cols+ifelse(type=="ext3", 1, 0)
            theta <- matrix(rnorm(N*cols2), N, cols2)
        }
        if(fixed.theta)
            theta <- theta[-1,]
        par <- c(betas, theta)
        
        extreme_ARS<-runif(1)
        # if(type %in% c("ext2", "ext3")){
        if(type %in% c("ext", "ext2", "ext3")){
            par <- c(par, extreme_ARS)
        }
        Jac <- numDeriv::jacobian(catprob_irtree_wrapper, x=par, 
                                  N=N, J=J, revItem=revItem, fixed.theta=fixed.theta,
                                  traitItem=traitItem, type=type)
        rk <- qr(Jac)$rank
        min.rank <- min(min.rank, rk)
        max.rank <- max(max.rank, rk)
    }
#     catprob_irtree (betas, theta, revItem=revItem, traitItem=traitItem, type=type, extreme_ARS)
#     catprob_irtree_wrapper(par, N, J, revItem=revItem, traitItem=traitItem, type=type, extreme_ARS)
    
    res <- paste("Model:", type,"(fixed.theta =", fixed.theta,")",
        "\nMinimum rank of Jacobian = ", min.rank, "(Maximum =", max.rank,")", 
        "\nNumber of parameters =", length(par), "(beta and theta)",
        "\n(tested at",rep,"locations in the parameter space)\n")
    
    cat(res)
}


catprob_irtree <- function(betas, theta, revItem=rep(1,J), traitItem=rep(1,J), 
                          type="ext", extreme_ARS){
    N <- nrow(theta)
    J <- nrow(betas)
    # multiple traits
    n.trait <- length(unique(traitItem))
    if(n.trait != 1 & (min(traitItem)!=1 | max(traitItem) != n.trait))
        warning("Check definition of traitItem!")
    S.type <- ifelse(as.character(type) != "2012", 3, 2)
    S <- S.type + n.trait + ifelse(type=="ext3", 1, 0)
    
    
    #     if(missing(theta_vcov) | is.null(theta_vcov)){
    #         theta_vcov <- diag(S)
    #     }else if (is.vector(theta_vcov)){
    #         theta_vcov <- theta_vcov * diag(S)
    #     }else if(any(dim(theta_vcov) != S)){
    #         warning(paste0("check definition of theta_vcov: wrong dimension (required: ", S,")"))
    #     }
    
    # decompose to person ability and item difficulty
    p <- array(NA, c(N, J, 5))
    m <- y <- e <- a <- matrix(0, N, J)
    
    # get predicted category probability
    for(i in 1:N){
        for(j in 1:J){      
            m[i,j] <- pnorm(theta[i,1] - betas[j,1])
            e[i,j] <- pnorm(theta[i,2] - betas[j,2])
            if(type != "2012"){
                a[i,j] <- pnorm(theta[i,3] - betas[j,3])
            }
            y[i,j] <- pnorm(theta[i,S.type+traitItem[j]] - betas[j,S.type+1])
            
            # reversed items: only relevant for trait-dimension/response process
            p.trait <- ifelse(revItem[j] == 1, 1-y[i,j], y[i,j])
            
            # type of ERS in case of ARS:
            extr_ars <- switch(as.character(type),
                               # "ext" =e[i,j], 
                               "ext" =pnorm(theta[i, 2] - qnorm(extreme_ARS)),
                               "ext2"=extreme_ARS, 
                               "ext3"=pnorm(theta[i, S] - qnorm(extreme_ARS)),
                               "2012"=.5)
            
            # response probabilities: MPT model from -2, -1, 0, 1, 2 // 0,1,2,3,4
            p[i,j,1] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*e[i,j]
            p[i,j,2] <- (1-a[i,j])*(1-m[i,j])*(1-p.trait)*(1-e[i,j])
            p[i,j,3] <- (1-a[i,j])*m[i,j]
            p[i,j,4] <- (1-a[i,j])*(1-m[i,j])*p.trait*(1-e[i,j]) +a[i,j]*(1-extr_ars)
            p[i,j,5] <- (1-a[i,j])*(1-m[i,j])*p.trait*e[i,j]     +a[i,j]*extr_ars
            
        }
    }
    
    return(p)
}

catprob_irtree_wrapper <- function(par, N, J, revItem=rep(1,J), traitItem=rep(1,J), 
                                   type="ext", fixed.theta=TRUE){
    n.trait <- length(unique(traitItem))
    S.type <- ifelse(as.character(type) != "2012", 3, 2)
    cols <- S.type+n.trait 
    idx.max <- J*(S.type+1)
    betas <- matrix(par[1:idx.max], J, S.type+1)
    # ll <- length(par) - ifelse(type %in% c("ext","2012"), 0, 1)
    ll <- length(par) - ifelse(type %in% c("2012"), 0, 1)
    if(!fixed.theta){
        theta <- matrix(par[(idx.max+1):ll], N, cols + ifelse(type=="ext3", 1, 0))
    }else{
        theta.tmp <- matrix(par[(idx.max+1):ll], N-1, cols + ifelse(type=="ext3", 1, 0))
        theta <- rbind(-colMeans(theta.tmp), theta.tmp)
    }
    
    extreme_ARS=par[length(par)] # will only be used for ext2 and ext3
    p <- catprob_irtree(betas, theta, revItem=revItem, traitItem=traitItem, 
                       type=type, extreme_ARS=extreme_ARS)
    return(c(p))
}