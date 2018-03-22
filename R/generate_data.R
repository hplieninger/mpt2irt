#' Generate data for Acquiescence Model.
#'
#' Function generates categorical data 1...5 for \code{N} persons and \code{J}
#' items given the item parameters \code{betas}.
#' 
#' @param N number of persons
#' @param J number of items
#' @param betas Jx4 matrix with item parameters on four response dimensions
#'   (middle, extreme, acquiescence, relevant trait defined by
#'   \code{traitItem}).
#' @param theta_vcov 4x4 covariance matrix for middle, extremity, acquiescence,
#'   trait(s) (can be a vector of length 4 with variances for uncorrelated
#'   processes).
#' @param prop.rev proportion of reversed items (rounded to next integer). Can
#'   be a vector if multiple traits are specified by \code{traitItem}. Only used
#'   if \code{revItem = NULL}.
#' @param revItem vector of length J specifying reversed items (1=reversed,
#'   0=not reversed). Overrides argument \code{prop.rev}.
#' @param genModel Character. Either \code{"2012"} (Boeckenholt Model without
#'   acquiescence) or \code{"ext"} (Acquiescence Model)
#' @param beta_ARS_extreme only for \code{genModel="ext"}: probability (on
#'   probit scale) of choosing category 5 (vs.4) in case of ARS
#' @param cat whether to return categorical data (response categories 1...5) or
#'   multinomial data (frequencies of 0 and 1)
#' @param theta Numeric. Optional matrix with \code{N} rows containing the true
#'   person parameters theta.
#' @inheritParams fit_irtree
#' @return The function returns a list containing the generated matrix of
#'   responses X, a vector revItem indicating reversed items and true, latent
#'   values of the parameters.
#' @examples
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' @export
# @import MASS
generate_irtree_ext <- function(N = NULL,
                                J = NULL,
                                betas = NULL,
                                traitItem = rep(1,J),
                                theta_vcov = NULL, 
                                prop.rev = .5,
                                revItem = NULL,
                                genModel = c("ext", "ext2"),
                                beta_ARS_extreme = NULL,
                                cat = TRUE,
                                theta = NULL) {
    genModel <- match.arg(genModel)
    
    checkmate::qassert(N, "X1[1,]")
    checkmate::qassert(J, "X>0[1,]")
    checkmate::assert_matrix(betas, mode = "double", any.missing = FALSE, 
                             nrows = J, ncols = ifelse(genModel == "ext", 4, 5))
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 len = J)
    # checkmate::qassert(prop.rev, "N1[0,1]")
    # checkmate::assert_numeric(prop.rev, lower = 0, upper = 1, any.missing = FALSE,
    #                           min.len = 1, max.len = length(unique(traitItem)))
    checkmate::assert_number(beta_ARS_extreme, finite = TRUE,
                             null.ok = ifelse(genModel == "ext", F, T))
    
    # multiple traits
    n.trait <- length(unique(traitItem))
    if (n.trait != 1 & (min(traitItem) != 1 | max(traitItem) != n.trait))
        warning("Check definition of traitItem!")
    S_style <- ifelse(as.character(genModel) != "2012", 3, 2)
    S <- S_style + n.trait + ifelse(genModel == "ext3", 1, 0)
    # if (length(prop.rev) == 1) {
    #     prop.rev <- rep(prop.rev, n.trait)
    # }
    if (is.null(theta)) {
        if (missing(theta_vcov) | is.null(theta_vcov)) {
            theta_vcov <- diag(S)
        } else if (is.vector(theta_vcov)) {
            theta_vcov <- theta_vcov * diag(S)
        } else if (any(dim(theta_vcov) != S)) {
            warning(paste0("check definition of theta_vcov: wrong dimension (required: ", S,")"))
        }
        
        theta.mu <- rep(0, S)
        theta <- MASS::mvrnorm(N, theta.mu, theta_vcov)
    }
    
    
    # decompose to person ability and item difficulty
    p <- X <- array(NA, c(N, J, 5))
    m <- y <- e <- a <- matrix(NA, N, J)
    
    # # reversed items
    revItem <- create_revItem(prop.rev = prop.rev, revItem = revItem, 
                              J = J, n.trait = n.trait,
                              traitItem = traitItem)
    # revItem <- rep(0, J)
    # for (tt in 1:n.trait) {
    #     Jtmp <-  ceiling(prop.rev[tt] * sum(traitItem == tt))
    #     tmp <- rep(0:1, c(sum(traitItem == tt) - Jtmp,Jtmp))
    #     revItem[traitItem == tt] <- sample(  tmp )
    # }
    
    # generate data
    for (i in 1:N) {
        for (j in 1:J) {      
            m[i, j] <- pnorm(theta[i,1] - betas[j,1])
            e[i, j] <- pnorm(theta[i,2] - betas[j,2])
            if (as.character(genModel) != "2012") {
                a[i, j] <- pnorm(theta[i,3] - betas[j,3])
            }
            y[i, j] <- pnorm(theta[i, S_style + traitItem[j]] - betas[j, S_style + 1])
            
            # reversed items: only relevant for trait-dimension/response process
            p.trait <- ifelse(revItem[j] == 1, 1 - y[i, j], y[i, j])
            
            # type of ERS in case of ARS:
            extr_ars <- switch(as.character(genModel),
                               "ext" = pnorm(theta[i, 2] - beta_ARS_extreme),
                               # "ext" = e[i,j], 
                               "ext2" = pnorm(theta[i, 2] - betas[j, 5]),
                               # "ext3" = pnorm(theta[i, S] - beta_ARS_extreme),
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
                betas = betas, theta_vcov = theta_vcov,
                p = p, middle = m, trait = y, extreme = e, genModel = genModel)
    if(genModel != "2012") res$beta_ARS_extreme <- beta_ARS_extreme
    return(res)
}


#' Generate data for Boeckenholt Model.
#'
#' Function generates categorical data 1...5 for \code{N} persons and \code{J} items given the item parameters \code{betas}.
#' 
#' @param betas Jx3 matrix with item parameters on three response dimensions (middle, extreme, target trait defined by \code{traitItem}).
#' @param theta_vcov 3x3 covariance matrix for middle, extremity, trait(s) (can be a vector of length 3 with variances for uncorrelated processes).
#' @return The function returns a list containing the generated matrix of responses X, a vector revItem indicating reversed items and true, latent values of the parameters.
#' @inheritParams fit_irtree
#' @inheritParams generate_irtree_ext
#' @examples
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 0))
#' dat <- generate_irtree_2012(N = N, J = J, betas = betas)
#' @export
generate_irtree_2012 <- function(N = NULL,
                                 J = NULL,
                                 betas = NULL,
                                 traitItem = rep(1, J),
                                 theta_vcov = NULL,
                                 prop.rev = .5,
                                 revItem = NULL,
                                 cat = TRUE){
    
    checkmate::qassert(N, "X1[2,]")
    checkmate::qassert(J, "X>0[1,]")
    checkmate::assert_matrix(betas, mode = "double", any.missing = FALSE, 
                             nrows = J, ncols = 3)
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 len = J)
    # checkmate::qassert(prop.rev, "N1[0,1]")
    # checkmate::assert_numeric(prop.rev, lower = 0, upper = 1, any.missing = FALSE,
    #                           min.len = 1, max.len = length(unique(traitItem)))
    
    
    n.trait <- length(unique(traitItem))
    if(n.trait != 1 & (min(traitItem)!=1 | max(traitItem) != n.trait))
        warning("Check definition of traitItem!")
    S <- 2 + n.trait
    
    if(missing(theta_vcov) | is.null(theta_vcov)){
        theta_vcov <- diag(S)
    }else if (is.vector(theta_vcov)){
        theta_vcov <- theta_vcov * diag(S)
    }
    # if(length(prop.rev) == 1)
    #     prop.rev <- rep(prop.rev, n.trait)
    
    theta.mu <- rep(0, S)
    theta <- MASS::mvrnorm(N, theta.mu, theta_vcov)
    
    # reversed items
    revItem <- create_revItem(prop.rev = prop.rev, revItem = revItem, 
                              J = J, n.trait = n.trait,
                              traitItem = traitItem)
    # revItem <- rep(0,J)
    # for(tt in 1:n.trait){
    #     Jtmp <-  ceiling(prop.rev[tt] * sum(traitItem == tt))
    #     tmp <- rep(0:1, c(sum(traitItem == tt)-Jtmp,Jtmp))
    #     revItem[traitItem == tt] <- sample(  tmp )
    # }
    
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
    if (cat) {
        X <- mult_to_cat(X)
    }
    
    return(list(X=X, revItem=revItem, traitItem=traitItem, theta=theta, betas=betas, theta_vcov=theta_vcov,
                p=p, middle=m, trait=y, extreme=e))
}

#' Generate data for Steps Model.
#'
#' Function generates categorical data 1...5 for \code{N} persons and \code{J}
#' items with paramters from \code{rnorm()}.
#' 
#' @section References:
#' 
#' Tutz, G. (1997). Sequential models for ordered responses. In W. J. van der
#' Linden & R. K. Hambleton (Eds.), Handbook of Modern Item Response Theory (pp.
#' 139-152). doi:10.1007/978-1-4757-2691-6_8
#' 
#' Verhelst, N. D., Glas, C. A. W., & de Vries, H. H. (1997). A steps model to 
#' analyze partial credit. In W. J. van der Linden & R. K. Hambleton (Eds.), 
#' Handbook of modern item response theory (pp. 123-138).
#' doi:10.1007/978-1-4757-2691-6_7
#' 
#' @inheritParams generate_irtree_ext
#' @return The function returns a list containing the generated matrix of
#'   responses X, a vector revItem indicating reversed items and true, latent
#'   values of the parameters.
#' @examples
#' N <- 20
#' J <- 10
#' dat <- generate_irtree_steps(N = N, J = J)
#' @export
generate_irtree_steps <- function(N = NULL,
                                  J = NULL,
                                  prop.rev = .5,
                                  revItem = NULL,
                                  traitItem = NULL) {
    
    checkmate::qassert(N, "X1[1,)")
    checkmate::qassert(J, "X1[1,)")
    checkmate::assert_integerish(revItem, lower = 0, upper = 1,
                                 any.missing = FALSE, len = J, null.ok = TRUE)
    checkmate::assert_integerish(traitItem, lower = 1, # upper = 1,
                                 any.missing = FALSE, len = J, null.ok = TRUE)
    if (is.null(traitItem)) traitItem <- rep(1, J)
    # if (is.null(revItem)) revItem <- rbinom(J, 1, .33)
    revItem <- create_revItem(prop.rev = prop.rev, revItem = revItem, 
                              J = J, n.trait = length(unique(traitItem)),
                              traitItem = traitItem)
    
    thres <- matrix(rnorm(J*4), J, 4) %>%
        apply(1, sort) %>%
        t
    theta <- matrix(rnorm(N*max(traitItem)), N, max(traitItem))
    
    p_catx <- array(NA, dim = c(N, J, 5))
    dat <- node1 <- node2 <- node3 <- node4 <- matrix(NA, N, J)
    
    for (i in 1:N) {	
        for (j in 1:J) {
            node1[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 1]);
            node2[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 2]);
            node3[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 3]);
            node4[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 4]);
            
            p_catx[i,j,1] = (1-node1[i, j]);
            p_catx[i,j,2] =    node1[i, j] *(1-node2[i, j]);
            p_catx[i,j,3] =    node1[i, j] *   node2[i, j] *(1-node3[i, j]);
            p_catx[i,j,4] =    node1[i, j] *   node2[i, j] *   node3[i, j] *(1-node4[i, j]);
            p_catx[i,j,5] =    node1[i, j] *   node2[i, j] *   node3[i, j] *   node4[i, j] ;
        }
    }
    
    p_cat <- p_catx
    
    p_cat[ , revItem == 1, 5] <- p_catx[ , revItem == 1, 1]
    p_cat[ , revItem == 1, 4] <- p_catx[ , revItem == 1, 2]
    p_cat[ , revItem == 1, 3] <- p_catx[ , revItem == 1, 3]
    p_cat[ , revItem == 1, 2] <- p_catx[ , revItem == 1, 4]
    p_cat[ , revItem == 1, 1] <- p_catx[ , revItem == 1, 5]
    
    dat[, ] <- p_cat %>%
        apply(1:2, function(x) rmultinom(1, 1, x)) %>%
        magrittr::equals(1) %>% 
        apply(2:3, which)
    
    return(list(X = dat,
                theta = theta,
                thres = thres,
                revItem = revItem,
                traitItem = traitItem))
}

#' @title Generate data for Shift-Tree Model.
#'
#' @description Function generates categorical data 1...5 for \code{N} persons
#'   and \code{J} items given the item parameters \code{betas}.
#'
#' @param betas Jx3 matrix with item parameters on three response dimensions
#'   (middle, extreme, and relevant trait defined by \code{traitItem}).
#' @inheritParams generate_irtree_ext
#' @inheritParams fit_irtree
#' @return The function returns a list containing the generated matrix of
#'   responses X, the true parameters, as well as further specifications such as
#'   \code{traitItem} or \code{revItem}.
#' @examples
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 0))
#' dat <- generate_irtree_shift(N = N, J = J, betas = betas)
#' @export
generate_irtree_shift <- function(N = NULL,
                                  J = NULL,
                                  betas = NULL,
                                  traitItem = rep(1, J),
                                  theta_vcov = NULL, 
                                  prop.rev = .5,
                                  revItem = NULL,
                                  genModel = "shift",
                                  # beta_ARS_extreme = NULL,
                                  cat = TRUE,
                                  theta = NULL) {
    
    checkmate::qassert(N, "X1[2,]")
    checkmate::qassert(J, "X>0[1,]")
    checkmate::assert_matrix(betas, mode = "double", any.missing = FALSE, 
                             nrows = J, ncols = 3)
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 len = J)
    # checkmate::assert_numeric(prop.rev, lower = 0, upper = 1, any.missing = FALSE,
    #                           len = length(unique(traitItem)))
    
    # multiple traits
    n.trait <- length(unique(traitItem))
    if (n.trait != 1 & (min(traitItem) != 1 | max(traitItem) != n.trait))
        warning("Check definition of traitItem!")
    S_style <- ifelse(as.character(genModel) != "2012", 3, 2)
    S <- S_style + n.trait
    # if (length(prop.rev) == 1) {
    #     prop.rev <- rep(prop.rev, n.trait)
    # }
    if (is.null(theta)) {
        if (missing(theta_vcov) | is.null(theta_vcov)) {
            theta_vcov <- diag(S)
        } else if (is.vector(theta_vcov)) {
            theta_vcov <- theta_vcov * diag(S)
        } else if (any(dim(theta_vcov) != S)) {
            warning(paste0("check definition of theta_vcov: wrong dimension (required: ", S,")"))
        }
        
        theta.mu <- rep(0, S)
        theta <- MASS::mvrnorm(N, theta.mu, theta_vcov)
    }
    
    
    # decompose to person ability and item difficulty
    p <- X <- array(NA, c(N, J, 5))
    m <- y <- e <- a <- matrix(NA, N, J)
    
    # reversed items
    revItem <- create_revItem(prop.rev = prop.rev, revItem = revItem, 
                              J = J, n.trait = n.trait,
                              traitItem = traitItem)
    # revItem <- rep(0, J)
    # for (tt in 1:n.trait) {
    #     Jtmp <-  ceiling(prop.rev[tt] * sum(traitItem == tt))
    #     tmp <- rep(0:1, c(sum(traitItem == tt) - Jtmp, Jtmp))
    #     revItem[traitItem == tt] <- sample(  tmp )
    # }
    
    # generate data
    for (i in 1:N) {
        for (j in 1:J) {      
            m[i, j] <- pnorm(theta[i, 1] - betas[j, 1])
            e[i, j] <- pnorm(theta[i, 2] - betas[j, 2])
            # if(as.character(genModel) != "2012"){
            #     a[i, j] <- pnorm(theta[i, 3] - betas[j, 3])
            # }
            y[i, j] <- pnorm(theta[i, 3] +
                                 (-1)^revItem[j] * (theta[i, S_style+traitItem[j]] - betas[j, 3]))
            
            # reversed items: only relevant for trait-dimension/response process
            # p.trait <- ifelse(revItem[j] == 1, 1-y[i, j], y[i, j])
            
            # response probabilities: MPT model from -2, -1, 0, 1, 2 // 0,1,2,3,4
            p[i, j, 1] <- (1-m[i, j])*(1-y[i, j])*   e[i, j]
            p[i, j, 2] <- (1-m[i, j])*(1-y[i, j])*(1-e[i, j])
            p[i, j, 3] <-    m[i, j]
            p[i, j, 4] <- (1-m[i, j])*   y[i, j] *(1-e[i, j])
            p[i, j, 5] <- (1-m[i, j])*   y[i, j] *   e[i, j]
            
            
            X[i, j, ] <- rmultinom(1, 1, p[i, j, 1:5])
        }
    }
    if (cat) {
        X <- mult_to_cat(X)
    }
    
    res <- list(X = X, revItem = revItem, traitItem = traitItem, theta = theta,
                betas = betas, theta_vcov = theta_vcov,
                p = p, middle = m, trait = y, extreme = e, genModel = genModel)
    # if(genModel != "2012") res$beta_ARS_extreme <- beta_ARS_extreme
    return(res)
}

#' Generate item parameters.
#' 
#' Function used internally in \code{\link{recovery_irtree}} to generate item parameters.
#' 
#' @inheritParams recovery_irtree
#' @return A matrix with \code{J} rows containing the item parameters for MRS,
#'   ERS, ARS (if \code{genModel != "2012"}), and the target trait.
#' @seealso \code{\link[truncnorm]{rtruncnorm}}
#' @examples
#' J <- 10
#' # use defaults
#' betapar <- mpt2irt:::gen_betas("ext", J = J, betas = NULL)
#' dat <- generate_irtree_ext(N = 20, J = J, betas = betapar, beta_ARS_extreme = 1)
#' 
#' # modify distribution (truncated normal) from which to draw betas, here for MRS
#' tmp1 <- list("beta.mrs" = list("mean" =  0,
#'                                "sd"   =  0.3,
#'                                "a"    = -2,
#'                                "b"    =  Inf))
#' betapar <- mpt2irt:::gen_betas("ext", J = J, betas = tmp1)
# @export
gen_betas <- function(genModel = c("ext", "2012", "pcm", "steps", "shift", "ext2"),
                      J = NULL,
                      betas = NULL) {
    
    genModel <- match.arg(genModel)
    checkmate::qassert(J, "X>0[1,]")
    checkmate::assert_list(betas, null.ok = TRUE)
    
    J <- sum(J)
    if (any(!names(betas) %in% c("beta.mrs", "beta.ers", "beta.trait", "beta.ars", "beta.aers"))) {
        stop("Argument 'betas' is an optional list with possible entries 'beta.mrs', ",
             "'beta.ers', 'beta.trait', 'beta.ars', and 'beta.aers'. Check ",
             "definition and spelling of 'betas'.")
    }
    for (iii in seq_len(length(betas))) {
        if (any(!names(betas[[iii]]) %in% c("mean", "sd", "a", "b"))) {
            stop("Possible entries for lists in argument 'betas' are 'mean', 'sd', 'a', and 'b'. Check definition and spelling of 'betas'.")
        }
    }
    S_b1 <- switch(genModel,
                   "2012"  = 3,
                   "ext"   = 4,
                   "pcm"   = 4,
                   "steps" = 4,
                   "shift" = 3,
                   "ext2"  = 5)
    
    res <- matrix(NA, J, S_b1)
    
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
    tmp1 <- ifelse(genModel %in% c("ext2"), 4, ncol(res))
    res[, tmp1] <- do.call(truncnorm::rtruncnorm, c(n = J, betas$beta.trait))
    
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
    
    # beta.aers
    if (genModel %in% c("ext2")) {
        if (is.null(betas$beta.aers$mean)) {
            betas$beta.aers$mean <- qnorm(.7)
        }
        if (is.null(betas$beta.aers$sd)) {
            betas$beta.aers$sd <- sqrt(.1)
        }
        if (is.null(betas$beta.aers$a)) {
            betas$beta.aers$a <- qnorm(.5)
        }
        if (is.null(betas$beta.aers$b)) {
            betas$beta.aers$b <- qnorm(.9)
        }
        res[, 5] <- do.call(truncnorm::rtruncnorm, c(n = J, betas$beta.aers))
    }
    
    return(res)
}

#' Convert to categorical data.
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

#' Convert to multinomial data.
#' 
#' categorical  (responses 1...5) to multinomial (1/0 frequencies) data
#' 
#' @param X a N x J   matrix with responses 1...5
#' @param C Number of categories.
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

#' @title Internal function to create revItem from prop.rev
#' @description Creates \code{revItem}  from \code{prop.rev} or simply returns \code{revItem} if supplied.
#' @param n.trait Integer, number of content-related traits.
#' @inheritParams generate_irtree_ext
#' @return vector \code{revItem} comprised of 0s and 1s.
create_revItem <- function(prop.rev = NULL, revItem = NULL, J = NULL, n.trait = NULL,
                           traitItem = NULL) {
    if (!is.null(revItem)) {
        checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                     len = J)
        # return(revItem)
    } else {
        checkmate::assert_numeric(prop.rev, lower = 0, upper = 1, any.missing = FALSE,
                                  min.len = 1, max.len = n.trait)
        if (length(prop.rev) == 1) {
            prop.rev <- rep(prop.rev, n.trait)
        }
        revItem <- rep(0, J)
        for (tt in 1:n.trait) {
            Jtmp <-  ceiling(prop.rev[tt] * sum(traitItem == tt))
            tmp <- rep(0:1, c(sum(traitItem == tt) - Jtmp, Jtmp))
            revItem[traitItem == tt] <- sample(tmp)
        }
    }
    return(revItem)
}
