#' Fit an mpt2irt model.
#' 
#' This function fits an mpt2irt model. Either the so-called Boeckenholt Model
#' can be fit (\code{fitModel = "2012"}) that assumes the three processes MRS,
#' ERS, and target trait; or the so-called Acquiescence Model can be fit that
#' additionally takes ARS into account (\code{fitModel = "ext"}).
#' 
#' The estimated parameters are arranged as follows:
#' \itemize{
#' \item Model "2012" (S = 2 + number of traits):
#' \itemize{
#'  \item theta[i, 1:S] = c(middle, extreme, trait(s))
#'  \item beta[j, 1:3] = c(middle, extreme, trait) (which trait depends on \code{traitItem}!)
#'  }
#' \item Model "ext" (S = 3 + number of traits):
#' \itemize{
#'  \item theta[i, 1:S] = c(middle, extreme, acquiesence, trait(s))
#'  \item beta[j, 1:4] = c(middle, extreme, acquiesence, trait) (which trait depends on \code{traitItem}!)
#'  }
#' }
#' If more than a single trait is measured, theta has more columns accordingly
#' (e.g., theta[i,1:6]=c(mid, extr, acq, trait1,..., trait3))
#' 
#' @param X an N x J matrix of observed responses for categories 1...5 (use
#'   \code{\link{mult_to_cat}} to transform a multinomial frequency matrix with 1s/0s to
#'   responses from 1...5)
#' @param revItem vector of length J specifying reversed items (1=reversed,
#'   0=not reversed)
#' @param traitItem vector of length J specifying the underlying traits (e.g.,
#'   indexed from 1...5). Standard: only a single trait is measured by all
#'   items. If the Big5 are measures, might be something like
#'   c(1,1,1,2,2,2,...,5,5,5,5)
#' @param df degrees of freedom for wishart prior on covariance of traits
#'   (standard/minimum: number of processes + 1)
# param items either "fixed" or "random" (with hierarchical normal-wishart structure estimated from the data)
#' @param V prior for wishart distribution (standard: diagonal matrix)
#' @param fitModel Character. Either \code{"2012"} (Boeckenholt Model without
#'   acquiescence) or \code{"ext"} (Acquiescence Model) or \code{"pcm"} (partial
#'   credit model).
#  or \code{"ext2"} (separate probability of choosing category 5 in case of ARS; requires at least two trait scales) or \code{"ext3"} (separate person-estimates theta[i,S] for latent tendency to choose cat.4/5 in case of ARS).
#' @param fitMethod whether to use JAGS or Stan
#' @param outFormat either "mcmc.list" (can be analyzed with coda package) or
#'   "stan" or "runjags"
#' @param startSmall Whether to use random starting values for
#'   beta sampled from "wide" (startSmall=F, generated within JAGS) or "narrow"
#'   priors (startSmall=T; beta and theta closer to 0; might solve problems with
#'   slow convergence of some chains for extreme starting values).
#' @param M number of MCMC samples (after warmup)
#' @param warmup number of samples for warmup (in JAGS: 3/4 for adaption, 1/4
#'   for burnin)
#' @param n.chains number of MCMC chains (and number of CPUs used)
#' @param thin thinning of MCMC samples
# @param mail email address to which a notification is sent when the simulation is finished (no dash "-" allowed!)
#' @param return_defaults Logical. Whether to return input specifications or not.
#' @param add2varlist Additional variables to monitor (e.g., \code{c("deviance",
#'   "pd", "popt", "dic")} for JAGS)
#' @param N2 Numeric. Number of persons for whom to draw posterior predictives.
#'   Specify equal to \code{nrow(X)} in order to draw values for all persons.
#'   This is mainly implemented for efficiency reasons in order to avoid
#'   massivly drawing samples which the user is not interested in.
#' @param ... further arguments passed to \code{\link[rstan]{sampling}} (for Stan) or \code{\link[runjags]{run.jags}} (for JAGS)
# @details  Note that the progress of Stan is shown in a text file in the
#'   working directory ("_Stanprogress.txt")
#' @inheritParams runjags::run.jags
# @import rstan
#' @return Returns a list where the output form either JAGS or Stan is stored in the entry \code{samples}.
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
#'                           traitItem = dat$traitItem, fitModel = res1$fitModel)
#' str(res3)
#' res3$plot
#' }
# import runjags
# @importFrom rstan stan sampling As.mcmc.list
# @import parallel
# @importFrom mail sendmail
#' @export
fit_irtree <- function(X,
                       revItem = NULL,
                       traitItem = rep(1, ncol(X)),
                       df = NULL,
                       V = NULL, 
                       fitModel = c("ext", "2012", "pcm"),
                       fitMethod = c("stan", "jags"), 
                       outFormat = NULL,
                       startSmall = TRUE,
                       M = 1000,
                       warmup = 1000,
                       n.chains = 2,
                       thin = 1,
                       # mail = NULL,
                       method = "parallel",
                       return_defaults = TRUE,
                       add2varlist = NULL,
                       cores = NULL,
                       summarise = FALSE,
                       N2 = 2, 
                       ...){
    
    checkmate::assert_matrix(X, mode = "integerish", any.missing = FALSE,
                             min.rows = 2, min.cols = 2)
    N <- nrow(X)
    J <- ncol(X)
    checkmate::assert_integerish(X, lower = 1, upper = 5, any.missing = FALSE)
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 len = J)
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 len = J)
    checkmate::assert_number(df, lower = 0, null.ok = TRUE)
    checkmate::assert_matrix(V, mode = "numeric", any.missing = FALSE,
                             min.rows = 2, min.cols = 2, null.ok = TRUE)
    fitModel <- match.arg(fitModel)
    fitMethod <- match.arg(fitMethod)
    checkmate::qassert(M, "X1")
    checkmate::qassert(warmup, "X1")
    checkmate::qassert(n.chains, "X1")
    checkmate::qassert(thin, "X1")
    checkmate::assert_character(add2varlist, any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_int(cores, lower = 1, null.ok = TRUE)
    checkmate::assert_int(N2, lower = 0, upper = N)
    
    
    n.trait <- length(table(traitItem))
    if(max(traitItem) != n.trait){
        stop("Check definition of traitItem")
    }
    if(length(unique(revItem)) > 1 & !any(cor(X) < 0)){
        stop("All items have positive bivariate correlations;",
             "however, data should be provided in 'input'-format such that",
             "reverse-coded and regular items correlate negatively.")
    }
    if (any(grepl("X_pred", add2varlist)) + (N2 > 2) == 1) {
        stop("If you want to draw posterior predictives for all persons,",
             "add 'X_pred' to argument 'add2varlist' and change argument 'N2'.")
    }

    # number of response processes/dimensions
    dimen <- switch(fitModel, 
                "ext" = 3,  
                # "ext2" = 3,
                # "ext3" = 4,
                # "ext4" = 3,
                # "ext5" = 3,
                "2012" = 2,
                "pcm" = 0) + 1
    S <- dimen - 1 + n.trait
    # if (!is.null(model2) & model2 == "HH") {
    #     S <- S + 1
    # }
    if(is.null(df)){
        df <- S+1
    }
    if(is.null(V)){
        V <- diag(S)
    }
    
    # adjust starting values
    inits <- list()
    if(startSmall == TRUE){
        inits <- lapply(1:n.chains, function(id) {
            get_inits(N = N, J = J, S = S, n.trait = n.trait,
                      fitMethod = fitMethod, fitModel = fitModel)
            })
    } else if (fitMethod == "stan") {
        inits <- vector("list", length = n.chains)
        for (iii in 1:n.chains) {
            inits[[iii]] <- list(# xi_beta  = runif(S, 3/4, 4/3),
                                 mu_beta = truncnorm::rtruncnorm(S, mean = 0, sd = 1,
                                                                 a = -2, b = 2),
                                 beta_raw = matrix(truncnorm::rtruncnorm(J*dimen, mean = 0, sd = 1,
                                                                         a = -2, b = 2),
                                                   nrow = J, ncol = dimen),
                                 beta_ARS_extreme = truncnorm::rtruncnorm(1, mean = 0, sd = 1,
                                                                          a = -2, b = 2),
                                 xi_theta = runif(S, 3/4, 4/3))
            if (fitModel == "pcm") {
                inits[[iii]]$beta_raw = matrix(truncnorm::rtruncnorm(J*4, mean = 0, sd = 1,
                                                                     a = -2, b = 2),
                                               nrow = J, ncol = 4)
                inits[[iii]]$beta_ARS_extreme <- NULL
            } else if (fitModel == "ext5") {
                inits[[iii]]$beta_ARS_extreme <- NULL
                inits[[iii]]$mu_beta <- truncnorm::rtruncnorm(S+1, mean = 0, sd = 1,
                                                             a = -2, b = 2)
                inits[[iii]]$beta_raw <- matrix(truncnorm::rtruncnorm(J*(dimen+1), mean = 0, sd = 1,
                                                                      a = -2, b = 2),
                                                nrow = J, ncol = dimen+1)
            }
        }
    }
    
    datalist <- list(X=X, S=S, df=df, V=V, J=J, N=N, theta_mu=rep(0,S), 
                     revItem=revItem, traitItem=traitItem, N2 = N2)
    varlist <-  c("theta", "beta", "Sigma", "sigma_beta", "mu_beta")
    if(fitModel %in% c("ext", "ext2", "ext3")){
        varlist <- c(varlist, "beta_ARS_extreme")
    }
    # if (fitMethod == "stan") {
    #     datalist$T1_CONST <- T1_CONST
    #     if (fitModel != "2012") {
    #         varlist <- c(varlist, "p_T1_item", "T1_item_obs", "T1_item_pred","p_T1_part", "T1_part_obs", "T1_part_pred")
    #     }
    # } else {
    #     varlist <- c(varlist,  "T_obs","T_pred","post_p"
    #                  # , "deviance", "pd", "popt", "dic"
    #                  # , "pd.i"
    #     )
    # }
    if (!is.null(add2varlist)) varlist <- c(varlist, add2varlist)
    
    time1 <- Sys.time()
    
    if(fitMethod == "jags"){
        ##################### fit JAGS ###############################################
        boeck.jags <- runjags::run.jags(model=paste0(.libPaths()[1],"/mpt2irt/models/jags_boeck_",
                                                     fitModel, ".txt"), 
                                        monitor = varlist, 
                                        data = datalist, inits = inits,
                                        n.chains = n.chains, sample = M, 
                                        burnin = ceiling(warmup*4/5), adapt = ceiling(warmup/5),
                                        thin = thin, method = method, 
                                        summarise = summarise, modules = c("glm", "dic"), ...)
        if (is.null(outFormat)) {
            boeck.samp <- boeck.jags
        } else if (outFormat == "mcmc.list") {
            try(boeck.samp <- runjags:::as.mcmc.list.runjags(boeck.samp))
        } else if (outFormat == "stan") {
            try(boeck.samp <- runjags:::as.mcmc.list.runjags(boeck.samp))
            try(boeck.samp <- mcmc.list2stan(boeck.samp))
        }
        time2 <- Sys.time()
        boeck.samp$date <- c(strftime(time1), strftime(time2), paste(round(difftime(time2, time1, units = "hours"), 2), "hours"))
        names(boeck.samp$date) <- c("Start", "End", "Difference")

    } else {
        ##################### fit Stan ###############################################
        
        # library("rstan")
        # loadNamespace("rstan")
        
        # if (!is.null(model2) & model2 == "HH") {
        #     # fitting "model" requires changing S2 to S2+1 (see above)
        #     stanExe <- boeck_stan_ext_HH
        # } else
        if (fitModel == "ext"){
            stanExe <- boeck_stan_ext
        } else if (fitModel == "pcm"){
            stanExe <- boeck_stan_pcm
        } else if (fitModel == "ext2"){
            stanExe <- boeck_stan_ext2
        } else if (fitModel == "ext3"){
            stanExe <- boeck_stan_ext3
        } else if (fitModel == "ext4"){
            stanExe <- boeck_stan_ext4
        } else if (fitModel == "ext5"){
            stanExe <- boeck_stan_ext5
        } else {
            stanExe <- boeck_stan_2012
        }
        if (is.null(cores)) cores <- min(n.chains, parallel::detectCores() - 1)
        boeck.samp <- rstan::sampling(stanExe,
                                      data = datalist, pars = varlist, 
                                      chains = n.chains, iter = warmup+M*thin, 
                                      warmup=warmup, thin = thin,
                                      cores = cores,
                                      init = inits, ...)
        try( unlink(paste0(getwd(), "\\_StanProgress.txt")), silent=T)
        if (is.null(outFormat)) {
            boeck.samp <- boeck.samp
        } else if (outFormat == "mcmc.list"){
            try(boeck.samp <- rstan::As.mcmc.list(boeck.samp))
        }
        time2 <- Sys.time()
        boeck.samp@date <- c(boeck.samp@date, strftime(time1), strftime(time2), paste(round(difftime(time2, time1, units = "hours"), 2), "hours"))
        names(boeck.samp@date) <- c("Stan", "Start", "End", "Difference")
    }
    # if(!missing(mail) && !is.null(mail))
    #     mail::sendmail(mail, subject=paste("IRT-MPT Model finished on", Sys.info()["nodename"]),
    #                    message="Calculation finished!", password="rmail")
    
    if (return_defaults == FALSE) {
        return(boeck.samp)
    } else {
        return(list(samples = boeck.samp, V = V, df = df, startSmall = startSmall,
                    fitModel = fitModel, fitMethod = fitMethod))
    }
    
}

# @export
get_inits <- function(N, J, S, n.trait, fitMethod, fitModel){
  if (fitMethod == "jags") {
    ll <- list(beta_raw = matrix(rnorm( J*(S-n.trait+1),0, 1), J),
               # xi_beta = runif(S, 3/4, 4/3),
               Tau.theta.raw=rWishart(1,S+3,diag(S))[,,1],
               xi_theta = runif(S, 3/4, 4/3))
  } else {
    ll <- list(mu_beta         = rnorm(S, 0, .5)
               , beta_raw        = matrix(rnorm(J*min(4, S-n.trait+1), 0, 1), J)
               , theta_raw       = matrix(rnorm(N*S, 0, 1), N)
               , Sigma_raw       = solve(rWishart(1, S+3, diag(S))[, , 1])
               , sigma2_beta_raw = runif(S, 3/4, 4/3)
               # , xi_beta         = runif(S, 3/4, 4/3)
               , xi_theta        = runif(S, 3/4, 4/3)
               )
    if (fitModel == "pcm") {
        ll$beta_raw <- matrix(rnorm(J*4, 0, 1), J, 4)
    }
  }
  return(ll)
}