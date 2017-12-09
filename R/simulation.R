#' Recovery simulation of mpt2irt models.
#' 
#' This function allows to run a simulation study of mpt2irt models. Data are
#' generated either from the Boeckenholt Model (\code{genModel = "2012"}) or
#' from the Acquiescence Model (\code{genModel = "ext"}). Subsequently, one or
#' both of these models are fit to the generated data using either JAGS or Stan.
#' The results are saved in an RData file in \code{dir}.
#' 
#' @param J number of items. Can be a vector for multiple traits (e.g.,
#'   J=c(10,10,10)).
#' @param prop.rev number of reversed items. Can be a vector for multiple
#'   traits(e.g., prop.rev=c(5,3,5)/10)
#' @param rrr Sequence of integers (e.g., \code{1:100}) of length greater or
#'   equal to 1 specifying the number of replications to run.
#' @param genModel Character. The data generating model (either "2012" or "ext").
#' @param fitModel Character. The model for data analysis ("2012", "ext", or both as vector
#'   c("2012", "ext")).
#' @param fitMethod Character. Whether to use "stan" or "jags".
# @param trait.sd standard deviation of normally distributed item properties
#   for trait. Can be a vector for multiple traits (e.g., trait.sd=c(2, 1.5,
#   2))
# @param RSmin minimum of range for uniformly distributed response styles (can
#   be a vector with separate bounds for middle/extreme/(acquiesence))
# @param RSmax maximum of range for response styles (can also be a vector)
#' @param theta_vcov true covariance matrix of response processes (order:
#'   middle, extreme, (acquiescence), trait). standard is diag(3) / diag(4). Can be a vector of variances (not SDs).
#' @param betas Optional list. May have entries \code{"beta.mrs"},
#'   \code{"beta.ers"}, \code{"beta.trait"}, and/or \code{"beta.ars"}. Each of
#'   those may have arguments passed to \code{\link[truncnorm]{rtruncnorm}}.
#' @param df_vcov Numeric. Degrees of freedom for wishart distribution from
#'   which the variance-covariance matrix for generating the data is drawn.
#' @param dir Path to directory where results should be stored,
#' @param keep_mcmc Logical indicating wheter to retain, besides a summary of the parameters, the raw mcmc samples.
#' @param savext_mcmc Logical indicating wheter to save the mcmc samples in
#'   an external RData file.
#' @param savext_all Logical indicating wheter to save the output from Stan/JAGS
#'   in an external RData file.
#' @param beta_ARS_extreme Numeric. Only for \code{genModel="ext"}: probability
#'   (on probit scale) of choosing category 5 (vs.4) in case of ARS. Defaults to
#'   \code{rtruncnorm(mean = qnorm(.7), sd = sqrt(.1), a = qnorm(.5), b =
#'   qnorm(.9))}.
# @param M number of MCMC samples (adaptation + burnin = M/2)
# @param n.chains number of chains
# @param thin thinning of MCMC samples
# @param printProgress whether to print the progress in a text file (either at the location defined by \code{path} or at the working directory)
# @param path where to save progress report (e.g., path="C:/"). disabled by using path=""
# @param saveTemp save temporary results to same location as progress report
# @param rng_seed random number seed
# @param mail an email address used when the simulation is finished (no dash "-" allowed!)
# @param ... further arguments passed to \code{\link[rstan]{sampling}} (for Stan) or \code{\link[runjags]{run.jags}} (for JAGS)
#' @inheritParams fit_irtree
#' @inheritParams generate_irtree_ext
#' @inheritParams runjags::run.jags
#' @return Function does not directly return anything but saves an external
#'   RData file to \code{dir}. This object is a list containing the generated
#'   parameters in \code{sim-results$param.sum$gen}, fitted parameters and other model fit
#'   information in \code{sim-results$param.sum$foo}, as well as a summary of the setup.
#' @details Note that a text file "progress.txt" is written (and updated) to \code{dir} informing you about the progress of the simulation.
#' @import coda
# @import runjags
# @import parallel
#' @examples
#' \dontrun{
#' recovery_irtree(rrr = 1:2, N = 20, J = 10, genModel = "ext", fitModel = "ext",
#'                 fitMethod = "stan", M = 200, n.chains = 2, warmup = 200,
#'                 dir = "~/")
#'                 
#' # run multiple simulations in parallel using the 'parallel' package
#' no_cores <- parallel::detectCores() - 1
#' cl <- parallel::makeCluster(no_cores)
#' parallel::clusterApplyLB(cl, x = 11:13, fun = recovery_irtree, cores = 1,
#'                          N = 20, J = 10, genModel = "ext", fitModel = "ext",
#'                          fitMethod = "stan", M = 200, n.chains = 2, warmup = 200,
#'                          dir = "~/")
#' parallel::stopCluster(cl = cl)
#'}
#' @export
recovery_irtree <- function(rrr = NULL,
                            N = NULL,
                            J = NULL,
                            prop.rev = .5,
                            genModel = c("ext", "2012"),
                            fitModel = c("ext", "2012", "pcm", "steps", "shift", "ext2"),
                            fitMethod = c("stan", "jags"), 
                            theta_vcov = NULL,
                            betas = NULL,
                            beta_ARS_extreme = NULL,
                            df = NULL,
                            V = NULL,
                            M = 500,
                            n.chains = 3,
                            thin = 1,
                            warmup = 500,
                            method = "simple",
                            outFormat = NULL,
                            startSmall = FALSE,
                            df_vcov = 50,
                            dir = NULL,
                            keep_mcmc = FALSE,
                            savext_all = FALSE,
                            savext_mcmc = TRUE,
                            add2varlist = c("deviance", "pd", "popt", "dic"), ...) {
    
    # It is assumed that parameters such as df and V are the same for all models
    # in fitModel. There's no problem in changing this, only remember to change
    # the output, e.g., for V from length 1 to length(fitModel).
    
    checkmate::qassert(rrr, "X>0[1,]")
    checkmate::qassert(N, "X1")
    checkmate::qassert(J, "X>0[1,]")
    genModel <- match.arg(genModel)
    fitMethod <- match.arg(fitMethod)
    if (fitMethod == "stan") {
        fitModel <- match.arg(fitModel, several.ok = TRUE)
    } else {
        # PCM and Steps not yet implemented in JAGS
        fitModel <- match.arg(fitModel, choices = c("ext", "2012"),
                              several.ok = TRUE)
    }
    checkmate::qassert(df_vcov, "N1[1,]")
    
    args <- c(as.list(environment()), list(...))
    
    ### INPUT WRANGLING -----
    
    #### multiple traits
    n.trait <- length(J)
    traitItem <- rep(1:n.trait, J)
    J <- sum(J)
    
    numRS.gen <- ifelse(genModel == "2012", 2, 3)
    S.gen <- numRS.gen + n.trait
    numRS.fit <- sapply(fitModel, function(x) ifelse(x == "2012", 2, 3))
    S.fit <- numRS.fit + n.trait
    
    ### DIRECTORIES -----
    
    # if (!is.null(dir)) {
    #     try(dir <- path.expand(dir))
    #     if (!dir.exists(dir)) {
    #         set.seed(Sys.Date())
    #         tmp1 <- paste0(sample(c(letters, LETTERS, 0:9), 12, T), collapse = "")
    #         dirx <- paste0(getwd(), "/", tmp1)
    #         if (!dir.exists(dirx)) {
    #             dir.create(dirx)
    #         }
    #         on.exit(message(paste0("Data saved in: ", dirx)))
    #     }
    # }
    # 
    # if (savext_mcmc == TRUE) {
    #     if (!dir.exists(paste0(ifelse(dir.exists(dir), dir, dirx), "/", "mcmc"))) {
    #         dir.create(paste0(ifelse(dir.exists(dir), dir, dirx), "/", "mcmc"))
    #     }
    # }
    
    ### LOOP OVER REPLICATIONS -----
    
    for (qqq in seq_along(rrr)) {
        
        time_1 <- Sys.time()
        
        ### GENERATE DATA -----
        cor2cov <- function(mat = NULL, sd = NULL) {
            # is in MBESS package, but MBESS has so many dependencies
            checkmate::qassert(mat, "M+[-1,1]")
            checkmate::qassert(sd, "N+(0,1]")
            return(diag(sd) %*% mat %*% diag(sd))
        }
        
        if (is.null(theta_vcov)) {
            xx1 <- diag(S.gen)
            xx1[1, 2] <- xx1[2, 1] <- -.2
            theta_vcov_i <- cov2cor(rWishart(1, df = df_vcov, xx1)[, , 1])
            # the default variances are 0.33 for RS and 1.0 for TRAIT
            theta_vcov_i <- cor2cov(theta_vcov_i, sqrt(c(rep(.33, S.gen - n.trait), rep(1, n.trait))))
            # theta_vcov_i <- diag(S.gen)  # middle  / extremity / acq / trait(s)
        } else if (is.vector(theta_vcov)) {
            repeat {
                vcov.0 <- cov2cor(rWishart(1, df = df_vcov, diag(theta_vcov))[, , 1])
                theta_vcov_i <- cor2cov(vcov.0, sqrt(theta_vcov))
                # theta_vcov <- diag(S.gen) * theta_vcov   # middle /  extremity / trait(s)
                if (det(theta_vcov_i) != 0) break
            }
        } else {
            theta_vcov_i <- theta_vcov
        }
        
        if (any(dim(theta_vcov_i)!= S.gen) | det(theta_vcov_i) == 0){
            warning("Check definition of theta_vcov!")
        }
        
        
        if (genModel == "2012") {
            betas_i <- gen_betas(genModel = genModel, J = J, betas = betas)
            gen <- generate_irtree_2012(N = N, J = J, betas = betas_i, theta_vcov = theta_vcov_i,
                                  prop.rev = prop.rev, traitItem = traitItem, cat = TRUE)
        } else if (genModel == "ext") {
            betas_i <- gen_betas(genModel = genModel, J = J, betas = betas)
            
            if (is.null(beta_ARS_extreme)) {
                beta_ARS_extreme_i <- truncnorm::rtruncnorm(n = 1, mean = qnorm(.7), sd = sqrt(.1),
                                                          a = qnorm(.5), b = qnorm(.9))
            } else {
                beta_ARS_extreme_i <- beta_ARS_extreme
            }
            gen <- generate_irtree_ext(N = N, J = J, betas = betas_i, theta_vcov = theta_vcov_i,
                                 prop.rev = prop.rev, traitItem = traitItem, cat = TRUE,
                                 beta_ARS_extreme = beta_ARS_extreme_i, genModel = genModel)
        }
        
        genNames <- c(paste0("theta[",apply(expand.grid(1:N, 1:S.gen),1, 
                                            paste0, collapse=","),"]"),
                      paste0("beta[",apply(expand.grid(1:J, 1:(numRS.gen+1)),1, 
                                           paste0, collapse=","),"]"), 
                      paste0("Sigma[",apply(expand.grid(1:S.gen, 1:S.gen),1, 
                                            paste0, collapse=","),"]"))
        
        
        if (genModel != "2012") {
            genNames <- c(genNames, "beta_ARS_extreme")
        }
        param.sum <- vector("list", length(S.fit) + 1)
        names(param.sum) <- c(fitModel, "gen")
        param.sum$gen <- matrix(c(gen$theta, gen$betas, gen$theta_vcov,
                                  switch(genModel, "ext" = gen$beta_ARS_extreme)),
                                ncol = 1, dimnames = list(genNames, NULL))
        
        returnlist <- list(param.sum = param.sum,
                           # genModel = genModel,
                           sim_args = args
                           # fitModel = fitModel,
                           # S = S.fit,
                           # revItem = gen$revItem,
                           # traitItem = gen$traitItem,
                           # X = gen$X,
                           # N = N,
                           # J = J,
                           # prop.rev = prop.rev,
                           # n.trait = n.trait,
                           # fitMethod = fitMethod,
                           # df = NULL,
                           # V = NULL,
                           # M = M,
                           # n.chains = n.chains,
                           # thin = thin,
                           # warmup = warmup,
                           # df_vcov = df_vcov,
                           # startSmall = startSmall,
                           # jagspath = runjags::runjags.getOption("jagspath")
                           )
        
        ### FIT MODELS -----
        
        for (sss in 1:length(S.fit)) {
            
            ### FIT IN JAGS -----
            
            if (fitMethod == "jags") {
                fit_jags <- fit_irtree(X = gen$X,
                                       revItem = gen$revItem,
                                       traitItem = gen$traitItem,
                                       fitModel = fitModel[sss], fitMethod = fitMethod, df = df, V = V,
                                       n.chains = n.chains, thin = thin, M = M, warmup = warmup,
                                       outFormat = outFormat, startSmall = startSmall,
                                       method = method,
                                       add2varlist = add2varlist, ...)
                
                # returnlist$sim_args$jagspath <- runjags::runjags.getOption("jagspath")
                
                fit2 <- summarize_irtree_fit(fit_jags)
                fit3 <- tidyup_irtree_fit(fit2, plot = FALSE)
                
                returnlist$param.sum[[fitModel[sss]]] <- fit3
                
                tmp0 <- sprintf("sim_%04d_%s_%s_mcmc.RData", rrr[qqq], fitMethod, fitModel[sss])
                
                dirx <- get_dir(dir = dir)
                dirmcmc <- paste0(dirx, "/mcmc")
                
                if (any(c(savext_mcmc, savext_all) == TRUE)) {
                    if (!dir.exists(dirmcmc)) {
                        dir.create(normalizePath(dirmcmc))
                    }
                }
                
                if (savext_mcmc == TRUE) {
                    
                    tmp11 <- gsub("_mcmc.RData", "", tmp0) %>% 
                        gsub("sim", "mcmc", x = .)
                    
                    assign(tmp11, fit2$mcmc)
                    
                    tmp13 <- paste0(dirmcmc, "/", tmp0)
                    do.call(save, list(tmp11, file = tmp13))
                }
                if (savext_all == TRUE) {
                    
                    tmp21 <- gsub("_mcmc.RData", "", tmp0) %>% 
                        gsub("sim", "res", x = .)
                    
                    tmp22 <- gsub("mcmc", "raw", tmp0)
                    tmp23 <- paste0(dirmcmc, "/", tmp22)
                    
                    assign(tmp21, fit_jags)
                    
                    do.call(save, list(tmp21, file = tmp23))
                }
                
                rm(fit_jags, fit2, fit3)
                
                # sum.jags <- runjags::add.summary(fit_jags$samples)
                # # sum.jags$mcmc <- NULL
                # if (savext_mcmc == TRUE) {
                #     returnlist_mcmc <- sum.jags$mcmc
                # }
                # if (keep_mcmc == FALSE) {
                #     for (iii in seq_along(sum.jags$mcmc)) {
                #         sum.jags$mcmc[[iii]] <- sum.jags$mcmc[[iii]][, 1, drop = F]
                #     }
                #     sum.jags$crosscorr <- NULL
                # }
                # returnlist$param.sum[[fitModel[sss]]] <- sum.jags
                # rm(fit_jags, sum.jags)
                
            } else if (fitMethod == "stan") {
                
                ### FIT IN STAN ----
                
                fit_stan <- fit_irtree(X = gen$X,
                                      revItem = gen$revItem,
                                      traitItem = gen$traitItem,
                                      fitModel = fitModel[sss], fitMethod = fitMethod, df = df, V = V,
                                      n.chains = n.chains, thin = thin, M = M, warmup = warmup,
                                      outFormat = outFormat, startSmall = startSmall,
                                      add2varlist = NULL, ...)
                
                
                fit2 <- summarize_irtree_fit(fit_stan)
                fit3 <- tidyup_irtree_fit(fit2, plot = FALSE)
                
                returnlist$param.sum[[fitModel[sss]]] <- fit3
                
                tmp0 <- sprintf("sim_%04d_%s_%s_mcmc.RData", rrr[qqq], fitMethod, fitModel[sss])
                
                dirx <- get_dir(dir = dir)
                dirmcmc <- paste0(dirx, "/mcmc")
                
                if (any(c(savext_mcmc, savext_all) == TRUE)) {
                    if (!dir.exists(dirmcmc)) {
                        dir.create(normalizePath(dirmcmc))
                    }
                }
                
                if (savext_mcmc == TRUE) {
                    
                    tmp11 <- gsub("_mcmc.RData", "", tmp0) %>% 
                        gsub("sim", "mcmc", x = .)
                    
                    assign(tmp11, fit2$mcmc)
                    
                    tmp13 <- paste0(dirmcmc, "/", tmp0)
                    do.call(save, list(tmp11, file = tmp13))
                }
                if (savext_all == TRUE) {
                    
                    tmp21 <- gsub("_mcmc.RData", "", tmp0) %>% 
                        gsub("sim", "res", x = .)
                    
                    tmp22 <- gsub("mcmc", "raw", tmp0)
                    tmp23 <- paste0(dirmcmc, "/", tmp22)
                    
                    assign(tmp21, fit_stan)
                    
                    do.call(save, list(tmp21, file = tmp23))
                }
                
                rm(fit_stan, fit2, fit3)
                
                # returnlist$param.sum[[fitModel[sss]]] <- fit_stan
                # returnlist$param.sum[[fitModel[sss]]]$mcmc <- rstan::As.mcmc.list(fit_stan$samples)
                # returnlist$param.sum[[fitModel[sss]]]$summary <-
                #     summary(returnlist$param.sum[[fitModel[sss]]]$mcmc)
                # if (savext_mcmc == TRUE) {
                #     returnlist_mcmc <- returnlist$param.sum[[fitModel[sss]]]$mcmc
                # }
                # if (keep_mcmc == FALSE) {
                #     returnlist$param.sum[[fitModel[sss]]]$mcmc <- NULL
                #     returnlist$param.sum[[fitModel[sss]]]$samples@sim$samples <- vector("list", length = n.chains)
                # }
                # rm(fit_stan)
            }
        }
        
        ### RETURN -----
        
        time_2 <- Sys.time()
        returnlist$date <- difftime(time_2, time_1, units = "hours")
        
        save_file <- paste0("sim-results-", sprintf("%04d", rrr[qqq]), "-", substr(Sys.time(), 1, 10), ".RData")
        save_file <- gsub(" CET", "", save_file)
        save_file <- gsub(" ", "-", save_file)
        save_file <- gsub(":", "-", save_file)
        
        returnName <- sprintf("sim_%04d", rrr[qqq])
        assign(returnName, returnlist)
        
        dirx <- get_dir(dir = dir)
        
        do.call(save, list(returnName, file = paste0(dirx, "/", save_file)))
        
        # if (savext_mcmc == TRUE) {
        #     save_file_mcmc <- paste0("mcmc-", sprintf("%04d", rrr[qqq]), ".RData")
        #     return_name_mcmc <- paste0("mcmc_", sprintf("%04d", rrr[qqq]))
        #     assign(return_name_mcmc, returnlist_mcmc)
        #     if (!dir.exists(paste0(ifelse(dir.exists(dir), dir, dirx), "/", "mcmc"))) {
        #         dir.create(paste0(ifelse(dir.exists(dir), dir, dirx), "/", "mcmc"))
        #     }
        #     do.call(save, list(return_name_mcmc, file = paste0(ifelse(dir.exists(dir), dir, dirx), "/mcmc/", save_file_mcmc)))
        # }
        
        if (!is.null(dir)) {
            if (dir.exists(dir)) {
                # cony <- file(paste0(ifelse(dir.exists(dir), dir, dirx), "/progress_", Sys.info()["nodename"], ".txt"), open = "a")
                cony <- file(paste0(dir, "/progress.txt"), open = "a")
                writeLines(con = cony,
                           text = paste0(sprintf("%04d", rrr[qqq]),
                                         "\t\tStarted at: ", time_1,
                                         "\t\tEnded at: ", time_2,
                                         "\t\tDifference of: ", format(round(difftime(time_2, time_1, units = "hours"), 3), nsmall = 3, width = 6), # " hours"
                                         "\t\t", Sys.info()["nodename"]
                           )
                )
                close(cony)
            }
        }
    }
}
