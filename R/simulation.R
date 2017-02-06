#' Recovery simulation of mpt2irt models.
#' 
#' This function allows to run a simulation study of mpt2irt models. Data are
#' generated either from the Boeckenholt Model (\code{genModel = "2012"}) or
#' from the Acquiescence Model (\code{genModel = "ext"}). Subsequently, one or
#' both of these models are fit to the generated data using eihter JAGS or Stan.
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
#' @param theta.vcov true covariance matrix of response processes (order:
#'   middle, extreme, (acquiescence), trait). standard is diag(3) / diag(4). Can be a vector of variances (not SDs).
#' @param betas Optional list. May have entries \code{"beta.mrs"},
#'   \code{"beta.ers"}, \code{"beta.trait"}, and/or \code{"beta.ars"}. Each of
#'   those may have arguments passed to \code{\link[truncnorm]{rtruncnorm}}.
#' @param df_vcov Numeric. Degrees of freedom for wishart distribution from
#'   which the variance-covariance matrix for generating the data is drawn.
#' @param dir Path to directory where results should be stored,
#' @param keep_mcmc Logical indicating wheter to retain, besides a summary of the parameters, the raw mcmc samples.
#' @param savext_mcmc Logical indicating wheter to save the raw mcmc samples in an external RData file.
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
#' @inheritParams runjags::run.jags
#' @return Function does not directly return anything but saves an external
#'   RData file to \code{dir}. This object is a list containing the generated
#'   parameters in \code{sim-results$param.sum$gen}, fitted parameters and other model fit
#'   information in \code{sim-results$param.sum$foo}, as well as a summary of the setup.
#' @details Note that a text file "progress.txt" is written (and updated) to \code{dir} informing you about the progress of the simulation.
# @importFrom coda
# @import doParallel
# @import runjags
# @import parallel
# @import foreach
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
                            fitModel = NULL,
                            fitMethod = c("stan", "jags"), 
                            theta.vcov = NULL,
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
                            savext_mcmc = FALSE,
                            add2varlist = c("deviance", "pd", "popt", "dic"),
                            ...) {
    
    # It is assumed that parameters such as df and V are the same for all models
    # in fitModel. There's no problem in changing this, only remember to change
    # the output, e.g., for V from length 1 to length(fitModel).
    
    checkmate::qassert(rrr, "X>0[1,]")
    checkmate::qassert(N, "X1")
    checkmate::qassert(J, "X>0[1,]")
    # checkmate::qassert(prop.rev, "X1[0,1]")
    genModel <- match.arg(genModel)
    checkmate::qassert(fitModel, "S>0")
    fitMethod <- match.arg(fitMethod)
    # checkmate::assert_list(betas, null.ok = TRUE)
    # checkmate::assert_number(beta_ARS_extreme, finite = TRUE, null.ok = TRUE)
    checkmate::qassert(df_vcov, "N1[1,]")
    
    ### INPUT WRANGLING --------------------------------------------------------
    
    #### multiple traits
    n.trait <- length(J)
    traitItem <- rep(1:n.trait, J)
    if(length(prop.rev) == 1)
        prop.rev <- rep(prop.rev, n.trait)
    #     if(length(trait.sd) ==1){
    #         trait.sd <- rep(trait.sd,  sum(J))
    #     }else if(length(trait.sd) == n.trait){
    #         trait.sd <- rep(trait.sd, J)
    #     }
    J <- sum(J)
    
    genModel <- as.character(genModel)
    numRS.gen <- ifelse(genModel == "2012", 2, 3)
    S.gen <- numRS.gen + n.trait
    numRS.fit <- sapply(fitModel, function(x) ifelse(x == "2012", 2, 3))
    S.fit <- numRS.fit + n.trait
    
    for (qqq in seq_along(rrr)) {
        
        time_1 <- Sys.time()
        
        ### GENERATE A SINGLE DATASET ----------------------------------------------
        
        if (missing(theta.vcov)) {
            xx1 <- diag(S.gen)
            xx1[1, 2] <- xx1[2, 1] <- -.2
            theta_vcov_i <- cov2cor(rWishart(1, df = df_vcov, xx1)[, , 1])
            # the default variances are 0.33 for RS and 1.0 for TRAIT
            theta_vcov_i <- MBESS::cor2cov(theta_vcov_i, sqrt(c(rep(.33, S.gen - n.trait), rep(1, n.trait))))
            # theta_vcov_i <- diag(S.gen)  # middle  / extremity / acq / trait(s)
        } else if (is.vector(theta.vcov)) {
            repeat {
                vcov.0 <- cov2cor(rWishart(1, df = df_vcov, diag(theta.vcov))[, , 1])
                theta_vcov_i <- MBESS::cor2cov(vcov.0, sqrt(theta.vcov))
                # theta.vcov <- diag(S.gen) * theta.vcov   # middle /  extremity / trait(s)
                if (det(theta_vcov_i) != 0) break
            }
        } else {
            theta_vcov_i <- theta_vcov
        }
        
        if (any(dim(theta_vcov_i)!= S.gen) | det(theta_vcov_i) == 0){
            warning("Check definition of theta.vcov!")
        }
        
        
        if (genModel == "2012") {
            # if (missing(betas)) {
            #     beta.mrs <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
            #                                       a = qnorm(.5), b = qnorm(.9))
            #     beta.ers <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
            #                                       a = qnorm(.5), b = qnorm(.9))
            #     beta.trait <- truncnorm::rtruncnorm(n = J, mean = 0, sd = sqrt(.5),
            #                                         a = qnorm(.3), b = qnorm(.7))
            #     betas_i <- cbind(beta.mrs, beta.ers, beta.trait)
            # } else {
            #     betas_i <- betas
            # }
            betas_i <- gen_betas(genModel = genModel, J = J, betas = betas)
            gen <- generate_irtree_2012(N = N, J = J, betas = betas_i, theta.vcov = theta_vcov_i,
                                  prop.rev = prop.rev, traitItem = traitItem, cat = TRUE)
        } else if (genModel == "ext") {
            # if (missing(betas)) {
            #     beta.mrs <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
            #                                       a = qnorm(.5), b = qnorm(.9))
            #     beta.ers <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
            #                                       a = qnorm(.5), b = qnorm(.9))
            #     beta.ars <- truncnorm::rtruncnorm(n = J, mean = qnorm(.95), sd = sqrt(.1),
            #                                       a = qnorm(.8), b = qnorm(.999))
            #     # beta.ars <- truncnorm::rtruncnorm(n = J, mean = qnorm(.9), sd = sqrt(.1),
            #     #                                   a = qnorm(.8), b = qnorm(.99))
            #     beta.trait <- truncnorm::rtruncnorm(n = J, mean = 0, sd = sqrt(.5),
            #                                         a = qnorm(.3), b = qnorm(.7))
            #     betas_i <- cbind(beta.mrs, beta.ers, beta.ars, beta.trait)
            # } else {
            #     betas_i <- betas
            # }
            betas_i <- gen_betas(genModel = genModel, J = J, betas = betas)
            
            if (is.null(beta_ARS_extreme)) {
                beta_ARS_extreme_i <- truncnorm::rtruncnorm(n = 1, mean = qnorm(.7), sd = sqrt(.1),
                                                          a = qnorm(.5), b = qnorm(.9))
            } else {
                beta_ARS_extreme_i <- beta_ARS_extreme
            }
            gen <- generate_irtree_ext(N = N, J = J, betas = betas_i, theta.vcov = theta_vcov_i,
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
        param.sum$gen <- matrix(c(gen$theta, gen$betas, gen$theta.vcov,
                                  switch(genModel, "ext" = gen$beta_ARS_extreme)),
                                ncol = 1, dimnames = list(genNames, NULL))
        
        ### FIT THAT SHIT IN JAGS ----------------------------------------------
        
        returnlist <- list(param.sum = param.sum, genModel = genModel, fitModel = fitModel,
                           S = S.fit, revItem = gen$revItem, traitItem = gen$traitItem,
                           X = gen$X, N = N, J = J, prop.rev = prop.rev, n.trait = n.trait,
                           fitMethod = fitMethod, df = NULL, V = NULL,
                           M = M, n.chains = n.chains, thin = thin,
                           warmup = warmup, df_vcov = df_vcov, startSmall = startSmall,
                           jagspath = runjags::runjags.getOption("jagspath"),
                           session = NULL)
        
        # fitpar <- list() ; dic <- list()
        for(sss in 1:length(S.fit)){
            
            # fitNames <- theta_mu <- V.ls <- list()
            # fitNames <- list()
            # for(sss in 1:length(S.fit)){
            #         fitNames[[sss]] <- c(paste0("theta[",apply(expand.grid(1:N, 1:S.fit[sss]),1, 
            #                                                    paste0, collapse=","),"]"),
            #                              paste0("beta[",apply(expand.grid(1:J, 1:(numRS.fit[sss]+1)),1, 
            #                                                   paste0, collapse=","),"]"), 
            #                              paste0("Sigma[",apply(expand.grid(1:S.fit[sss], 1:S.fit[sss]),1, 
            #                                                    paste0, collapse=","),"]"),
            #                              ttt)
            # theta_mu[[sss]] <- rep(0,S.fit[sss])
            
            ### set up parameters
            #         if(missing(V)){
            #             V.ls[[sss]] <- diag(S.fit[sss])
            #         }else if(length(S.fit) == 1){
            #             V.ls[[sss]] <- unlist(V)
            #         }else{
            #             V.ls[[sss]] <- V[[sss]]
            #         }
            # }
            #         if(missing(df)){
            #             df <- S.fit+1
            #         }
            
            if(fitMethod == "jags") {
                fit.jags <- fit_irtree(X = gen$X,
                                      revItem = gen$revItem,
                                      traitItem = gen$traitItem,
                                      fitModel = fitModel[sss], fitMethod = fitMethod, df = df, V = V,
                                      n.chains = n.chains, thin = thin, M = M, warmup = warmup,
                                      outFormat = outFormat, startSmall = startSmall,
                                      return_defaults = TRUE, method = method,
                                      add2varlist = add2varlist, ...)
                sum.jags <- runjags::add.summary(fit.jags$samples)
                # sum.jags$mcmc <- NULL
                if (savext_mcmc == TRUE) {
                    returnlist_mcmc <- sum.jags$mcmc
                }
                if (keep_mcmc == FALSE) {
                    for (iii in seq_along(sum.jags$mcmc)) {
                        sum.jags$mcmc[[iii]] <- sum.jags$mcmc[[iii]][, 1, drop = F]
                    }
                    sum.jags$crosscorr <- NULL
                }
                returnlist$param.sum[[fitModel[sss]]] <- sum.jags
                returnlist$df <- fit.jags$df
                returnlist$V <- fit.jags$V
                returnlist$session <- sessionInfo()
                returnlist$nodename <- Sys.info()["nodename"]
                rm(fit.jags, sum.jags)
                
                #             rjags::load.module("glm", quiet=T)
                #             #         boeck.jags <- run.jags(model=paste0( modelPath,"/mpt2irt/models/jags_boeck_",fitModel[sss],"_1d.txt"), 
                #             #                               monitor = c("theta", "beta", "Sigma", "T_obs","T_pred", "post_p",
                #             #                                           'deviance', 'pd', 'pd.i', 'popt', 'dic'), 
                #             #                               data=datalist, n.chains=n.chains,  
                #             #                               burnin = M/4, sample = M, adapt=M/4, 
                #             #                               summarise = T, thin = thin, method="simple")
                #             boeck.jags <- rjags::jags.model(file=paste0( modelPath,"/mpt2irt/models/jags_boeck_",
                #                                                          fitModel[sss],".txt"), 
                #                                             data=datalist,  n.chains=n.chains, n.adapt=M/4)
                #             adapt <- F
                #             while(!adapt){
                #                 adapt <- rjags::adapt(boeck.jags,M/4,end.adaptation = FALSE)
                #             }
                #             rjags::update(boeck.jags, M/4)
                #             boeck.samp <- coda.samples.dic(boeck.jags, 
                #                                            variable.names = c("theta", "beta", "Sigma",
                #                                                               "T_obs","T_pred", "post_p"),
                #                                            n.iter=M*thin, thin = thin)
                #             dic[[sss]] <- boeck.samp$dic
                #             boeck.stan <- mcmc.list2stan(boeck.samp$samples)
                #             rm(boeck.samp)
            } else if (fitMethod == "stan") {
                ### Stan
                fit.stan <- fit_irtree(X = gen$X,
                                      revItem = gen$revItem,
                                      traitItem = gen$traitItem,
                                      fitModel = fitModel[sss], fitMethod = fitMethod, df = df, V = V,
                                      n.chains = n.chains, thin = thin, M = M, warmup = warmup,
                                      outFormat = outFormat, startSmall = startSmall,
                                      return_defaults = TRUE, add2varlist = NULL, ...)
                returnlist$param.sum[[fitModel[sss]]] <- fit.stan
                returnlist$param.sum[[fitModel[sss]]]$mcmc <- rstan::As.mcmc.list(fit.stan$samples)
                returnlist$param.sum[[fitModel[sss]]]$summary <-
                    coda:::summary.mcmc.list(returnlist$param.sum[[fitModel[sss]]]$mcmc)
                returnlist$df <- fit.stan$df
                returnlist$V <- fit.stan$V
                returnlist$session <- sessionInfo()
                returnlist$nodename <- Sys.info()["nodename"]
                if (savext_mcmc == TRUE) {
                    returnlist_mcmc <- returnlist$param.sum[[fitModel[sss]]]$mcmc
                }
                if (keep_mcmc == FALSE) {
                    returnlist$param.sum[[fitModel[sss]]]$mcmc <- NULL
                    returnlist$param.sum[[fitModel[sss]]]$samples@sim$samples <- vector("list", length = n.chains)
                    # for (iii in seq_along(fit.stan$samples@sim$samples)[-1]) {
                    #     returnlist$param.sum[[fitModel[sss]]]$samples@sim$samples[[iii]] <-
                    #         lapply(returnlist$param.sum[[fitModel[sss]]]$samples@sim$samples[[iii]], function(x) x <- NULL)
                    #     attributes(returnlist$param.sum[[fitModel[sss]]]$samples@sim$samples[[iii]]) <-
                    #         lapply(attributes(returnlist$param.sum[[fitModel[sss]]]$samples@sim$samples[[iii]]), function(x) x <- NULL)
                    # }
                }
                rm(fit.stan)
            }
            #         else {
            #             ### Stan
            #             data(boeck_stan_models)
            #             if(fitModel == "ext"){
            #                 stanExe <- boeck_stan_ext
            #             }else{
            #                 stanExe <- boeck_stan_2012
            #             }
            #             boeck.stan <- rstan::sampling(stanExe, 
            #                                           # model_name=paste0("Boeckenholt_", fitModel[sss]),
            #                                           pars=c("theta","beta", "Sigma","T_obs","T_pred","post_p"),
            #                                           data=datalist, chains = n.chains, thin=thin, 
            #                                           iter=M*thin)
            #             dic[[sss]] <- NULL
            #         }
            
            # fitpar[[sss]] <- rstan::monitor(boeck.stan, print=T)[fitNames[[sss]],] 
        }
        
        
        ### RETURN SOMETHING -------------------------------------------------------
        
        # if (is.null(save.dir)) {
        #     exist.dirs <- list.dirs(getwd(), recursive = F, full.names = T)
        #     exist.dir2 <- grepl("sim-results", exist.dirs)
        #     if (length(exist.dirs) > 0) {
        #         save.dir <- tail(exist.dirs[exist.dir2], 1)
        #         save.di2 <- paste0(save.dir, "/", Sys.info()[["nodename"]])
        #         if (!dir.exists(save.di2)) {
        #             dir.create(save.di2)
        #         } else {
        #             save.dir <- paste0(getwd(), "/sim-results-",
        #                                sprintf("%02d", sum(exist.dir2) + 1))
        #             dir.create(save.dir)
        #             save.di2 <- paste0(save.dir, "/", Sys.info()[["nodename"]])
        #             dir.create(save.di2)
        #         }
        #     } else {
        #         save.dir <- paste0(getwd(), "/sim-results-", sprintf("%02d", 1))
        #         dir.create(save.dir)
        #         save.di2 <- paste0(save.dir, "/", Sys.info()[["nodename"]])
        #         dir.create(save.di2)
        #     }
        # }
        time_2 <- Sys.time()
        returnlist$date <- difftime(time_2, time_1, units = "hours")
        try(dir <- path.expand(dir))
        if (!dir.exists(dir)) {
            set.seed(Sys.Date())
            tmp1 <- paste0(sample(c(letters, LETTERS, 0:9), 12, T), collapse = "")
            dirx <- paste0(getwd(), "/", tmp1)
            if (!dir.exists(dirx)) {
                dir.create(dirx)
            }
            on.exit(warning(paste0("Data saved in: ", dirx)))
        }
        save_file <- paste0("sim-results-", sprintf("%04d", rrr[qqq]), "-", substr(Sys.time(), 1, 10), ".RData")
        save_file <- gsub(" CET", "", save_file)
        save_file <- gsub(" ", "-", save_file)
        save_file <- gsub(":", "-", save_file)
        
        returnName <- paste0("sim_sum_", sprintf("%04d", rrr[qqq]))
        assign(returnName, returnlist)
        
        do.call(save, list(returnName, file = paste0(ifelse(dir.exists(dir), dir, dirx), "/", save_file)))
        
        if (savext_mcmc == TRUE) {
            save_file_mcmc <- paste0("mcmc-", sprintf("%04d", rrr[qqq]), ".RData")
            return_name_mcmc <- paste0("mcmc_", sprintf("%04d", rrr[qqq]))
            assign(return_name_mcmc, returnlist_mcmc)
            if (!dir.exists(paste0(ifelse(dir.exists(dir), dir, dirx), "/", "mcmc"))) {
                dir.create(paste0(ifelse(dir.exists(dir), dir, dirx), "/", "mcmc"))
            }
            do.call(save, list(return_name_mcmc, file = paste0(ifelse(dir.exists(dir), dir, dirx), "/mcmc/", save_file_mcmc)))
        }
        
        
        if (!is.null(dir)) {
            # cony <- file(paste0(ifelse(dir.exists(dir), dir, dirx), "/progress_", Sys.info()["nodename"], ".txt"), open = "a")
            cony <- file(paste0(ifelse(dir.exists(dir), dir, dirx), "/progress.txt"), open = "a")
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
    
    ### REST OLD CODE ------------------------------------------------------------------
    
    
    
    #### response styles
    #     if(length(RSmin) == 1)
    #         RSmin <- rep(RSmin, numRS.gen)
    #     if(length(RSmax) == 1)
    #         RSmax <- rep(RSmax, numRS.gen)
    #   if(length(RSmin)!=(S.gen-1) | length(RSmax)!=(S.gen-1))
    #     warning("Check definition of RSmin and RSmax")
    
    #     genNames <- c(paste0("theta[",apply(expand.grid(1:N, 1:S.gen),1, 
    #                                         paste0, collapse=","),"]"),
    #                   paste0("beta[",apply(expand.grid(1:J, 1:(numRS.gen+1)),1, 
    #                                        paste0, collapse=","),"]"), 
    #                   paste0("Sigma[",apply(expand.grid(1:S.gen, 1:S.gen),1, 
    #                                         paste0, collapse=","),"]"))
    # ttt <- c("T_obs","T_pred","post_p")
    
    # separate lists for fitting fitModel
    #     fitNames <- theta_mu <- V.ls <- list()
    #     for(sss in 1:length(S.fit)){
    #         fitNames[[sss]] <- c(paste0("theta[",apply(expand.grid(1:N, 1:S.fit[sss]),1, 
    #                                                    paste0, collapse=","),"]"),
    #                              paste0("beta[",apply(expand.grid(1:J, 1:(numRS.fit[sss]+1)),1, 
    #                                                   paste0, collapse=","),"]"), 
    #                              paste0("Sigma[",apply(expand.grid(1:S.fit[sss], 1:S.fit[sss]),1, 
    #                                                    paste0, collapse=","),"]"),
    #                              ttt)
    #         theta_mu[[sss]] <- rep(0,S.fit[sss])
    #         
    #         ### set up parameters
    #         if(missing(V)){
    #             V.ls[[sss]] <- diag(S.fit[sss])
    #         }else if(length(S.fit) == 1){
    #             V.ls[[sss]] <- unlist(V)
    #         }else{
    #             V.ls[[sss]] <- V[[sss]]
    #         }
    #     }
    #     if(missing(df)){
    #         df <- S.fit+1
    #     }
    
    
    
    
    #### function for a single CPU ####################
    
    #     boeckrep <- function(rr){
    
    #         #### generate data
    #         if(genModel == "2012"){
    #             beta.mrs <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
    #                                               a = qnorm(.5), b = qnorm(.9))
    #             beta.ers <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
    #                                               a = qnorm(.5), b = qnorm(.9))
    #             beta.trait <- truncnorm::rtruncnorm(n = J, mean = 0, sd = sqrt(.5),
    #                                                 a = qnorm(.3), b = qnorm(.7))
    #             betas_i <- cbind(beta.mrs, beta.ers, beta.trait)
    #             gen <- generate_irtree_2012(N = N, J = J, betas = betas_i, theta.vcov = theta_vcov_i,
    #                                   prop.rev = prop.rev, traitItem = traitItem, cat = TRUE)
    #         }else if (genModel == "ext"){
    #             beta.mrs <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
    #                                               a = qnorm(.5), b = qnorm(.9))
    #             beta.ers <- truncnorm::rtruncnorm(n = J, mean = qnorm(.7), sd = sqrt(.1),
    #                                               a = qnorm(.5), b = qnorm(.9))
    #             beta.ars <- truncnorm::rtruncnorm(n = J, mean = qnorm(.975), sd = sqrt(.1),
    #                                               a = qnorm(.9), b = qnorm(.999))
    #             beta.trait <- truncnorm::rtruncnorm(n = J, mean = 0, sd = sqrt(.5),
    #                                                 a = qnorm(.3), b = qnorm(.7))
    #             betas_i <- cbind(beta.mrs, beta.ers, beta.ars, beta.trait)
    # #             betas_i <- cbind(runif(J,RSmin[1],RSmax[1]),
    # #                           runif(J,RSmin[2],RSmax[2]), 
    # #                           runif(J,RSmin[3],RSmax[3]),
    # #                           rnorm(J,0,trait.sd ))
    #             gen <- generate_irtree_ext(N = N, J = J, betas = betas_i, theta.vcov = theta_vcov_i,
    #                                  prop.rev = prop.rev, traitItem = traitItem, cat = TRUE)
    #         }
    
    
    
    #         ##################### fit JAGS / Stan
    #         fitpar <- list() ; dic <- list()
    #         for(sss in 1:length(S.fit)){
    #             datalist <- list(S=S.fit[sss], df=df[sss], V=V.ls[[sss]], N=N, J=J, revItem=gen$revItem,
    #                              traitItem=gen$traitItem, X=gen$X, theta_mu=theta_mu[[sss]])
    #             if(fitMethod == "jags"){
    #                 rjags::load.module("glm", quiet=T)
    #                 #         boeck.jags <- run.jags(model=paste0( modelPath,"/mpt2irt/models/jags_boeck_",fitModel[sss],"_1d.txt"), 
    #                 #                               monitor = c("theta", "beta", "Sigma", "T_obs","T_pred", "post_p",
    #                 #                                           'deviance', 'pd', 'pd.i', 'popt', 'dic'), 
    #                 #                               data=datalist, n.chains=n.chains,  
    #                 #                               burnin = M/4, sample = M, adapt=M/4, 
    #                 #                               summarise = T, thin = thin, method="simple")
    #                 boeck.jags <- rjags::jags.model(file=paste0( modelPath,"/mpt2irt/models/jags_boeck_",
    #                                                              fitModel[sss],".txt"), 
    #                                                 data=datalist,  n.chains=n.chains, n.adapt=M/4)
    #                 adapt <- F
    #                 while(!adapt){
    #                     adapt <- rjags::adapt(boeck.jags,M/4,end.adaptation = FALSE)
    #                 }
    #                 rjags::update(boeck.jags, M/4)
    #                 boeck.samp <- coda.samples.dic(boeck.jags, 
    #                                                variable.names = c("theta", "beta", "Sigma",
    #                                                                   "T_obs","T_pred", "post_p"),
    #                                                n.iter=M*thin, thin = thin)
    #                 dic[[sss]] <- boeck.samp$dic
    #                 boeck.stan <- mcmc.list2stan(boeck.samp$samples)
    #                 rm(boeck.samp)
    #             }else{
    #                 ### Stan
    #                 data(boeck_stan_models)
    #                 if(fitModel == "ext"){
    #                     stanExe <- boeck_stan_ext
    #                 }else{
    #                     stanExe <- boeck_stan_2012
    #                 }
    #                 boeck.stan <- rstan::sampling(stanExe, 
    #                                               # model_name=paste0("Boeckenholt_", fitModel[sss]),
    #                                               pars=c("theta","beta", "Sigma","T_obs","T_pred","post_p"),
    #                                               data=datalist, chains = n.chains, thin=thin, 
    #                                               iter=M*thin)
    #                 dic[[sss]] <- NULL
    #             }
    #             
    #             fitpar[[sss]] <- rstan::monitor(boeck.stan, print=T)[fitNames[[sss]],] 
    #         }
    #         genpar <- c(gen$theta, gen$betas, gen$theta.vcov)
    #         names(genpar) <- genNames
    #         
    #         if(path != ""){
    #             try(write(paste0("\n|", paste0(rep("#",floor(rr/R*30)), collapse=""),
    #                              paste0(rep(" ",30-floor(rr/R*30)), collapse=""),
    #                              "| rr=",rr,"/",R," ; ",Sys.time()), append=T,
    #                       file=paste0(path,"/recovery_irtree_progress.txt")))
    #             if(saveTemp)
    #                 try(save(genpar, fitpar, dic, 
    #                          file=paste0(path,"/boeck_recov_",fitMethod,"_gen-",genModel,
    #                                      "_fit-",paste0(fitModel, collapse="-") ,
    #                                      "_N=",N,"_J=",J,"_rep",rr,".RData")))
    #         }
    #         
    #         # clean up to get RAM
    #         gc(T, verbose=F)
    #         res <- list(gen=genpar, fit=fitpar, dic=dic)
    #         return(res)
    #     }
    #     
    #     # paths to save temporary results
    #     modelPath <- .libPaths()[1]
    #     if(path != "")
    #         try(write(paste0("Boeckenholt Recovery Simulation: \ngenModel = ", genModel, "; fitModel = ", 
    #                          paste0(fitModel, collapse="-"), " ; Start: ",Sys.time()),
    #                   file=paste0(path,"/recovery_irtree_progress.txt")))
    #     
    #     ######## parallel computing
    #     
    #     cl <- makeCluster(nCPU)
    #     registerDoParallel(cl)
    #     
    #     runtime <- system.time({
    #         res <- foreach::foreach(rr=1:R, .packages = c("mpt2irt","rjags","rstan","coda"), 
    #                                 .combine="c") %dopar% {boeckrep(rr)}
    #     })
    #     cat("\nTotal Runtime: ",round(runtime[3]/3600, 2), "Hours\n")
    #     
    #     stopCluster(cl)
    #     
    #     ###### collect results
    #     
    #     estimate1 <-  array(NA, c(R, N*S.fit[1] + J*(numRS.fit[1]+1) + S.fit[1]*S.fit[1] + 3 ,10), 
    #                         dimnames= list(NULL, fitNames[[1]], colnames(res[[2]][[1]])))
    #     if(length(S.fit) == 2){
    #         estimate2 <-  array(NA, c(R, N*S.fit[2] + J*(numRS.fit[2]+1) + S.fit[2]*S.fit[2] + 3 ,10), 
    #                             dimnames= list(NULL, fitNames[[2]], colnames(res[[2]][[2]])))
    #     }else{
    #         estimate2 <- NULL
    #     }
    #     true <- matrix(NA, R, length(res[[1]]),  dimnames=list(NULL, names(res[[1]])))
    #     dic <- data.frame(dev1=NA, penalty1=NA, dic1=NA, dev2=NA, penalty2=NA, dic2=NA)
    #     
    #     for(r in 1:R){
    #         true[r,] <- res[[(r-1)*3+1]]
    #         estimate1[r,,] <- res[[(r-1)*3+2]][[1]]
    #         ddd <- res[[(r-1)*3+3]]
    #         try(dic[r,1:3] <- c(ddd[[1]]$dev, ddd[[1]]$pen, ddd[[1]]$dev+ ddd[[1]]$pen), silent=T)
    #         if(length(S.fit) == 2){
    #             estimate2[r,,] <- res[[(r-1)*3+2]][[2]]
    #             try(dic[r,4:6] <- c(ddd[[2]]$dev, ddd[[2]]$pen, ddd[[2]]$dev+ ddd[[2]]$pen), silent=T)
    #         }
    #     }
    #     
    #     if(path != ""){
    #         try(unlink(paste0(path,"/recovery_irtree_progress.txt")), silent=T)
    #         for(rr in 1:R){
    #             try(unlink(paste0(path,"/boeck_recov_",fitMethod,"_gen-",genModel,
    #                               "_fit-",paste0(fitModel, collapse="-") ,"_N=",N,"_J=",J,"_rep",rr,".RData")), silent=T)
    #         }
    #     }
    #     
    #     res <- list(estimate1=estimate1, estimate2=estimate2, true=true, dic=dic, N=N, J=J, S.gen=S.gen, 
    #                 S.fit=S.fit, R=R, genModel=genModel, fitModel=fitModel, fitMethod=fitMethod, 
    #                 RSmin=RSmin, RSmax=RSmax, theta.vcov=theta_vcov_i, numRS.gen=numRS.gen, numRS.fit=numRS.fit,
    #                 prop.rev=prop.rev, traitItem=traitItem, trait.sd=trait.sd, thin=thin, M=M)
    #     class(res) <- "boecksim"
    #     return(res)
}
