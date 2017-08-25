if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Summarize and tidy up a fitted model.
#' 
#' Function takes a fitted model returned from \code{\link{fit_irtree}} and
#' returns a tidy summary of the fitted parameters.
#' 
#' @param N number of persons
#' @param fit a fitted object from \code{\link{fit_irtree}} or preferably from \code{\link{summarize_irtree_fit}}.
#' @param plot Logical. Whether a plot should be produced using \code{\link{plot_irtree}} and returned.
#' @param ... Further arguments passed to \code{\link{plot_irtree}}.
#' @inheritParams fit_irtree
#' @inheritParams plot_irtree
#' @return The function returns a list for each of the core parameters (e.g., theta, beta), which are each summarized with the posterior mean, median and 95% CI.
#' @import coda
#' @examples
#' \dontrun{
#' # generate data
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' 
#' # fit model
#' res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200)
#' res2 <- summarize_irtree_fit(res1)
#' res3 <- tidyup_irtree_fit(res2)
#' names(res3)
#' res3$plot
#' cor(res3$theta$Median, dat$theta)
#' cor(res3$beta$Median, dat$betas)
#' }
#' @export
tidyup_irtree_fit <- function(fit,
                              S = NULL,
                              N = NULL,
                              J = NULL,
                              revItem = NULL,
                              traitItem = NULL,
                              fitModel = NULL,
                              fitMethod = NULL, 
                              plot = TRUE,
                              ...){
    
    checkmate::qassert(fit, "L+")
    checkmate::assert_int(S, lower = 1,
                          null.ok = !is.null(fit$args$S))
    checkmate::assert_int(N, lower = 1,
                          null.ok = !is.null(fit$args$N))
    checkmate::assert_int(J, lower = 1,
                          null.ok = !is.null(fit$args$J))
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 min.len = 1,
                                 null.ok = !is.null(fit$args$revItem))
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 min.len = 1,
                                 null.ok = !is.null(fit$args$traitItem))
    checkmate::assert_character(fitModel, min.chars = 1, any.missing = FALSE,
                                len = 1,
                                null.ok = !is.null(fit$args$fitModel))
    checkmate::assert_character(fitMethod, min.chars = 1, any.missing = FALSE,
                                len = 1, 
                                null.ok = !is.null(fit$args$fitMethod))
    checkmate::qassert(plot, "B1")
    
    if (!any(names(fit) %in% c("args", "V"))) {
        fit <- list("samples" = fit)
    }
    
    if (is.null(S)) S <- fit$args$S
    if (is.null(J)) J <- fit$args$J
    if (is.null(N)) N <- fit$args$N
    if (is.null(revItem)) revItem <- fit$args$revItem
    if (is.null(traitItem)) traitItem <- fit$args$traitItem
    if (is.null(fitMethod)) fitMethod <- fit$args$fitMethod
    if (is.null(fitModel)) fitModel <- fit$args$fitModel
    
    
    if (!("summary" %in% names(fit))) {
        if (!("mcmc" %in% names(fit))) {
            if (fitMethod == "jags") {
                fit$mcmc <- fit$samples$mcmc
            } else {
                fit$mcmc <- rstan::As.mcmc.list(fit$samples)
            }
        }
        if (class(fit$mcmc) != "mcmc.list") stop("Unable to find or create object of class 'mcmc.list' in 'fit$mcmc'.")
        # fit$summary <- coda:::summary.mcmc.list(fit$mcmc)
        fit$summary <- summary(fit$mcmc)
        flag1 <- FALSE
    } else {
        flag1 <- TRUE
    }
    
    # S_t:  Number of theta dimensions
    # S_b1: Number of columns in beta matrix
    # S_b2: Number of beta dimensions
    
    S_t  <- switch(fitModel,
                   S)
    S_b1 <- switch(fitModel,
                   "2012"  = 3,
                   "ext"   = 4,
                   "pcm"   = 4,
                   "steps" = 4,
                   "shift" = 3,
                   S)
    # S_b2 <- switch(fitModel,
    #                # "2012"  = 3,
    #                # "ext"   = 4,
    #                "pcm"   = 4*length(unique(traitItem)),
    #                "steps" = 4*length(unique(traitItem)),
    #                "shift" = S - 1,
    #                S)
    S_b2 <- switch(fitModel,
                   # "2012"  = 3,
                   # "ext"   = 4,
                   "pcm"   = c(length(unique(traitItem)), 4),
                   "steps" = c(length(unique(traitItem)), 4),
                   "shift" = c(1, S - 1),
                   c(1, S))

    ### Initialize list of empty lists to store estimates in
    
    return_list <- list("beta" = list(),
                        "theta" = list(),
                        "sigma_vec" = list(),
                        "Sigma" = list(),
                        "Corr" = list(),
                        "mu_beta" = list(),
                        "sigma_beta" = list(),
                        "beta_ARS_extreme" = list())
    
    for (iii in seq_along(return_list)) {
        tmp1 <- switch(names(return_list)[iii],
                       "beta" = matrix(NA_real_, J, S_b1),
                       "theta" = matrix(NA_real_, N, S_t),
                       "sigma_vec" = rep(NA_real_, sum(upper.tri(diag(S_t), diag = T))),
                       "Sigma" = matrix(NA_real_, S_t, S_t),
                       "Corr" = matrix(NA_real_, S_t, S_t),
                       # "mu_beta" = rep(NA_real_, ifelse(fitModel == "ext5", S_t + 1, S_t)),
                       "mu_beta" = matrix(NA_real_, S_b2[1], S_b2[2]),
                       # "sigma_beta" = rep(NA_real_, ifelse(fitModel == "ext5", S_t + 1, S_t)),
                       "sigma_beta" = matrix(NA_real_, S_b2[1], S_b2[2]),
                       "beta_ARS_extreme" = NA_real_)
        return_list[[iii]] <- setNames(list(tmp1, tmp1, tmp1, tmp1),
                                       nm = c("Mean", "Median", "Q_025", "Q_975"))
    }
    return_list$args <- fit$args
    
    ### Retrive estimates and save them in 'return_list'
    
    for (i in 1:S_t){
        return_list$theta$Mean[, i] <- fit$summary$statistics[paste0("theta[",1:N,",",i,"]"), "Mean"]
        return_list$theta$Median[, i] <- fit$summary$quantiles[paste0("theta[",1:N,",",i,"]"), "50%"]
        return_list$theta$Q_025[, i] <- fit$summary$quantiles[paste0("theta[",1:N,",",i,"]"), "2.5%"]
        return_list$theta$Q_975[, i] <- fit$summary$quantiles[paste0("theta[",1:N,",",i,"]"), "97.5%"]
    }
    for (i in 1:S_b1){
        return_list$beta$Mean[, i] <- fit$summary$statistics[paste0("beta[",1:J,",",i,"]"), "Mean"]
        return_list$beta$Median[, i] <- fit$summary$quantiles[paste0("beta[",1:J,",",i,"]"), "50%"]
        return_list$beta$Q_025[, i] <- fit$summary$quantiles[paste0("beta[",1:J,",",i,"]"), "2.5%"]
        return_list$beta$Q_975[, i] <- fit$summary$quantiles[paste0("beta[",1:J,",",i,"]"), "97.5%"]
    }
    sigma_idx <- which(upper.tri(diag(S_t), diag = T) == TRUE, arr.ind = T)
    return_list$sigma_vec$Mean <- fit$summary$statistics[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "Mean"]
    return_list$sigma_vec$Median <- fit$summary$quantiles[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "50%"]
    return_list$sigma_vec$Q_025 <- fit$summary$quantiles[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "2.5%"]
    return_list$sigma_vec$Q_975 <- fit$summary$quantiles[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "97.5%"]
    # rm(sigma_idx)
    
    return_list$Sigma$Mean[upper.tri(diag(S_t), T)] <- return_list$sigma_vec$Mean
    return_list$Sigma$Mean[lower.tri(diag(S_t))] <- t(return_list$Sigma$Mean)[lower.tri(diag(S_t))]
    return_list$Sigma$Median[upper.tri(diag(S_t), T)] <- return_list$sigma_vec$Median
    return_list$Sigma$Median[lower.tri(diag(S_t))] <- t(return_list$Sigma$Median)[lower.tri(diag(S_t))]
    return_list$Sigma$Q_025[upper.tri(diag(S_t), T)] <- return_list$sigma_vec$Q_025
    return_list$Sigma$Q_025[lower.tri(diag(S_t))] <- t(return_list$Sigma$Q_025)[lower.tri(diag(S_t))]
    return_list$Sigma$Q_975[upper.tri(diag(S_t), T)] <- return_list$sigma_vec$Q_975
    return_list$Sigma$Q_975[lower.tri(diag(S_t))] <- t(return_list$Sigma$Q_975)[lower.tri(diag(S_t))]
    return_list$sigma_vec <- NULL
    
    # Correlations
    if (flag1 == TRUE) {
        return_list$corr_vec$Mean <- fit$summary$statistics[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "Mean"]
        return_list$corr_vec$Median <- fit$summary$quantiles[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "50%"]
        return_list$corr_vec$Q_975 <- fit$summary$quantiles[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "97.5%"]
        return_list$corr_vec$Q_025 <- fit$summary$quantiles[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "2.5%"]
        rm(sigma_idx)
        
        return_list$Corr$Mean[upper.tri(diag(S_t), T)] <- return_list$corr_vec$Mean
        return_list$Corr$Mean[lower.tri(diag(S_t))] <- t(return_list$Corr$Mean)[lower.tri(diag(S_t))]
        return_list$Corr$Median[upper.tri(diag(S_t), T)] <- return_list$corr_vec$Median
        return_list$Corr$Median[lower.tri(diag(S_t))] <- t(return_list$Corr$Median)[lower.tri(diag(S_t))]
        return_list$Corr$Q_025[upper.tri(diag(S_t), T)] <- return_list$corr_vec$Q_025
        return_list$Corr$Q_025[lower.tri(diag(S_t))] <- t(return_list$Corr$Q_025)[lower.tri(diag(S_t))]
        return_list$Corr$Q_975[upper.tri(diag(S_t), T)] <- return_list$corr_vec$Q_975
        return_list$Corr$Q_975[lower.tri(diag(S_t))] <- t(return_list$Corr$Q_975)[lower.tri(diag(S_t))]
        return_list$corr_vec <- NULL
    } else {
        message("Output will not contain latent correlations but only covariances.\n",
                "You may run summarize_irtree_fit() first.")
    }
    
    if (S_b2[1] == 1) {
        tmp1 <- grep("sigma_beta\\[.*[0-9]+\\]", rownames(fit$summary$statistics))
        tmp2 <- grep("mu_beta\\[.*[0-9]+\\]", rownames(fit$summary$statistics))
        # return_list$sigma_beta$Mean <- fit$summary$statistics[paste0("sigma_beta[", 1:S_b2[2] ,"]"), "Mean"]
        # return_list$sigma_beta$Median <- fit$summary$quantiles[paste0("sigma_beta[", 1:S_b2[2] ,"]"), "50%"]
        # return_list$sigma_beta$Q_025 <- fit$summary$quantiles[paste0("sigma_beta[", 1:S_b2[2] ,"]"), "2.5%"]
        # return_list$sigma_beta$Q_975 <- fit$summary$quantiles[paste0("sigma_beta[", 1:S_b2[2] ,"]"), "97.5%"]
        return_list$sigma_beta$Mean   <- fit$summary$statistics[tmp1, "Mean"]
        return_list$sigma_beta$Median <- fit$summary$quantiles[ tmp1, "50%"]
        return_list$sigma_beta$Q_025  <- fit$summary$quantiles[ tmp1, "2.5%"]
        return_list$sigma_beta$Q_975  <- fit$summary$quantiles[ tmp1, "97.5%"]
        # return_list$mu_beta$Mean[1,]   <- fit$summary$statistics[paste0("mu_beta[", 1:S_b2[2] ,"]"), "Mean"]
        # return_list$mu_beta$Median[1,] <- fit$summary$quantiles[paste0("mu_beta[", 1:S_b2[2] ,"]"), "50%"]
        # return_list$mu_beta$Q_025[1,]  <- fit$summary$quantiles[paste0("mu_beta[", 1:S_b2[2] ,"]"), "2.5%"]
        # return_list$mu_beta$Q_975[1,]  <- fit$summary$quantiles[paste0("mu_beta[", 1:S_b2[2] ,"]"), "97.5%"]
        return_list$mu_beta$Mean[1,]   <- fit$summary$statistics[tmp2, "Mean"]
        return_list$mu_beta$Median[1,] <- fit$summary$quantiles[ tmp2, "50%"]
        return_list$mu_beta$Q_025[1,]  <- fit$summary$quantiles[ tmp2, "2.5%"]
        return_list$mu_beta$Q_975[1,]  <- fit$summary$quantiles[ tmp2, "97.5%"]
    } else {
        for (i in seq_len(S_b2[1])) {
            return_list$sigma_beta$Mean[i,]   <- fit$summary$statistics[paste0("sigma_beta[", i, ",", 1:S_b2[2], "]"), "Mean"]
            return_list$sigma_beta$Median[i,] <- fit$summary$quantiles[ paste0("sigma_beta[", i, ",", 1:S_b2[2], "]"), "50%"]
            return_list$sigma_beta$Q_025[i,]  <- fit$summary$quantiles[ paste0("sigma_beta[", i, ",", 1:S_b2[2], "]"), "2.5%"]
            return_list$sigma_beta$Q_975[i,]  <- fit$summary$quantiles[ paste0("sigma_beta[", i, ",", 1:S_b2[2], "]"), "97.5%"]
            return_list$mu_beta$Mean[i,]   <- fit$summary$statistics[paste0("mu_beta[", i, ",", 1:S_b2[2], "]"), "Mean"]
            return_list$mu_beta$Median[i,] <- fit$summary$quantiles[ paste0("mu_beta[", i, ",", 1:S_b2[2], "]"), "50%"]
            return_list$mu_beta$Q_025[i,]  <- fit$summary$quantiles[ paste0("mu_beta[", i, ",", 1:S_b2[2], "]"), "2.5%"]
            return_list$mu_beta$Q_975[i,]  <- fit$summary$quantiles[ paste0("mu_beta[", i, ",", 1:S_b2[2], "]"), "97.5%"]
        }
        
    }
    # if (fitModel == "ext5") {
    #     return_list$sigma_beta$Mean <- fit$summary$statistics[paste0("sigma_beta[", 1:(S_b2+1) ,"]"), "Mean"]
    #     return_list$sigma_beta$Median <- fit$summary$quantiles[paste0("sigma_beta[", 1:(S_b2+1) ,"]"), "50%"]
    #     return_list$sigma_beta$Q_025 <- fit$summary$quantiles[paste0("sigma_beta[", 1:(S_b2+1) ,"]"), "2.5%"]
    #     return_list$sigma_beta$Q_975 <- fit$summary$quantiles[paste0("sigma_beta[", 1:(S_b2+1) ,"]"), "97.5%"]
    #     return_list$mu_beta$Mean <- fit$summary$statistics[paste0("mu_beta[", 1:(S_b2+1) ,"]"), "Mean"]
    #     return_list$mu_beta$Median <- fit$summary$quantiles[paste0("mu_beta[", 1:(S_b2+1) ,"]"), "50%"]
    #     return_list$mu_beta$Q_025 <- fit$summary$quantiles[paste0("mu_beta[", 1:(S_b2+1) ,"]"), "2.5%"]
    #     return_list$mu_beta$Q_975 <- fit$summary$quantiles[paste0("mu_beta[", 1:(S_b2+1) ,"]"), "97.5%"]
    # }
    
    if (grepl("ext", fitModel) & fitModel != "ext5") {
        return_list$beta_ARS_extreme$Mean <- fit$summary$statistics["beta_ARS_extreme", "Mean"]
        return_list$beta_ARS_extreme$Median <- fit$summary$quantiles["beta_ARS_extreme", "50%"]
        return_list$beta_ARS_extreme$Q_025 <- fit$summary$quantiles["beta_ARS_extreme", "2.5%"]
        return_list$beta_ARS_extreme$Q_975 <- fit$summary$quantiles["beta_ARS_extreme", "97.5%"]
    } else {
        return_list$beta_ARS_extreme <- NULL
    }
    
    ### Save further settings in 'return_list'

    try(return_list$dic <- fit$samples$dic, silent = TRUE)
    return_list$date <- switch(fitMethod,
                               "stan" = fit$samples@date,
                               "jags" = fit$samples$date)
    # if (fitMethod == "stan") {
    #     return_list$date <- fit$samples@date
    # } else {
    #     if ("date" %in% names(fit$samples)) {
    #         return_list$date <- fit$samples$date
    #     } else {
    #         return_list$date <- fit$date
    #     }
    # }
    
    if (plot == TRUE) {
        return_list$plot <- invisible(plot_irtree(fit, S = S_b1, J = J, revItem = revItem,
                                                  traitItem = traitItem,
                                                  return_data = FALSE, 
                                                  fitMethod = fitMethod, ...))
    }
    
    return(return_list)
}

#' Calculate summary for a fitted model.
#' 
#' Function takes a fitted model returned from \code{\link{fit_irtree}} and
#' calculates summaries for all parameters using coda's
#' \code{\link[coda]{summary.mcmc.list}}.
#' 
#' The difference between the present function and directly calling coda's 
#' \code{\link[coda]{summary.mcmc.list}} is that, herein, the correlations of
#' the thetas are summarized in addition to the covariances. This is useful if
#' \code{fitMethod = "jags"}, because JAGS does only save the covariances; in
#' Stan, the correlations are computed directly during sampling.
#' 
#' @param fit a fitted object from \code{\link{fit_irtree}}.
#' @param interact logical. If set to \code{TRUE}, the function may be prompt for input.
#' @param ... Further arguments passed to \code{\link[coda]{summary.mcmc.list}}.
#' @return Returns a list containing the input \code{fit} as well as an MCMC list and \code{summary}.
#' @inheritParams fit_irtree
# @importFrom magrittr %>% 
#' @import coda
#' @examples
#' \dontrun{
#' # generate data
#' N <- 20
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)
#' 
#' # fit model
#' res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200)
#' res2 <- summarize_irtree_fit(res1)
#' head(res2$summary$statistics)
#' head(res2$summary$quantiles)
#' }
#' @export
summarize_irtree_fit <- function(fit,
                                 fitMethod = NULL,
                                 interact = FALSE,
                                 ...) {
    checkmate::qassert(fit, "L+")
    checkmate::assert_character(fitMethod, min.chars = 1, any.missing = FALSE,
                                len = 1, 
                                null.ok = !is.null(fit$args$fitMethod))
    checkmate::qassert(interact, "B1")
    
    if (interact == TRUE & object.size(fit) > 1000000000) {
        proceed <- switch(menu(c("Yes, proceed", "No, thank you."),
                               title = paste0("The object 'fit' is larger than ",
                                              format(object.size(fit), "Gb"),
                                              ", you may run out of memory. Do you want to proceed?")) + 1,
                          FALSE, TRUE, FALSE)
        if (proceed == FALSE) {
            message("I'm terminating the function call and invisibly returning 'fit'.")
            return(invisible(fit))
        }
    }
    
    if (!any((names(fit) %in% c("V", "args")))) {
        fit <- list("samples" = fit)
    }
    
    if (is.null(fitMethod)) {
        fitMethod <- fit$args$fitMethod
        # if ("samples" %in% names(fit)) {
        #     if (is.list(fit$samples)) {
        #         fitMethod <- "jags"
        #     } else {
        #         fitMethod <- "stan"
        #     }
        # } else if ("dic" %in% names(fit)) {
        #     fitMethod <- "jags"
        # }
    }
    if (!("mcmc" %in% names(fit))) {
        if (fitMethod == "stan") {
            fit$mcmc <- rstan::As.mcmc.list(fit$samples)
        } else if (fitMethod == "jags") {
            fit$mcmc <- fit$samples$mcmc
            fit$samples$mcmc <- window(fit$samples$mcmc,
                                       start = attr(fit$samples$mcmc[[1]], "mcpar")[2] - 1,
                                       end = attr(fit$samples$mcmc[[1]], "mcpar")[2],
                                       thin = 1)
        }
    }
    
    if (!("Corr[1,1]" %in% colnames(fit$mcmc[[1]]))) {
        # for every chain do (via purrr::map) cov2cor and save correlations in tmp1
        tmp1 <- fit$mcmc %>%
            purrr::map(~ runjags::combine.mcmc(., vars = "^Sigma", collapse.chains = F) %>% 
                           # magrittr::extract(1:2,)) %>%
                           split(., 1:nrow(.)) %>% 
                           purrr::map(~ matrix(., sqrt(length(.)), sqrt(length(.)))) %>% 
                           purrr::map(cov2cor) %>% 
                           purrr::map(as.vector) %>%
                           data.frame %>% 
                           t %>% 
                           magrittr::set_rownames(NULL) %>% {
                               tmp1 <- sqrt(ncol(.))
                               tmp2 <- paste0("Corr[", 
                                              rep(seq(1, tmp1), each = tmp1),
                                              ",",
                                              rep(seq(1, tmp1), times = tmp1),
                                              "]")
                               magrittr::set_colnames(., tmp2)
                           })
        
        # cbind every chain with correponding correlations
        fit$mcmc <- tmp1  %>% 
            purrr::map2(fit$mcmc, ., cbind) %>% 
            purrr::map(~ do.call(coda::mcmc, args = c(list(data = .), as.list(attr(fit$mcmc[[1]], "mcpar"))))) %>% 
            # purrr::map(~ coda::mcmc(., 301, 600, 1)) %>% 
            coda::mcmc.list()
    }
    
    # for (iii in 1:length(fit$mcmc)) {
    #     tmp1 <- fit$mcmc %>% 
    #         extract2(iii) %>% 
    #         runjags::combine.mcmc(., vars = "^Sigma", collapse.chains = F)
    #     tmp2 <- tmp1 %>% 
    #         split(., 1:nrow(.)) %>% 
    #         purrr::map(~ matrix(., sqrt(length(.)), sqrt(length(.)))) %>% 
    #         purrr::map(cov2cor) %>% 
    #         purrr::map(as.vector) %>% 
    #         data.frame %>% 
    #         t %>%
    #         set_colnames(gsub("Sigma", "Corr", colnames(tmp1))) %>% 
    #         set_rownames(NULL)
    #     # fit_mcmc[[iii]] <- cbind(fit_mcmc[[iii]], tmp2)
    #     tmp3 <- tmp2 %>%
    #         cbind(fit$mcmc[[iii]], .) %>% 
    #         coda::mcmc() %T>% 
    #         {attr(., "mcpar") <- attr(fit$mcmc[[iii]], "mcpar")}
    #     fit$mcmc[[iii]] <- tmp3
    # }
    
    # fit$summary <- coda:::summary.mcmc.list(fit$mcmc, ...)
    fit$summary <- list()
    try({
        if (fitMethod == "stan") {
            sum1 <- shinystan::as.shinystan(fit$samples)
        } else {
            sum1 <- shinystan::as.shinystan(fit$mcmc)
        }
        fit$summary$statistics <- sum1@summary[, c("mean", "se_mean", "sd", "n_eff", "Rhat")]
        colnames(fit$summary$statistics)[1] <- "Mean"
        fit$summary$quantiles  <- sum1@summary[, !(colnames(sum1@summary) %in% c("mean", "se_mean", "sd", "n_eff", "Rhat"))]
        fit$summary$start  <- attr(fit$mcmc[[1]], "mcpar")[1]
        fit$summary$end    <- attr(fit$mcmc[[1]], "mcpar")[2]
        fit$summary$thin   <- attr(fit$mcmc[[1]], "mcpar")[3]
        fit$summary$nchain <- length(fit$mcmc)
    }, silent = TRUE)
    if (length(fit$summary) == 0) {
        fit$summary <- summary(fit$mcmc, ...)
    }
    
    return(fit)
}
