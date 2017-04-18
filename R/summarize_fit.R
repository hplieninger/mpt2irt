if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

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
#' res3 <- tidyup_irtree_fit(res2, N = N, J = J, revItem = dat$revItem,
#'                           traitItem = dat$traitItem, fitModel = res1$fitModel)
#' str(res3)
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
                              fitModel = c("ext", "2012", "pcm"),
                              fitMethod = c("stan", "jags"), 
                              measure = c("Median", "Mean"),
                              tt_names = NULL, 
                              # rs_names = c("m", "e", "a", "t"), trait = NULL,
                              plot = TRUE,
                              ...){
    fitModel <- match.arg(fitModel)
    fitMethod <- match.arg(fitMethod)
    measure <- match.arg(measure)
    
    checkmate::assert_int(S, lower = 1, null.ok = TRUE)
    checkmate::assert_numeric(N, lower = 1, upper = Inf, finite = TRUE,
                              any.missing = FALSE, len = 1,
                              null.ok = FALSE)
    checkmate::assert_numeric(J, lower = 1, upper = Inf, finite = TRUE,
                              any.missing = FALSE, len = 1,
                              null.ok = FALSE)
    checkmate::assert_numeric(traitItem, lower = 1, upper = Inf, finite = TRUE,
                              any.missing = FALSE, min.len = 1,
                              null.ok = FALSE)
    
    if (!("V" %in% names(fit))) {
        fit <- list("samples" = fit)
    }
    if (is.null(fitMethod)) {
        if ("samples" %in% names(fit)) {
            if (is.list(fit$samples)) {
                fitMethod <- "jags"
            } else {
                fitMethod <- "stan"
            }
        } else if ("dic" %in% names(fit)) {
            fitMethod <- "jags"
        }
    }
    
    if (!("summary" %in% names(fit))) {
        if (!("mcmc" %in% names(fit))) {
            if (fitMethod == "jags") {
                fit$mcmc <- fit$samples$mcmc
            } else {
                fit$mcmc <- rstan::As.mcmc.list(fit$samples)
            }
        }
        if(class(fit$mcmc) != "mcmc.list") stop("Unable to find or create object of class 'mcmc.list' in 'fit$mcmc'.")
        # fit$summary <- coda:::summary.mcmc.list(fit$mcmc)
        fit$summary <- summary(fit$mcmc)
    }
    
    if (!is.null(S)) {
        S2 <- S - 1 + length(unique(traitItem))
    } else if (fitModel == "ext") {
        S <- 4
        S2 <- S - 1 + length(unique(traitItem))
    } else if (fitModel == "2012") {
        S <- 3
        S2 <- S - 1 + length(unique(traitItem))
    } else if (fitModel == "ext5") {
        S <- 5
        S2 <- 3 + length(unique(traitItem))
    } else if (fitModel == "pcm") {
        S <- 4
        S2 <- length(unique(traitItem))
    }
    
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
                       "beta" = matrix(NA_real_, J, S),
                       "theta" = matrix(NA_real_, N, S2),
                       "sigma_vec" = rep(NA_real_, sum(upper.tri(diag(S2), diag = T))),
                       "Sigma" = matrix(NA_real_, S2, S2),
                       "Corr" = matrix(NA_real_, S2, S2),
                       "mu_beta" = rep(NA_real_, ifelse(fitModel == "ext5", S2+1, S2)),
                       "sigma_beta" = rep(NA_real_, ifelse(fitModel == "ext5", S2+1, S2)),
                       "beta_ARS_extreme" = NA_real_)
        return_list[[iii]] <- setNames(list(tmp1, tmp1, tmp1, tmp1),
                                       nm = c("Mean", "Median", "Q_025", "Q_975"))
    }
    
    ### Retrive estimates and save them in 'return_list'
    
    for(i in 1:S2){
        return_list$theta$Mean[, i] <- fit$summary$statistics[paste0("theta[",1:N,",",i,"]"), "Mean"]
        return_list$theta$Median[, i] <- fit$summary$quantiles[paste0("theta[",1:N,",",i,"]"), "50%"]
        return_list$theta$Q_025[, i] <- fit$summary$quantiles[paste0("theta[",1:N,",",i,"]"), "2.5%"]
        return_list$theta$Q_975[, i] <- fit$summary$quantiles[paste0("theta[",1:N,",",i,"]"), "97.5%"]
    }
    for(i in 1:S){
        return_list$beta$Mean[, i] <- fit$summary$statistics[paste0("beta[",1:J,",",i,"]"), "Mean"]
        return_list$beta$Median[, i] <- fit$summary$quantiles[paste0("beta[",1:J,",",i,"]"), "50%"]
        return_list$beta$Q_025[, i] <- fit$summary$quantiles[paste0("beta[",1:J,",",i,"]"), "2.5%"]
        return_list$beta$Q_975[, i] <- fit$summary$quantiles[paste0("beta[",1:J,",",i,"]"), "97.5%"]
    }
    sigma_idx <- which(upper.tri(diag(S2), diag = T) == TRUE, arr.ind = T)
    return_list$sigma_vec$Mean <- fit$summary$statistics[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "Mean"]
    return_list$sigma_vec$Median <- fit$summary$quantiles[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "50%"]
    return_list$sigma_vec$Q_025 <- fit$summary$quantiles[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "2.5%"]
    return_list$sigma_vec$Q_975 <- fit$summary$quantiles[paste0("Sigma[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "97.5%"]
    # rm(sigma_idx)
    
    return_list$Sigma$Mean[upper.tri(diag(S2), T)] <- return_list$sigma_vec$Mean
    return_list$Sigma$Mean[lower.tri(diag(S2))] <- t(return_list$Sigma$Mean)[lower.tri(diag(S2))]
    return_list$Sigma$Median[upper.tri(diag(S2), T)] <- return_list$sigma_vec$Median
    return_list$Sigma$Median[lower.tri(diag(S2))] <- t(return_list$Sigma$Median)[lower.tri(diag(S2))]
    return_list$Sigma$Q_025[upper.tri(diag(S2), T)] <- return_list$sigma_vec$Q_025
    return_list$Sigma$Q_025[lower.tri(diag(S2))] <- t(return_list$Sigma$Q_025)[lower.tri(diag(S2))]
    return_list$Sigma$Q_975[upper.tri(diag(S2), T)] <- return_list$sigma_vec$Q_975
    return_list$Sigma$Q_975[lower.tri(diag(S2))] <- t(return_list$Sigma$Q_975)[lower.tri(diag(S2))]
    return_list$sigma_vec <- NULL
    
    # Correlations
    return_list$corr_vec$Mean <- fit$summary$statistics[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "Mean"]
    return_list$corr_vec$Median <- fit$summary$quantiles[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "50%"]
    return_list$corr_vec$Q_975 <- fit$summary$quantiles[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "97.5%"]
    return_list$corr_vec$Q_025 <- fit$summary$quantiles[paste0("Corr[", sigma_idx[, 1], ",", sigma_idx[, 2] ,"]"), "2.5%"]
    rm(sigma_idx)
    
    return_list$Corr$Mean[upper.tri(diag(S2), T)] <- return_list$corr_vec$Mean
    return_list$Corr$Mean[lower.tri(diag(S2))] <- t(return_list$Corr$Mean)[lower.tri(diag(S2))]
    return_list$Corr$Median[upper.tri(diag(S2), T)] <- return_list$corr_vec$Median
    return_list$Corr$Median[lower.tri(diag(S2))] <- t(return_list$Corr$Median)[lower.tri(diag(S2))]
    return_list$Corr$Q_025[upper.tri(diag(S2), T)] <- return_list$corr_vec$Q_025
    return_list$Corr$Q_025[lower.tri(diag(S2))] <- t(return_list$Corr$Q_025)[lower.tri(diag(S2))]
    return_list$Corr$Q_975[upper.tri(diag(S2), T)] <- return_list$corr_vec$Q_975
    return_list$Corr$Q_975[lower.tri(diag(S2))] <- t(return_list$Corr$Q_975)[lower.tri(diag(S2))]
    return_list$corr_vec <- NULL
    
    if (fitModel == "ext5") {
        return_list$sigma_beta$Mean <- fit$summary$statistics[paste0("sigma_beta[", 1:(S2+1) ,"]"), "Mean"]
        return_list$sigma_beta$Median <- fit$summary$quantiles[paste0("sigma_beta[", 1:(S2+1) ,"]"), "50%"]
        return_list$sigma_beta$Q_025 <- fit$summary$quantiles[paste0("sigma_beta[", 1:(S2+1) ,"]"), "2.5%"]
        return_list$sigma_beta$Q_975 <- fit$summary$quantiles[paste0("sigma_beta[", 1:(S2+1) ,"]"), "97.5%"]
        return_list$mu_beta$Mean <- fit$summary$statistics[paste0("mu_beta[", 1:(S2+1) ,"]"), "Mean"]
        return_list$mu_beta$Median <- fit$summary$quantiles[paste0("mu_beta[", 1:(S2+1) ,"]"), "50%"]
        return_list$mu_beta$Q_025 <- fit$summary$quantiles[paste0("mu_beta[", 1:(S2+1) ,"]"), "2.5%"]
        return_list$mu_beta$Q_975 <- fit$summary$quantiles[paste0("mu_beta[", 1:(S2+1) ,"]"), "97.5%"]
    } else {
        return_list$sigma_beta$Mean <- fit$summary$statistics[paste0("sigma_beta[", 1:S2 ,"]"), "Mean"]
        return_list$sigma_beta$Median <- fit$summary$quantiles[paste0("sigma_beta[", 1:S2 ,"]"), "50%"]
        return_list$sigma_beta$Q_025 <- fit$summary$quantiles[paste0("sigma_beta[", 1:S2 ,"]"), "2.5%"]
        return_list$sigma_beta$Q_975 <- fit$summary$quantiles[paste0("sigma_beta[", 1:S2 ,"]"), "97.5%"]
        return_list$mu_beta$Mean <- fit$summary$statistics[paste0("mu_beta[", 1:S2 ,"]"), "Mean"]
        return_list$mu_beta$Median <- fit$summary$quantiles[paste0("mu_beta[", 1:S2 ,"]"), "50%"]
        return_list$mu_beta$Q_025 <- fit$summary$quantiles[paste0("mu_beta[", 1:S2 ,"]"), "2.5%"]
        return_list$mu_beta$Q_975 <- fit$summary$quantiles[paste0("mu_beta[", 1:S2 ,"]"), "97.5%"]
    }
    
    if (grepl("ext", fitModel) & fitModel != "ext5") {
        return_list$beta_ARS_extreme$Mean <- fit$summary$statistics["beta_ARS_extreme", "Mean"]
        return_list$beta_ARS_extreme$Median <- fit$summary$quantiles["beta_ARS_extreme", "50%"]
        return_list$beta_ARS_extreme$Q_025 <- fit$summary$quantiles["beta_ARS_extreme", "2.5%"]
        return_list$beta_ARS_extreme$Q_975 <- fit$summary$quantiles["beta_ARS_extreme", "97.5%"]
    } else {
        return_list$beta_ARS_extreme <- NULL
    }

    ### Save further settings in 'return_list'
    
    return_list$fitMethod = fitMethod
    return_list$fitModel = fitModel
    return_list$revItem = revItem
    return_list$traitItem = traitItem

    try(return_list$dic <- fit$samples$dic, silent = TRUE)
    return_list$setup <- list("V" = fit$V,
                              "df" = fit$df,
                              "startSmall" = fit$startSmall)
    if (fitMethod == "stan") {
        return_list$setup$warmup <- fit$samples@stan_args[[1]]$warmup
        return_list$setup$M <- fit$samples@stan_args[[1]]$iter - fit$samples@stan_args[[1]]$warmup
        return_list$setup$thin <- fit$samples@stan_args[[1]]$thin
        return_list$setup$chains <- length(fit$samples@stan_args)
        return_list$date <- fit$samples@date
    } else {
        return_list$setup$warmup <- unique(c(fit$samples$burnin,
                                             fit$summary$start - 1))
        return_list$setup$M <- unique(c(fit$samples$sample,
                                        length(seq(fit$summary$start, fit$summary$end, fit$summary$thin))))
        return_list$setup$thin <- unique(c(fit$samples$thin,
                                           fit$summary$thin))
        return_list$setup$chains <- fit$summary$nchain
        if ("date" %in% names(fit$samples)) {
            return_list$date <- fit$samples$date
        } else {
            return_list$date <- fit$date
        }
    }
    
    if (plot == TRUE) {
        return_list$plot <- invisible(plot_irtree(fit, S = S, J = J, revItem = revItem,
                                                  traitItem = traitItem,
                                                  # trait = trait, rs_names = rs_names,
                                                  return_data = FALSE, tt_names = tt_names,
                                                  measure = measure, fitMethod = fitMethod,
                                                  ...))
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
#' \code{\link[coda]{summary.mcmc.list}} is that, herein, the correlations of the thetas
#' are summarized in addition to the covariances.
#' 
#' @param fit a fitted object from \code{\link{fit_irtree}}.
#' @param interact logical. If set to \code{TRUE}, the function may be prompt for input.
#' @param ... Further arguments passed to \code{\link[coda]{summary.mcmc.list}}.
#' @return Returns a list containing the input \code{fit} as well as an MCMC list and \code{summary}.
#' @inheritParams fit_irtree
#' @importFrom magrittr %>% 
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
    
    if (!("V" %in% names(fit))) {
        fit <- list("samples" = fit)
    }
    
    if (is.null(fitMethod)) {
        if ("samples" %in% names(fit)) {
            if (is.list(fit$samples)) {
                fitMethod <- "jags"
            } else {
                fitMethod <- "stan"
            }
        } else if ("dic" %in% names(fit)) {
            fitMethod <- "jags"
        }
    }
    if (!("mcmc" %in% names(fit))) {
        if (fitMethod == "stan") {
            fit$mcmc <- rstan::As.mcmc.list(fit$samples)
        } else if (fitMethod == "jags") {
            fit$mcmc <- fit$samples$mcmc
        }
    }
    
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
    fit$summary <- summary(fit$mcmc, ...)
    
    return(fit)
}
