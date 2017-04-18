if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Draw posterior predictive values.
#' 
#' Function takes an MCMC list and draws posterior predictive values.
#' 
#' Given an MCMC list, this function draws theta values for \code{N} persons
#' from the "current" multivariate normal distribution and---given the current
#' item parameters---calculates the probabilities (given the fitted model) for
#' every item--person combination. Posterior predictive responses (i..e, an
#' integer between 1 and 5) are then drawn based on these probabilities and this
#' is done \code{iter} times per chain. In each iteration and for every item,
#' the responses are aggregated across persons resulting in a frequency table
#' (i.e., a distribution), and---across iterations---this finally results in a
#' posterior predictive distribution for every category for every item. This
#' distribution is summarized with the quantiles given in \code{probs} and a
#' data frame containing these quantiles for every item is returned.
#' 
#' @param iter Numeric. Number of iterations per chain to use for generating
#'   posterior predictives.
#' @param probs Numeric. Vector of probabilities (passed to
#'   \code{\link[stats]{quantile}}) that is used to calculate quantiles of the
#'   posterior predictive distribution. Argument may be named (used to generate
#'   colnames). Defaults to \code{c(.025, .975, .16, .84)}.
#' @param N Numeric. Number of persons for whom posterior predictives should be
#'   drawn.
#' @return Returns a data frame containing the specified quantiles of the
#'   posterior predictive distribution for every category of every item.
#' @inheritParams fit_irtree
#' @inheritParams runjags::combine.mcmc
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import coda
#' @examples 
#' \dontrun{
#' J <- 10
#' betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
#' dat <- generate_irtree_ext(N = 20, J = J, betas = betas, beta_ARS_extreme = .5)
#' 
#' # fit model
#' res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200)
#' res2 <- summarize_irtree_fit(res1)
#' 
#' # posterior predictives for 512 hypothetical persons
#' res3 <- pp_irtree(res2$mcmc, iter = 10, N = 512, traitItem = dat$traitItem,
#'                   revItem = dat$revItem, fitModel = res2$fitModel)
#' }
#' @export
pp_irtree <- function(mcmc.objects,
                     iter = 100,
                     probs = NULL,
                     N = NULL,
                     traitItem = NULL,
                     revItem = NULL,
                     fitModel = NULL) {
    
    checkmate::assert_numeric(iter, lower = 1, upper = Inf, finite = TRUE,
                              any.missing = FALSE, len = 1,
                              null.ok = FALSE)
    checkmate::assert_numeric(probs, lower = 0, upper = 1, finite = TRUE,
                              any.missing = FALSE,  min.len = 1, unique = TRUE,
                              null.ok = TRUE)
    checkmate::assert_numeric(N, lower = 1, upper = Inf, finite = TRUE,
                              any.missing = FALSE, len = 1,
                              null.ok = FALSE)
    checkmate::assert_numeric(traitItem, lower = 1, upper = Inf, finite = TRUE,
                              any.missing = FALSE, min.len = 1,
                              null.ok = FALSE)
    checkmate::assert_numeric(revItem, lower = 0, upper = 1, finite = TRUE,
                              any.missing = FALSE, len = length(traitItem),
                              null.ok = FALSE)
    checkmate::assert_string(fitModel, na.ok = FALSE, min.chars = 1, 
                             null.ok = FALSE)
    
    if (is.null(probs)) {
        probs <- c(.025, .975, .16, .84)
        names(probs) <- c("q025", "q975", "q16", "q84")
    } else if (is.null(attr(probs, "names"))) {
        names(probs) <- paste0("q", probs)
    }
    
    J <- length(traitItem)
    M <- mcmc.objects %>% magrittr::extract2(1) %>% nrow
    chains <- mcmc.objects %>% length
    
    arsModel <- switch(fitModel, 
                       "ext"  = TRUE,  
                       "ext2" = TRUE,
                       "ext3" = TRUE,
                       "ext4" = TRUE,
                       "ext5" = TRUE,
                       "2012" = FALSE,
                       "pcm"  = FALSE)
    
    # tmp1 <- M %>% magrittr::divide_by(iter) %>% floor
    # reps <- sapply(seq(1, by = M, length = chains), . %>% seq(by = tmp1, length = iter)) %>% as.vector
    
    thin <- M %>% magrittr::divide_by(iter) %>% round

    dimen <- switch(fitModel, 
                    "ext" = 3,  
                    "ext2" = 3,
                    "ext3" = 4,
                    "ext4" = 3,
                    "ext5" = 3,
                    "2012" = 2,
                    "pcm" = 0) + 1
    n.trait <- length(table(traitItem))
    S <- dimen - 1 + n.trait
    
    vars <- c("^Sigma", "^beta")
    if (arsModel == TRUE)
        vars <- c(vars, "^beta_ARS_extreme")
    # vars <- switch(fitModel, 
    #                 "ext"  = c(vars, "^beta_ARS_extreme"),  
    #                 "ext2" = c(vars, "^beta_ARS_extreme"),
    #                 "ext3" = c(vars, "^beta_ARS_extreme"),
    #                 "ext4" = c(vars, "^beta_ARS_extreme"),
    #                 "ext5" = c(vars, "^beta_ARS_extreme"),
    #                 "2012" = vars,
    #                 "pcm"  = vars)

    # fit_mcmc2 <- coda:::window.mcmc.list(mcmc.objects, thin = thin) %>%
    fit_mcmc2 <- window(mcmc.objects, thin = thin) %>% 
        runjags::combine.mcmc(., vars = vars)
    
    # fit_mcmc2 <- runjags::combine.mcmc(mcmc.objects, vars = vars)
    rm(mcmc.objects)
    
    reps <- seq_len(nrow(fit_mcmc2))
    
    Sigma1 <- runjags::combine.mcmc(fit_mcmc2, vars = "^Sigma", collapse.chains = F)
    
    betas <- runjags::combine.mcmc(fit_mcmc2, vars = "^beta", collapse.chains = F) %>% 
        tibble::as_tibble() %>% 
        dplyr::select(., dplyr::matches("beta\\[")) %>% 
        dplyr::select(sapply(1:ifelse(fitModel == "pcm", 4, dimen), . %>%
                                 seq(., by = ifelse(fitModel == "pcm", 4, dimen), length = J)) %>%
                          as.vector) %>%
        t %>% 
        matrix %>% 
        array(dim = c(J, ifelse(fitModel == "pcm", 4, dimen), length(reps)))
    
    if (arsModel == TRUE)
        beta_ARS_extreme <- runjags::combine.mcmc(fit_mcmc2, vars = "^beta_ARS_extreme",
                                                  collapse.chains = F)
    
    pp <- array(NA_real_, dim = c(N, J, length(reps)))
    
    message(paste0("I'm going to loop over ", length(reps),
                   " replications, this may take a little time, go get a coffee."))
    
    pb <- txtProgressBar(style = 3, char = "zzz ", min = min(reps), max = max(reps))
    
    if (fitModel != "pcm") {
        for(rrr in reps) {
            
            theta <- Sigma1[rrr, ] %>%
                matrix(., S, S) %>% 
                MASS::mvrnorm(n = N, mu = rep(0, S), Sigma = .)
            
            # trait <- acquies <- extreme <- middle <- matrix(NA_real_, N, J)
            p_cat <- array(NA_real_, dim = c(N, J, 5))
            # extreme_a <- rep(NA_real_, N)
            
            middle <- matrix(theta[, 1], N, J) %>% 
                magrittr::subtract(matrix(betas[, 1, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            extreme <- matrix(theta[, 2], N, J) %>% 
                magrittr::subtract(matrix(betas[, 2, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            if (arsModel == TRUE) {
                acquies <- matrix(theta[, 3], N, J) %>%
                    magrittr::subtract(matrix(betas[, 3, rrr], N, J, byrow = T)) %>% 
                    pnorm
            } else {
                acquies <- matrix(0, N, J)
            }
            
            trait <- theta[, (traitItem + dimen - 1)] %>% 
                magrittr::subtract(matrix(betas[, dimen, rrr], N, J, byrow = T)) %>% 
                magrittr::multiply_by(matrix((-1)^revItem, N, J, byrow = T)) %>% 
                pnorm
            
            if (arsModel == TRUE) {
                extreme_a <- matrix(theta[, 2], N, J) %>% 
                    magrittr::subtract(matrix(beta_ARS_extreme[rrr, ], N, J)) %>% 
                    pnorm
            } else {
                extreme_a <- matrix(0, N, J)
            }
            
            p_cat[, , 1] <- (1-acquies)*(1-middle)*(1-trait)*extreme
            p_cat[, , 2] <- (1-acquies)*(1-middle)*(1-trait)*(1-extreme)
            p_cat[, , 3] <- (1-acquies)*   middle
            p_cat[, , 4] <- (1-acquies)*(1-middle)*trait*(1-extreme) + acquies*(1-extreme_a)
            p_cat[, , 5] <- (1-acquies)*(1-middle)*trait*   extreme  + acquies*   extreme_a
            
            # pp[, , which(reps == rrr)] <- p_cat %>%
            pp[, , rrr] <- p_cat %>%
                apply(1:2, function(x) rmultinom(1, 1, x)) %>%
                magrittr::equals(1) %>% 
                apply(2:3, which)
            
            setTxtProgressBar(pb, rrr)
            
            # for(nnn in 1:N){	
            #     
            #     extreme_a[nnn] <- pnorm(theta[nnn,2] - beta_ARS_extreme[rrr, ])
            #     
            #     # loop across items
            #     for(jjj in 1:J){
            #         middle[nnn,jjj]  <- pnorm(theta[nnn,1]-beta[jjj,1,rrr])   #pnorm()
            #         extreme[nnn,jjj] <- pnorm(theta[nnn,2]-beta[jjj,2,rrr])
            #         acquies[nnn,jjj] <- pnorm(theta[nnn,3]-beta[jjj,3,rrr])
            #         # standard items: theta-beta  // reversed items: beta-theta
            #         trait[nnn,jjj]   <- pnorm( (-1)^revItem[jjj] *(theta[nnn,3+traitItem[jjj]]-beta[jjj,4,rrr]) )
            #         
            #         # response probabilities: MPT model for response categories 1, 2, 3, 4, 5
            #         p_cat[nnn,jjj,1] <- (1-acquies[nnn,jjj])*(1-middle[nnn,jjj])*(1-trait[nnn,jjj])*extreme[nnn,jjj]
            #         p_cat[nnn,jjj,2] <- (1-acquies[nnn,jjj])*(1-middle[nnn,jjj])*(1-trait[nnn,jjj])*(1-extreme[nnn,jjj])
            #         p_cat[nnn,jjj,3] <- (1-acquies[nnn,jjj])*   middle[nnn,jjj]
            #         p_cat[nnn,jjj,4] <- (1-acquies[nnn,jjj])*(1-middle[nnn,jjj])*trait[nnn,jjj]*(1-extreme[nnn,jjj]) + acquies[nnn,jjj]*(1-extreme_a[nnn])
            #         p_cat[nnn,jjj,5] <- (1-acquies[nnn,jjj])*(1-middle[nnn,jjj])*trait[nnn,jjj]*extreme[nnn,jjj]     + acquies[nnn,jjj]*   extreme_a[nnn]
            #         
            #         pp[nnn, jjj, which(reps == rrr)] <- rmultinom(1, 1, p_cat[nnn, jjj, ]) %>% equals(1) %>% which
            #     }
            # }
        }
    } else {
        for(rrr in reps) {
            
            theta <- Sigma1[rrr, ] %>%
                matrix(., S, S) %>% 
                MASS::mvrnorm(n = N, mu = rep(0, S), Sigma = .)
            
            # p_cat <- array(NA_real_, dim = c(N, J, 5))
            
            theta_x <- theta[, rep((traitItem + dimen - 1), each = 4)]
            
            betas_x <- betas[, , rrr] %>% 
                t %>% 
                matrix(nrow = N, ncol = J*4, byrow = T)
            
            log_arg <- (theta_x - betas_x) %>% 
                t %>% 
                matrix(ncol = 4, byrow = T)
            p1 <- log_arg %>% 
                cbind(1, .) %>% 
                apply(1, . %>% cumsum %>% exp) %>%
                t %>% 
                # exp %>% 
                magrittr::divide_by(., rowSums(.))
            
            p1[revItem == 1, ] <- p1[revItem == 1, ] %>% 
                apply(1, rev) %>% 
                t
            
            p_cat <- p1 %>%
                array(dim = c(J, N, 5)) %>%
                apply(c(1, 3), t)
            
            pp[, , rrr] <- p_cat %>%
                apply(1:2, function(x) rmultinom(1, 1, x)) %>%
                magrittr::equals(1) %>% 
                apply(2:3, which)
            
            setTxtProgressBar(pb, rrr)
        }
    }
    
    close(pb)
    
    message("I drew posterior predictives, I need a little bit more time to restructure the results, I ask for your patience.")
    
    tmp1 <- pp %>% 
        reshape2::melt() %>%
        dplyr::transmute(Person = factor(.data$Var1),
                         Item = factor(.data$Var2),
                         Iter = .data$Var3,
                         value = .data$value)
    
    tmp2 <- tmp1 %>% 
        dplyr::group_by_("Item", "Iter") %>% 
        dplyr::summarise(Categ = list(.data$value %>%
                                          factor(levels = 1:5) %>%
                                          table %>%
                                          names %>%
                                          factor),
                  pp = list(.data$value %>%
                                factor(levels = 1:5) %>%
                                table %>%
                                prop.table)) %>%
        # summarise(pp = list(quantile(value, probs = c(.025, .975, .16, .84))),
        #           Q = list(c("q025", "q975", "q16", "q84"))) %>%
        tidyr::unnest()
    
    tmp3 <- tmp2 %>% 
        dplyr::group_by_("Item", "Categ") %>% 
        dplyr::summarise(pp = list(pp %>%
                                       quantile(., probs = probs)),
                         Q = list(names(probs))) %>% 
        tidyr::unnest()
    
    # tmp3$Item <- factor(tmp3$Item, labels = levels(dat_obs1$Item))
    
    tmp4 <- tmp3 %>% 
        reshape2::dcast(value.var = "pp", formula = Item + Categ ~ Q)
    
    return(tmp4)
}
