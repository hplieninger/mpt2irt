if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

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
#' @param fit_sum List. Output from \code{\link{summarize_irtree_fit}} that
#'   contains \code{mcmc.objects}.
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
# @importFrom magrittr %>%
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
#' res3 <- pp_irtree(res2, iter = 10, N = 512)
#' }
#' @export
pp_irtree <- function(fit_sum = NULL,
                      mcmc.objects = NULL,
                      iter = 100,
                      probs = NULL,
                      N = NULL,
                      traitItem = NULL,
                      revItem = NULL,
                      fitModel = NULL) {
    
    checkmate::assert_list(fit_sum, all.missing = FALSE, min.len = 1,
                           null.ok = TRUE)
    checkmate::assert_list(mcmc.objects, types = "numeric", any.missing = FALSE,
                           min.len = 1, null.ok = TRUE)
    if (is.null(fit_sum) & is.null(mcmc.objects)) 
        stop("At least one of 'fit_sum' and 'mcmc.objects' must be specified.")
    if (is.null(mcmc.objects)) {
        if ("mcmc" %in% names(fit_sum)) {
            mcmc.objects <- fit_sum$mcmc
        } else {
            stop("'fit_sum' must contain an element 'mcmc'; or 'mcmc.objects' ",
                 "must be specified.")
        }
        # flag1 <- ifelse("args" %in% names(fit_sum), TRUE, FALSE)
    } # else {
    #     flag1 <- FALSE
    # }
    
    checkmate::assert_int(iter, lower = 1, upper = nrow(mcmc.objects[[1]]))
    checkmate::assert_numeric(probs, lower = 0, upper = 1,
                              any.missing = FALSE,  min.len = 1, unique = TRUE,
                              null.ok = TRUE)
    checkmate::qassert(N, "X1[1,)")
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 min.len = 1, null.ok = !is.null(fit_sum$args$revItem))
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 min.len = 1, null.ok = !is.null(fit_sum$args$traitItem))
    checkmate::assert_character(fitModel, min.chars = 1, any.missing = FALSE,
                                len = 1, null.ok = !is.null(fit_sum$args$fitModel))
    
    if (is.null(revItem)) revItem <- fit_sum$args$revItem
    if (is.null(traitItem)) traitItem <- fit_sum$args$traitItem
    if (is.null(fitModel)) fitModel <- fit_sum$args$fitModel
    # if (fitModel == "steps") stop("fitModel 'steps' not implemented yet.")
    
    if (is.null(probs)) {
        probs <- c(.025, .975, .16, .84, .50)
        names(probs) <- c("q025", "q975", "q16", "q84", "q50")
    } else if (is.null(attr(probs, "names"))) {
        names(probs) <- paste0("q", probs)
    }
    
    J <- length(traitItem)
    M <- nrow(mcmc.objects[[1]])
    # chains <- length(mcmc.objects)
    thin <- round(M/iter)
    
    arsModel <- switch(fitModel, 
                       "ext"   = TRUE,  
                       "ext2"  = TRUE,
                       # "ext3"  = TRUE,
                       # "ext4"  = TRUE,
                       # "ext5"  = TRUE,
                       "2012"  = FALSE,
                       "pcm"   = FALSE,
                       "steps" = FALSE,
                       "shift" = FALSE)
    
    # tmp1 <- M %>% magrittr::divide_by(iter) %>% floor
    # reps <- sapply(seq(1, by = M, length = chains), . %>% seq(by = tmp1, length = iter)) %>% as.vector

    dimen <- switch(fitModel, 
                    "ext"   = 3,  
                    "ext2"  = 3,
                    # "ext3"  = 4,
                    # "ext4"  = 3,
                    # "ext5"  = 3,
                    "2012"  = 2,
                    "pcm"   = 0,
                    "steps" = 0,
                    "shift" = 3) + 1
    if (is.null(fit_sum$args$S)) {
        n.trait <- length(table(traitItem))
        S <- dimen - 1 + n.trait
    } else {
        S <- fit_sum$args$S
    }
    
    # S_t:  Number of theta dimensions
    # S_b1: Number of columns in beta matrix
    # S_b2: Number of beta dimensions

    S_b1 <- switch(fitModel,
                   "2012"  = 3,
                   "ext"   = 4,
                   "ext2"  = 5,
                   "pcm"   = 4,
                   "steps" = 4,
                   "shift" = 3,
                   S)

    vars <- c("^Sigma", "^beta")
    if (fitModel %in% c("ext")) {
        vars <- c(vars, "^beta_ARS_extreme")
    }

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
        dplyr::select(sapply(1:S_b1, . %>% seq(., by = S_b1, length = J)) %>%
                          as.vector) %>%
        t %>% 
        matrix %>% 
        array(dim = c(J, S_b1, length(reps)))
    
    if (fitModel %in% c("ext")) {
        beta_ARS_extreme <- runjags::combine.mcmc(fit_mcmc2, vars = "^beta_ARS_extreme",
                                                  collapse.chains = F)
    }
    
    pp <- array(NA_real_, dim = c(N, J, length(reps)))
    
    message(paste0("I'm going to loop over ", length(reps),
                   " replications, this may take a little time, go get a coffee."))
    
    # If user has dplyr installed, use dplyr::progress_estimated(), otherwise
    # use txtProgressBar()
    p <- FALSE
    try(p <- dplyr::progress_estimated(length(reps), 5), silent = T)
    if (!is.environment(p)) {
        pb <- txtProgressBar(style = 3, char = "zzz ", min = min(reps), max = max(reps))
    }
    # pb <- txtProgressBar(style = 3, char = "zzz ", min = min(reps), max = max(reps))
    
    ##### pp for 'ext' and '2012' #####
    if (fitModel %in% c("ext", "2012", "ext2")) {
        for (rrr in reps) {
            
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
            
            if (fitModel == "ext") {
                extreme_a <- matrix(theta[, 2], N, J) %>% 
                    magrittr::subtract(matrix(beta_ARS_extreme[rrr, ], N, J)) %>% 
                    pnorm
            } else if (fitModel == "ext2") {
                extreme_a <- matrix(theta[, 2], N, J) %>% 
                    magrittr::subtract(matrix(betas[, 5, rrr], N, J, byrow = T)) %>% 
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
            
            if (is.environment(p)) {
                p$tick()$print()
            } else {
                setTxtProgressBar(pb, rrr)
            }
        }
    } else if (fitModel == "steps") {
        ##### pp for 'steps' #####
        for (rrr in reps) {
            
            theta <- Sigma1[rrr, ] %>%
                matrix(., S, S) %>% 
                MASS::mvrnorm(n = N, mu = rep(0, S), Sigma = .)
            
            p_catx <- array(NA_real_, dim = c(N, J, 5))
            
            node1 <- matrix(theta[, traitItem], N, J) %>% 
                magrittr::subtract(matrix(betas[, 1, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            node2 <- matrix(theta[, traitItem], N, J) %>% 
                magrittr::subtract(matrix(betas[, 2, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            node3 <- matrix(theta[, traitItem], N, J) %>% 
                magrittr::subtract(matrix(betas[, 3, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            node4 <- matrix(theta[, traitItem], N, J) %>% 
                magrittr::subtract(matrix(betas[, 4, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            p_catx[ , , 1] = (1-node1);
            p_catx[ , , 2] =    node1 *(1-node2);
            p_catx[ , , 3] =    node1 *   node2 *(1-node3);
            p_catx[ , , 4] =    node1 *   node2 *   node3 *(1-node4);
            p_catx[ , , 5] =    node1 *   node2 *   node3 *   node4 ;
            
            p_cat <- p_catx
            
            p_cat[ , revItem == 1, 5] <- p_catx[ , revItem == 1, 1]
            p_cat[ , revItem == 1, 4] <- p_catx[ , revItem == 1, 2]
            p_cat[ , revItem == 1, 3] <- p_catx[ , revItem == 1, 3]
            p_cat[ , revItem == 1, 2] <- p_catx[ , revItem == 1, 4]
            p_cat[ , revItem == 1, 1] <- p_catx[ , revItem == 1, 5]
            
            pp[, , rrr] <- p_cat %>%
                apply(1:2, function(x) rmultinom(1, 1, x)) %>%
                magrittr::equals(1) %>% 
                apply(2:3, which)
            
            if (is.environment(p)) {
                p$tick()$print()
            } else {
                setTxtProgressBar(pb, rrr)
            }
        }
    } else if (fitModel == "pcm") {
        ##### pp for 'pcm' #####
        for (rrr in reps) {
            
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
            
            if (is.environment(p)) {
                p$tick()$print()
            } else {
                setTxtProgressBar(pb, rrr)
            }
        }
    } else if (fitModel == "shift") {
        ##### pp for 'shift' #####
        for (rrr in reps) {
            
            theta <- Sigma1[rrr, ] %>%
                matrix(., S, S) %>% 
                MASS::mvrnorm(n = N, mu = rep(0, S), Sigma = .)
            
            p_cat <- array(NA_real_, dim = c(N, J, 5))
            
            middle <- matrix(theta[, 1], N, J) %>% 
                magrittr::subtract(matrix(betas[, 1, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            extreme <- matrix(theta[, 2], N, J) %>% 
                magrittr::subtract(matrix(betas[, 2, rrr], N, J, byrow = T)) %>% 
                pnorm
            
            # acquies <- matrix(theta[, 3], N, J) %>%
            #     magrittr::subtract(matrix(betas[, 3, rrr], N, J, byrow = T)) %>% 
            #     pnorm
            
            trait <- theta[, (traitItem + dimen - 1)] %>% 
                magrittr::subtract(matrix(betas[, 3, rrr], N, J, byrow = T)) %>% 
                magrittr::multiply_by(matrix((-1)^revItem, N, J, byrow = T)) %>% 
                magrittr::add(matrix(theta[, 3], N, J)) %>% 
                pnorm
            
            p_cat[, , 1] <- (1-middle)*(1-trait)*extreme
            p_cat[, , 2] <- (1-middle)*(1-trait)*(1-extreme)
            p_cat[, , 3] <-    middle
            p_cat[, , 4] <- (1-middle)*trait*(1-extreme)
            p_cat[, , 5] <- (1-middle)*trait*   extreme 
            
            # pp[, , which(reps == rrr)] <- p_cat %>%
            pp[, , rrr] <- p_cat %>%
                apply(1:2, function(x) rmultinom(1, 1, x)) %>%
                magrittr::equals(1) %>% 
                apply(2:3, which)
            
            if (is.environment(p)) {
                p$tick()$print()
            } else {
                setTxtProgressBar(pb, rrr)
            }
        }
    }
    
    if (!is.environment(p)) close(pb)
    
    message("I drew posterior predictives, I need a little bit more time to ",
            "restructure the results, I ask for your patience.")
    
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
                                prop.table),
                  Persons = dplyr::n_distinct(Person)) %>%
        # summarise(pp = list(quantile(value, probs = c(.025, .975, .16, .84))),
        #           Q = list(c("q025", "q975", "q16", "q84"))) %>%
        tidyr::unnest()
    
    tmp3 <- tmp2 %>% 
        dplyr::group_by_("Item", "Categ") %>% 
        dplyr::summarise(pp = list(pp %>%
                                       quantile(., probs = probs)),
                         Q = list(names(probs)),
                         Samples = dplyr::n_distinct(Iter),
                         Persons = dplyr::first(Persons)) %>% 
        tidyr::unnest()
    
    # tmp3$Item <- factor(tmp3$Item, labels = levels(dat_obs1$Item))
    
    tmp4 <- tmp3 %>% 
        reshape2::dcast(value.var = "pp",
                        formula = Item + Categ + Persons + Samples ~ Q)
    
    return(tmp4)
}
