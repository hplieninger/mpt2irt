if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Calculate posterior predictive probabilities.
#' 
#' Function takes an MCMC list of posterior samples and calculates the
#' model-predicted probabilities. This can either be done for the persons in the
#' sample or for out-of-sample predictions for new persons (\code{new_theta =
#' TRUE}).
#' 
#' @param fit_sum List. Output from \code{\link{summarize_irtree_fit}} that
#'   contains \code{mcmc.objects}.
#' @param iter Numeric. Number of iterations to use, the maximum is the total
#'   number of retained iterations (via \code{\link{fit_irtree}}).
#' @param N Numeric. Number of persons for whom posterior predictives should be
#'   drawn. Should be equal to the number of persons in the sample if
#'   \code{new_theta = FALSE}.
#' @param new_theta Logical. Wheter to calculate the probabilities for the
#'   persons in the sample or for out-of-sample predictions for \code{N} 'new'
#'   persons (\code{new_theta = TRUE}).
#' @return Returns an array of probabilities of dimension \code{iter} x \code{N}
#'   x J x 5 (for J items with 5 response categories).
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
#' # posterior predictive probabilities
#' res3 <- post_prob_irtree(res2)
#' dim(res3)
#' }
#' @export
post_prob_irtree <- function(fit_sum = NULL,
                             mcmc.objects = NULL,
                             iter = 100,
                             N = NULL,
                             traitItem = NULL,
                             revItem = NULL, 
                             fitModel = NULL,
                             new_theta = FALSE) {
    
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
            stop("'fit_sum' must contain an element 'mcmc'. Or 'mcmc.objects' ",
                 "must be specified.")
        }
    }
    
    checkmate::assert_int(iter, lower = 1, upper = nrow(mcmc.objects[[1]])*length(mcmc.objects))
    # checkmate::qassert(N, "X1[1,)")
    checkmate::assert_int(N, lower = 1, null.ok = !is.null(fit_sum$args$N))
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 min.len = 1, null.ok = !is.null(fit_sum$args$revItem))
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 min.len = 1, null.ok = !is.null(fit_sum$args$traitItem))
    checkmate::assert_character(fitModel, min.chars = 1, any.missing = FALSE,
                                len = 1, null.ok = !is.null(fit_sum$args$fitModel))
    checkmate::qassert(new_theta, "B1")
    
    # if (is.null(revItem)) revItem <- fit_sum$revItem
    # if (is.null(traitItem)) traitItem <- fit_sum$traitItem
    # if (is.null(fitModel)) fitModel <- fit_sum$fitModel
    # if (is.null(N)) N <- fit_sum$N
    if (is.null(revItem)) revItem <- fit_sum$args$revItem
    if (is.null(traitItem)) traitItem <- fit_sum$args$traitItem
    if (is.null(fitModel)) fitModel <- fit_sum$args$fitModel
    if (is.null(N)) N <- fit_sum$args$N
    
    J <- length(traitItem)
    M <- nrow(mcmc.objects[[1]])
    
    mcmc_attr <- attr(mcmc.objects[[1]], "mcpar")
    
    mcmc.objects <- coda::as.mcmc.list(
        lapply(mcmc.objects, coda::mcmc, start = mcmc_attr[1], thin = 1))
    
    thin <- floor(M/iter*length(mcmc.objects))
    
    arsModel <- fitModel %in% c("ext", "ext2")
    
    dimen <- 1 + switch(fitModel, 
                        "ext"   = 3,  
                        "ext2"  = 3,
                        "2012"  = 2,
                        "pcm"   = 0,
                        "steps" = 0,
                        "shift" = 3)
    if (is.null(fit_sum$args$S)) {
        n.trait <- length(table(traitItem))
        S <- dimen - 1 + n.trait
    } else {
        S <- fit_sum$args$S
    }
    
    # S_t = S:  Number of theta dimensions
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
    
    vars <- c("^Sigma", "^beta", "^theta")
    if (fitModel %in% c("ext")) {
        vars <- c(vars, "^beta_ARS_extreme")
    }
    
    fit_mcmc2 <- window(mcmc.objects, thin = thin) %>% 
        runjags::combine.mcmc(., vars = vars)
    rm(mcmc.objects)
    
    reps <- seq_len(nrow(fit_mcmc2))
    R <- length(reps)
    
    # beta <- runjags::combine.mcmc(fit_mcmc2, vars = "^beta", collapse.chains = F) %>% 
    #     tibble::as_tibble() %>% 
    #     dplyr::select(., dplyr::matches("beta\\[")) %>% 
    #     dplyr::select(sapply(1:S_b1, . %>% seq(., by = S_b1, length = J)) %>%
    #                       as.vector) %>%
    #     t %>% 
    #     matrix %>% 
    #     array(dim = c(J, S_b1, length(reps)))
    
    tmp1 <- runjags::combine.mcmc(fit_mcmc2, vars = "^beta", collapse.chains = F) %>% 
        as.data.frame %>% 
        dplyr::select(., dplyr::matches("beta\\[")) %>% 
        {suppressMessages(reshape2::melt(.))}
    tmp2 <- regmatches(tmp1$variable, regexpr("[0-9]+,[0-9]+", tmp1$variable)) %>%
        strsplit(",") %>%
        unlist %>% as.numeric %>% matrix(ncol = 2, byrow = T) %>% 
        data.frame(tmp1, .)
    tmp3 <- tmp2[order(tmp2$X1, tmp2$X2), ]
    
    beta <- apply(
        array(
            matrix(tmp3$value, ncol = length(reps), byrow = T),
            dim = c(S_b1, J, length(reps))),
        c(1, 3), t)
    
    if (new_theta == TRUE) {
        message("New person parameters are sampled from covariance matrix Sigma.")
        Sigma <- runjags::combine.mcmc(fit_mcmc2, vars = "^Sigma", collapse.chains = F)
        theta <- array(NA, c(N, S, R))
        for (rrr in reps)
            theta[, , rrr] <-  MASS::mvrnorm(n = N, mu = rep(0, S), 
                                           Sigma = matrix(Sigma[rrr, ], S, S))
        
    } else {
        # # sel_theta <- grepl("theta\\[[0-9]+,[0-9]+\\]", colnames(fit_mcmc2))
        # theta <- runjags::combine.mcmc(fit_mcmc2, vars = "^theta", collapse.chains = F) %>%
        #     tibble::as_tibble() %>% 
        #     dplyr::select(., dplyr::matches("theta\\[")) %>% 
        #     dplyr::select(sapply(1:S, . %>% seq(., by = S, length = N)) %>%
        #                       as.vector) %>%
        #     t %>% 
        #     matrix %>% 
        #     array(dim = c(N, S, R))
        
        tmp1 <- runjags::combine.mcmc(fit_mcmc2, vars = "^theta", collapse.chains = F) %>% 
            as.data.frame %>% 
            dplyr::select(., dplyr::matches("theta\\[")) %>% 
            {suppressMessages(reshape2::melt(.))}
        tmp2 <- regmatches(tmp1$variable, regexpr("[0-9]+,[0-9]+", tmp1$variable)) %>%
            strsplit(",") %>%
            unlist %>% as.numeric %>% matrix(ncol = 2, byrow = T) %>% 
            data.frame(tmp1, .)
        tmp3 <- tmp2[order(tmp2$X1, tmp2$X2), ]
        
        theta <- apply(
            array(
                matrix(tmp3$value, ncol = length(reps), byrow = T),
                dim = c(S, N, length(reps))),
            c(1, 3), t)    
    }
    
    if (fitModel %in% c("ext")) {
        beta_ARS_extreme <- runjags::combine.mcmc(fit_mcmc2, vars = "^beta_ARS_extreme",
                                                  collapse.chains = F)
    }
    
    # If user has dplyr installed, use dplyr::progress_estimated(), otherwise
    # use txtProgressBar()
    p <- FALSE
    try(p <- dplyr::progress_estimated(R, 5), silent = TRUE)
    if (!is.environment(p)) {
        pb <- txtProgressBar(style = 3, char = "zzz ", min = min(reps), max = max(reps))
    }
    
    prob <- array(NA_real_, dim = c(length(reps), N, J, 5))
    for (rrr in reps) {
        
        if (fitModel %in% c("ext", "2012", "ext2")) {
            ##### probs for 'ext' and '2012' #####
            
            middle <- pnorm(outer(theta[, 1, rrr], beta[, 1, rrr], "-"))
            extreme <- pnorm(outer(theta[, 2, rrr], beta[, 2, rrr], "-"))
            acquies <- matrix(0, N, J)
            if (arsModel == TRUE)
                acquies <- pnorm(outer(theta[, 3, rrr], beta[, 3, rrr], "-"))
            
            latent <- theta[, traitItem + dimen - 1, rrr] - 
                matrix(beta[, dimen, rrr], N, J, byrow = TRUE)
            latent[, revItem == 1] <- -latent[, revItem == 1]
            trait <- pnorm(latent)
            
            extreme_a <- matrix(0, N, J)
            if (fitModel == "ext") {
                extreme_a <- matrix(pnorm(theta[, 2, rrr] - beta_ARS_extreme[rrr, ]), N, J)
            }
            if (fitModel == "ext2") {
                extreme_a <- pnorm(outer(theta[, 2, rrr], beta[, 5, rrr], "-"))
            }
            
            prob[rrr, , , 1] <- (1-acquies)*(1-middle)*(1-trait)*extreme
            prob[rrr, , , 2] <- (1-acquies)*(1-middle)*(1-trait)*(1-extreme)
            prob[rrr, , , 3] <- (1-acquies)*   middle
            prob[rrr, , , 4] <- (1-acquies)*(1-middle)*trait*(1-extreme) + acquies*(1-extreme_a)
            prob[rrr, , , 5] <- (1-acquies)*(1-middle)*trait*   extreme  + acquies*   extreme_a
            
            
        } else if (fitModel == "steps") {
            ##### probs for 'steps' #####
            node1 <- pnorm(theta[, traitItem, rrr] - 
                               matrix(beta[, 1, rrr], N, J, byrow = TRUE))
            node2 <- pnorm(theta[, traitItem, rrr] - 
                               matrix(beta[, 2, rrr], N, J, byrow = TRUE))
            node3 <- pnorm(theta[, traitItem, rrr] - 
                               matrix(beta[, 3, rrr], N, J, byrow = TRUE))
            node4 <- pnorm(theta[, traitItem, rrr] - 
                               matrix(beta[, 4, rrr], N, J, byrow = TRUE))
            
            p_catx <- array(NA_real_, dim = c(N, J, 5))
            p_catx[ , , 1] = (1-node1);
            p_catx[ , , 2] =    node1 *(1-node2);
            p_catx[ , , 3] =    node1 *   node2 *(1-node3);
            p_catx[ , , 4] =    node1 *   node2 *   node3 *(1-node4);
            p_catx[ , , 5] =    node1 *   node2 *   node3 *   node4 ;
            
            prob[rrr, , ,] <- p_catx
            
            prob[rrr, , revItem == 1, 5] <- p_catx[ , revItem == 1, 1]
            prob[rrr, , revItem == 1, 4] <- p_catx[ , revItem == 1, 2]
            prob[rrr, , revItem == 1, 3] <- p_catx[ , revItem == 1, 3]
            prob[rrr, , revItem == 1, 2] <- p_catx[ , revItem == 1, 4]
            prob[rrr, , revItem == 1, 1] <- p_catx[ , revItem == 1, 5]
            
            
        } else if (fitModel == "pcm") {
            ##### probs for 'pcm' #####
            theta_x <- theta[, rep((traitItem + dimen - 1), each = 4),rrr]
            
            betas_x <- beta[, , rrr] %>% 
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
            
            prob[rrr, , ,] <- p1 %>%
                array(dim = c(J, N, 5)) %>%
                apply(c(1, 3), t)
            
            
        } else if (fitModel == "shift") {
            ##### probs for 'shift' #####
            
            middle <- pnorm(outer(theta[, 1, rrr], beta[, 1, rrr], "-"))
            extreme <- pnorm(outer(theta[, 2, rrr], beta[, 2, rrr], "-"))
            
            trait <- theta[, (traitItem + dimen - 1),rrr] %>% 
                magrittr::subtract(matrix(beta[, 3, rrr], N, J, byrow = T)) %>% 
                magrittr::multiply_by(matrix((-1)^revItem, N, J, byrow = T)) %>% 
                magrittr::add(matrix(theta[, 3,rrr], N, J)) %>% 
                pnorm
            
            prob[rrr, , , 1] <- (1-middle)*(1-trait)*extreme
            prob[rrr, , , 2] <- (1-middle)*(1-trait)*(1-extreme)
            prob[rrr, , , 3] <-    middle
            prob[rrr, , , 4] <- (1-middle)*trait*(1-extreme)
            prob[rrr, , , 5] <- (1-middle)*trait*   extreme 
            
        }
        if (is.environment(p)) {
            p$tick()$print()
        } else {
            setTxtProgressBar(pb, rrr)
        }
    }
    return(prob)
}

#' Calculate expected responses on the basis of posterior predicitive
#' probabilities.
#'
#' The function takes the posterior predicitve probabilities (from
#' \code{\link{post_prob_irtree}}) and calculates for every
#' iteration-person-item combination the expected response. This is needed for
#' Yen's Q3 (see \code{\link{statistic_Q3}}).
#'
#' @param prob Numeric array of dimension R x N x J x 5 (for N persons, J items
#'   with 5 categories, and R iterations).
#' @return Array of expected responses of dimension R x N x J (for N persons, J
#'   items, and R iterations).
#'
# @export
posterior_expected <- function(prob = NULL) {
    # checkmate::assert_array(prob, mode = "numeric", any.missing = FALSE, d = 4)
    # checkmate::qassert(prob, "N+[0,1]")
    # 
    # # apply is slow here:
    # apply(prob, MARGIN = 1:3, function(p) sum(p * (1:5))) - 
    d <- dim(prob)
    r <- array(matrix(prob, ncol = 5) %*% 1:5, d[1:3])
    return(r)
}

#' Draw responses on the basis of posterior predicitive probabilities.
#'
#' The function takes the posterior predicitve probabilities (from
#' \code{\link{post_prob_irtree}}) and draws for every
#' iteration-person-item combination a predicted response.
#'
#' @param prob Numeric array of dimension R x N x J x 5 (for N persons, J items
#'   with 5 categories, and R iterations).
#' @return Array of predicted responses of dimension R x N x J (for N persons, J
#'   items, and R iterations).
#'
# @export
posterior_predicted <- function(prob = NULL) {
    checkmate::assert_array(prob, mode = "numeric", any.missing = FALSE, d = 4)
    checkmate::qassert(prob, "N+[0,1]")
    
    ### apply is unefficient here:
    # apply(prob, MARGIN = 1:3, function(p) sample.int(5, 1, prob = p))
    d <- dim(prob)
    array(sample_pp(matrix(prob, ncol = 5)), d[1:3])
    # pp <- array(NA, dim(prob)[1:3])
    # for (r in 1:dim(pp)[1]){
    #     for (n in 1:dim(pp)[2]){
    #         for (j in 1:dim(pp)[3]){
    #             pp[r,n,j] <- sample.int(5, 1, prob = prob[r,n,j,])
    #         }
    #     }
    # }
    # pp
}

##### low-level implementation of statistics (input: 1 replication) #####
# 
# resp: N x J matrix      [observed / predicted responses]
# exp: N x J matrix       [expected response]
# prob: N x J x 5 array   [predicted probabilities]

#' Discrepency measure for PPC: Item-total correlation.
#'
#' This function calculates the polyserial correlation between a person's total
#' score and his or her item response. This discrepency measure is useful to
#' detect misfit due to a missing item-slope/-discrimination parameter. It may
#' not be suitable for the response style models discussed herein, because these
#' models assume that the total score as well as the item response is a
#' composite of target trait and response styles.
#' 
#' @references Li, T., Xie, C., & Jiao, H. (2017). Assessing fit of alternative
#'   unidimensional polytomous IRT models using posterior predictive model
#'   checking. Psychological Methods, 22, 397-408. doi:10.1037/met0000082
#' @references Zhu, X., & Stone, C. A. (2012). Bayesian comparison of
#'   alternative graded response models for performance assessment applications.
#'   Educational and Psychological Measurement, 72, 774-799.
#'   doi:10.1177/0013164411434638
#'
#' @param resp Numeric matrix of dimension N x J (for N persons and J items)
#'   containing the observed item responses.
#' @inheritParams fit_irtree
#' @return Vector of length J of polyserial correlations.
#'
# @export
statistic_item_cor <- function(resp = NULL,
                               revItem = NULL,
                               traitItem = NULL) {
    # if (is.null(revItem)) revItem <- rep(0, ncol(resp))
    # if (is.null(traitItem)) traitItem <- rep(1, ncol(resp))
    
    resp[, revItem == 1] <- 6 - resp[, revItem == 1]
    scores <- vapply(1:6, function(x) rowSums(resp[, traitItem == x]),
                     vector("numeric", length = nrow(resp)))
    # scores <- rowSums(resp)
    # apply(resp, 2, function(r) polycor::polyserial(scores, r))
    s <-  vapply(1:ncol(resp),
                 function(r) polycor::polyserial(scores[, traitItem[r]], resp[, r]), 1)
    return(s)
}

#' Discrepency measure for PPC: Yen's Q3.
#'
#' This function calculates Yen's Q3 statistic.
#' 
#' @references Li, T., Xie, C., & Jiao, H. (2017). Assessing fit of alternative
#'   unidimensional polytomous IRT models using posterior predictive model
#'   checking. Psychological Methods, 22, 397-408. doi:10.1037/met0000082
#' @references Yen, W. M. (1993). Scaling performance assessments: Strategies
#'   for managing local item dependence. Journal of Educational Measurement, 30,
#'   187-213. doi:10.1111/j.1745-3984.1993.tb00423.x
#' @references Zhu, X., & Stone, C. A. (2012). Bayesian comparison of
#'   alternative graded response models for performance assessment applications.
#'   Educational and Psychological Measurement, 72, 774-799.
#'   doi:10.1177/0013164411434638
#'
#' @param resp Numeric matrix of dimension N x J (for N persons and J items)
#'   containing the responses.
#' @param exp Numeric matrix of dimension N x J (for N persons and J items)
#'   containing the expected responses.
# @inheritParams ppc_irtree
#' @return Vector of length J of polyserial correlations.
#'
# @export
statistic_Q3 <- function(resp = NULL,
                         exp = NULL){
    d <- resp - exp
    # # plot correlations of observed/expected/diffs
    # hist(mapply(cor, data.frame(t(resp)), data.frame(t(exp))))  # per person
    # hist(mapply(cor, data.frame(resp), data.frame(exp)))        # per item
    # hist(cor(d)) 
    
    r <- cor(d)
    r1 <- r[lower.tri(r)]
    return(r1)
    
    # GDDM (Levy & Svetina, 2011)
    # MBC <- combn(ncol(d), 2, FUN = function(idx)
    #     mean(d[,idx[1]] * d[,idx[2]]) )
    # GDDM <- mean(abs(MBC))
}

#' Discrepency measure for PPC: Global Odds ratio for two items.
#'
#' This function calculates the (natural logarithm of the) global odds ratio of
#' a pair of items, which are dichotomized beforehand.
#' 
#' @references Li, T., Xie, C., & Jiao, H. (2017). Assessing fit of alternative
#'   unidimensional polytomous IRT models using posterior predictive model
#'   checking. Psychological Methods, 22, 397-408. doi:10.1037/met0000082
#' @references Zhu, X., & Stone, C. A. (2012). Bayesian comparison of
#'   alternative graded response models for performance assessment applications.
#'   Educational and Psychological Measurement, 72, 774-799.
#'   doi:10.1177/0013164411434638
#'
#' @param x Numeric vector containing the responses to the first item.
#' @param y Numeric vector containing the responses to the second item.
#' @param threshold Numeric vector of length 1. The threshold used for dichotomization.
#' @param constant Numeric vector of length 1. A constant added to all four
#'   frequencies in order to avoid dividing by zero.
# @inheritParams ppc_irtree
#' @return The odds ratio.
#'
# @export
global_OR <- function(x,
                      y,
                      threshold = 3,
                      constant = .1) {
    # odds ratio: (as log, numerically more stable and better to plot)
    ln <- log(table(factor(x > threshold, levels = c(FALSE, TRUE)), 
                    factor(y > threshold, levels = c(FALSE, TRUE))) + constant)
    or <- ln[1,1] + ln[2,2] - (ln[1,2] + ln[2,1])
    return(or)
}

#' Discrepency measure for PPC: Global Odds ratio for three or more items.
#'
#' This function calculates the (natural logarithm of the) global odds ratio of
#' pairs of items, which are dichotomized beforehand.
#' 
#' @references Li, T., Xie, C., & Jiao, H. (2017). Assessing fit of alternative
#'   unidimensional polytomous IRT models using posterior predictive model
#'   checking. Psychological Methods, 22, 397-408. doi:10.1037/met0000082
#' @references Zhu, X., & Stone, C. A. (2012). Bayesian comparison of
#'   alternative graded response models for performance assessment applications.
#'   Educational and Psychological Measurement, 72, 774-799.
#'   doi:10.1177/0013164411434638
#'
#' @param resp Numeric matrix of dimension N x J (for N persons and J items)
#'   containing the responses.
# @inheritParams ppc_irtree
#' @return Vector of length 'J choose 2' odds ratios.
#'
# @export
# 
statistic_OR <- function(resp) {
    combn(ncol(resp), 2, FUN = function(idx)
        global_OR(resp[,idx[1]], resp[,idx[2]]) )
}

# group: N vector         [group membership for Yen's Q1, e.g., ability]
# not Q1: deviation between response and expected value != expected frequency!
# statistic_Q1 <- function(resp, exp, group){
#     sum((resp - exp)^2 / exp)
# }

#' Discrepency measure for PPC: Item score distribution X2.
#'
#' This function calculates the item score distribution, i.e., the discrepancy
#' between observed frequencies and posterior predictive frequencies for each
#' item using Pearson's X2 statistic.
#' 
#' @references Zhu, X., & Stone, C. A. (2012). Bayesian comparison of
#'   alternative graded response models for performance assessment applications.
#'   Educational and Psychological Measurement, 72, 774-799.
#'   doi:10.1177/0013164411434638
#'
#' @param resp Numeric matrix of dimension N x J (for N persons and J items)
#'   containing the responses.
#' @inheritParams ppc_irtree
#' @return Vector of length J containing X2 values.
#'
# @export
statistic_item_score <- function(resp = NULL,
                                 prob = NULL){
    # Zhu & Stone (2012, p. 779)
    x1 <- apply(prob, c(2:3), sum)
    x2 <- t(apply(resp, 2, function(x) table(factor(x, 1:5))))
    
    x3 <- rowSums(((x2 - x1)^2)/x1)
    return(x3)
}

##### high-level usage of statistics (input: all replication) #####

#' Posterior predicitve checks and p-values for several discrepency measures.
#'
#' This function takes posterior predictive probabilities (from
#' \code{\link{post_prob_irtree}}), and compares---for several discrepency
#' measures---observations and predictions both descriptively and via posterior
#' predicitive p-values.
#' 
#' @references Levy, R. (2011). Posterior predictive model checking for
#'   conjunctive multidimensionality in item response theory. Journal of
#'   Educational and Behavioral Statistics, 36, 672-694.
#'   doi:10.3102/1076998611410213
#' @references Li, T., Xie, C., & Jiao, H. (2017). Assessing fit of alternative
#'   unidimensional polytomous IRT models using posterior predictive model
#'   checking. Psychological Methods, 22, 397-408. doi:10.1037/met0000082
#' @references Sinharay, S., Johnson, M. S., & Stern, H. S. (2006). Posterior
#'   predictive assessment of item response theory models. Applied Psychological
#'   Measurement, 30, 298-321. doi:10.1177/0146621605285517
#' @references Zhu, X., & Stone, C. A. (2012). Bayesian comparison of
#'   alternative graded response models for performance assessment applications.
#'   Educational and Psychological Measurement, 72, 774-799.
#'   doi:10.1177/0013164411434638
#'
#' @param prob Numeric array of dimension R x N x J x 5 (for N persons, J items
#'   with 5 categories, and R iterations).
#' @param statistics Character vector of length >= 1 specifying the discrepency
#'   measures to use.
#' @param X Numeric matrix of dimension N x J containing the observed item responses.
#' 
#' \itemize{
#' \item \code{item_cor}: Item-total correlation, i.e., polyserial correlation between response and total score.
#' \item \code{OR}: Odds ratio for (dichotomized) pairs of items.
#' \item \code{Q3}: Yen's Q3 statistic.
#' \item \code{ISD}: Item-score-distribution, i.e., Pearson's X2 of residuals.
#' \item \code{resp}: Observed responses (integer), predicted responses (integer), and expected responses (numeric).
#' }
#' @inheritParams tidyup_irtree_fit
#' @return Returns an object of class \code{ppc}
#' @seealso \code{\link{print.ppc}} for summarizing the results and
#'   \code{\link{ppc_resp_irtree}} for summarizing posterior predictive response
#'   frequencies.
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
#' # posterior predictive checking
#' res3 <- post_prob_irtree(res2)
#' res4 <- ppc_irtree(prob = res3, fit = res1)
#' res4
#' }
#'
#' @export
ppc_irtree <- function(prob = NULL,
                       statistics = c("item_cor", "OR", "Q3", "ISD", "resp"),
                       fit = NULL,
                       X = NULL,
                       revItem = NULL,
                       traitItem = NULL) {
    
    checkmate::assert_array(prob, mode = "numeric", any.missing = FALSE, d = 4)
    checkmate::qassert(prob, "N+[0,1]")
    statistics <- match.arg(statistics, several.ok = TRUE)
    checkmate::qassert(fit, "l+")
    checkmate::assert_matrix(X, mode = "integerish", any.missing = FALSE,
                             min.rows = 2, min.cols = 2,
                             null.ok = !is.null(fit$args$traitItem))
    checkmate::assert_integerish(X, lower = 1, upper = 5, any.missing = FALSE,
                                 null.ok = !is.null(fit$args$traitItem))
    checkmate::assert_integerish(revItem, lower = 0, upper = 1, any.missing = FALSE,
                                 min.len = 1,
                                 null.ok = !is.null(fit$args$revItem))
    checkmate::assert_integerish(traitItem, lower = 1, any.missing = FALSE,
                                 min.len = 1,
                                 null.ok = !is.null(fit$args$traitItem))
    
    if (is.null(X)) X <- fit$args$X
    if (is.null(revItem)) revItem <- fit$args$revItem
    if (is.null(traitItem)) traitItem <- fit$args$traitItem
    
    R <- dim(prob)[1]
    J <- dim(prob)[3]
    
    message("Drawing samples and computing statistics for ", R,
            " replications.\nThis may take a little time, go get a coffee.")
    X.pred <- posterior_predicted(prob)
    X.exp  <- posterior_expected(prob)
    
    ppc <- list()
    
    # item-total correlation
    if ("item_cor" %in% statistics) {
        time1 <- Sys.time()
        item_cor <- list(obs = statistic_item_cor(X, revItem = revItem,
                                                  traitItem = traitItem),
                         pred = apply(X.pred, 1, statistic_item_cor,
                                      revItem = revItem, traitItem = traitItem))
        item_cor$ppp <- rowMeans(item_cor$obs <= item_cor$pred)
        ppc$item_cor <- item_cor
        time2 <- Sys.time()
        
        message("Finished: Item-total correlation (", round(difftime(time2, time1, units = "mins"), 1), " minutes)")
    }
    
    # global odds ratio (OR)
    if ("OR" %in% statistics) {
        time1 <- Sys.time()
        OR <- list(obs = c(statistic_OR(X)),
                   pred = apply(X.pred, 1, statistic_OR))
        OR$ppp <- rowMeans(OR$obs <= OR$pred)
        ppc$OR <- OR
        time2 <- Sys.time()
        
        message("Finished: Odds ratio (", round(difftime(time2, time1, units = "mins"), 1), " minutes)")
    }
    
    # Yen's Q3
    if ("Q3" %in% statistics) {
        time1 <- Sys.time()
        Q3 <- list(obs = apply(X.exp, 1, function(e) statistic_Q3(X, e)),
                   pred = matrix(NA, choose(J, 2), R))
        for (r in 1:R) {
            Q3$pred[, r] <- statistic_Q3(X.pred[r, , ], X.exp[r, , ])
        }
        Q3$ppp <- rowMeans(Q3$obs <= Q3$pred)
        ppc$Q3 <- Q3
        time2 <- Sys.time()
        
        message("Finished: Yen's Q3 (", round(difftime(time2, time1, units = "mins"), 1), " minutes)")
    }
    
    # # GDDM (Levy & Svetina, 2011)
    # GDDM <- list(obs = apply(X.exp, 1, function(e) statistic_Q3(X, e)),
    #              pred = rep(NA, R))
    # for (r in 1:R) {
    #     GDDM$pred[r] <- statistic_Q3(X.pred[r,,], X.exp[r,,])
    # }
    # GDDM$ppp <- GDDM$obs <= GDDM$pred
    # 
    # message("Finished: GDDM")
    
    # Item score distribution
    if ("ISD" %in% statistics) {
        time1 <- Sys.time()
        ISD <- list(obs = apply(prob, 1, function(e) statistic_item_score(X, e)),
                    # pred = matrix(NA, J, R)
                    pred = sapply(1:R, function(x) statistic_item_score(X.pred[x, , ], prob[x, , , ]))
        )
        # for (r in 1:R){
        #     ISD$pred[, r] <- statistic_item_score(X.pred[r,,], prob[r,,,])
        # }
        ISD$ppp <- rowMeans(ISD$obs <= ISD$pred)
        ppc$ISD <- ISD
        time2 <- Sys.time()
        
        message("Finished: Item score distribution (", round(difftime(time2, time1, units = "mins"), 1), " minutes)")
    }
    
    # ppc <- list(item_cor = item_cor, OR = OR, Q3 = Q3
    #             , ISD = ISD
    #             # , GDDM = GDDM
    #             # , X = list(obs = X, pred = X.pred, exp = X.exp)
    #             )
    if ("resp" %in% statistics) {
        ppc$X <- list(obs = X, pred = X.pred, exp = X.exp)
    }
    class(ppc) <- "ppc"
    return(ppc)
}

#' Print posterior predictive p-values
#'
#' This function calculates and prints PPP-values for objects of class
#' \code{ppc}. More precisely, it prints---across all items---the average of the
#' PPP-values as well as the percentage of extreme PPP-values. For example, for
#' a data set comprised of 20 items, 190 odds ratios can be calculated, namely,
#' one for each pair of items. The present print method calculates the avarage
#' PPP-value (should probably be close to 0.5) as well as the percentage of
#' those 190 PPP-values that are extreme (should hopefully be close to 0).
#' 
#' @param x List as returned from \code{\link{ppc_irtree}}
#' @param na.rm Logical argument passed to \code{\link[base]{mean}}
#'   \code{\link[stats]{median}}. If \code{na.rm = FALSE}, the present print
#'   method will return \code{NA} if at least one PPP-value could not be
#'   calculated for a discrepency measure. Thus this is extremely helpful for
#'   identifying problems; set it to \code{TRUE} only if you know what you are
#'   doing.
#' @param ... Further arguments passed to \code{\link{print.table}}.
#' @return Table of PPP-values
#'
#' @export
print.ppc <- function(x, na.rm = FALSE, ...) {
    # ppp <- sapply(x[1:3], "[[", "ppp")
    ppp <- lapply(x[names(x) != "X"], "[[", "ppp")
    tab <- sapply(ppp, function(p) 
        c(average = mean(p, na.rm = na.rm), 
          median  = median(p, na.rm = na.rm),
          below05 = mean(p <= .05, na.rm = na.rm), 
          above95 = mean(p >= .95, na.rm = na.rm)))
    print.table(tab, ...)
}


#' Posterior predicitve checks of response frequencies.
#' 
#' This function summarizes posterior predicted responses (1, ..., 5) for every
#' item and category. For a single iteration/replication, this results in a
#' response distribution (barplot); across multiple iterations, quantiles (using
#' \code{probs}) are calculated. These quantiles should overlap the observed
#' frequencies.
#' 
#' @param ppc List. Output from \code{\link{ppc_irtree}} that contains
#'   posterior predicted responses with dimensions R x N x J (for N persons, J
#'   items, and R replications/iterations).
#' @param probs Numeric. Vector of probabilities (passed to
#'   \code{\link[stats]{quantile}}) that is used to calculate quantiles of the
#'   posterior predictive distribution. Argument may be named (used to generate
#'   colnames). Defaults to \code{c(.025, .975, .16, .84)}.
#' @return Returns a data frame containing the specified quantiles of the
#'   posterior predictive distribution for every category of every item.
# @inheritParams fit_irtree
# @inheritParams runjags::combine.mcmc
# @importFrom magrittr %>%
# @importFrom rlang .data
# @import coda
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
#' # posterior predictive checking
#' res3 <- post_prob_irtree(res2)
#' res4 <- ppc_irtree(prob = res3, statistics = "resp", fit = res1)
#' res5 <- ppc_resp_irtree(res4)
#' 
#' library(ggplot2)
#' ggplot(res5, aes(x = Categ, y = Obs, ymin = q16, ymax = q84)) +
#'     geom_col() +
#'     geom_errorbar() +
#'     geom_point(aes(y = q50)) + 
#'     facet_wrap(~ Item)
#' }
#' @export
ppc_resp_irtree <- function(ppc = NULL, 
                            probs = NULL) {
    
    if ("ppc" %in% class(ppc)) {
        try(pred <- ppc$X$pred)
        try(obs <- ppc$X$obs)
    } else if (is.list(ppc)) {
        try(pred <- ppc$pred)
        try(obs <- ppc$obs)
    } else {
        pred <- ppc
    }
    
    checkmate::assert_array(pred, mode = "integerish", any.missing = FALSE, d = 3)
    checkmate::qassert(pred, "X+[1,5]")
    
    checkmate::assert_numeric(probs, lower = 0, upper = 1,
                              any.missing = FALSE,  min.len = 1, unique = TRUE,
                              null.ok = TRUE)
    
    if (is.null(probs)) {
        probs <- c(.025, .975, .16, .84, .50)
        names(probs) <- c("q025", "q975", "q16", "q84", "q50")
    } else if (is.null(attr(probs, "names"))) {
        names(probs) <- paste0("q", probs)
    }
    
    a1 <- apply(pred, c(1,3), function(x) prop.table(table(factor(x, levels = 1:5))))
    a2 <- apply(a1, c(1, 3), quantile, probs = probs, names = F)
    names(dimnames(a2)) <- c("Q", "Categ", "Item")
    dimnames(a2)$Q <- names(probs)
    a3 <- reshape2::melt(a2)
    a4 <- reshape2::dcast(a3, value.var = "value",
                          formula = Item + Categ ~ Q)
    a4$Persons <- dim(pred)[2]
    a4$Samples <- dim(pred)[1]
    
    if (exists("obs", inherits = FALSE)) {
        a4 <- apply(obs, 2, function(x) prop.table(table(factor(x, levels = 1:5)))) %>% 
            reshape2::melt(varnames = c("Categ", "Item"), value.name = "Obs") %>% 
            dplyr::full_join(a4, by = c("Categ", "Item"))
    }
    
    return(a4)
}
