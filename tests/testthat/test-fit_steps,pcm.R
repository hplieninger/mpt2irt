library("mpt2irt")
library("magrittr")

generate_pcm <- function(N = NULL, J = NULL, revItem = NULL) {
    if (is.null(revItem)) revItem <- rbinom(J, 1, .33)
    tmp1 <- matrix(rnorm(J*4), J, 4) %>%
        apply(1, sort) %>% 
        rbind(0, .)
    thres <- tmp1 %>% 
        apply(2, cumsum) %>% 
        matrix %>% 
        {do.call(what = cbind, args = rep(list(.), N))}
    tmp2 <- rnorm(N) 
    theta <- tmp2 %>% 
        outer(0:4, .) %>% 
        {do.call(what = rbind, args = rep(list(.), J))}
    
    
    num <- exp(theta - thres)
    num <- array(num, dim = c(5, J, N))
    den <- array(rep(colSums(num), each = 5), dim = c(5, J, N))
    p <- num / den
    
    dat <- 1 + t(apply(p, c(2, 3), function(i) {
        as.integer(findInterval(runif(1), cumsum(i)))
    }))
    
    dat[, revItem == 1] <- 6 - dat[, revItem == 1]
    
    return(list(X = dat,
                theta = tmp2,
                thres = tmp1[-1, ],
                revItem = revItem))
}

# generate_irtree_steps <- function(N = NULL,
#                                   J = NULL,
#                                   revItem = NULL,
#                                   traitItem = NULL) {
#     
#     checkmate::qassert(N, "X1[1,)")
#     checkmate::qassert(J, "X1[1,)")
#     checkmate::assert_integerish(revItem, lower = 0, upper = 1,
#                                  any.missing = FALSE, len = J, null.ok = TRUE)
#     checkmate::assert_integerish(traitItem, lower = 1, upper = 1,
#                                  any.missing = FALSE, len = J, null.ok = TRUE)
#     if (is.null(revItem)) revItem <- rbinom(J, 1, .33)
#     if (is.null(traitItem)) traitItem <- rep(1, J)
#     
#     thres <- matrix(rnorm(J*4), J, 4) %>%
#         apply(1, sort) %>%
#         t
#     theta <- matrix(rnorm(N*max(traitItem)), N, max(traitItem))
#     
#     p_catx <- array(NA, dim = c(N, J, 5))
#     dat <- node1 <- node2 <- node3 <- node4 <- matrix(NA, N, J)
#     
#     for (i in 1:N) {	
#         for (j in 1:J) {
#             node1[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 1]);
#             node2[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 2]);
#             node3[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 3]);
#             node4[i, j] = pnorm(theta[i, traitItem[j]] - thres[j, 4]);
#             
#             p_catx[i,j,1] = (1-node1[i, j]);
#             p_catx[i,j,2] =    node1[i, j] *(1-node2[i, j]);
#             p_catx[i,j,3] =    node1[i, j] *   node2[i, j] *(1-node3[i, j]);
#             p_catx[i,j,4] =    node1[i, j] *   node2[i, j] *   node3[i, j] *(1-node4[i, j]);
#             p_catx[i,j,5] =    node1[i, j] *   node2[i, j] *   node3[i, j] *   node4[i, j] ;
#         }
#     }
#     
#     p_cat <- p_catx
#     
#     p_cat[ , revItem == 1, 5] <- p_catx[ , revItem == 1, 1]
#     p_cat[ , revItem == 1, 4] <- p_catx[ , revItem == 1, 2]
#     p_cat[ , revItem == 1, 3] <- p_catx[ , revItem == 1, 3]
#     p_cat[ , revItem == 1, 2] <- p_catx[ , revItem == 1, 4]
#     p_cat[ , revItem == 1, 1] <- p_catx[ , revItem == 1, 5]
#     
#     dat[, ] <- p_cat %>%
#         apply(1:2, function(x) rmultinom(1, 1, x)) %>%
#         magrittr::equals(1) %>% 
#         apply(2:3, which)
#     
#     return(list(X = dat,
#                 theta = theta,
#                 thres = thres,
#                 revItem = revItem,
#                 traitItem = traitItem))
# }


# DATA GENERATION ---------------------------------------------------------

N <- sample(10:20, 1)
J <- sample(5:20, 1)
# N <- 10
# J <- 5

cond2 <- FALSE
while (cond2 == FALSE) {
    dat1 <- generate_pcm(N = N, J = J)
    cond1 <- suppressWarnings(cor(dat1$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat1$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

cond2 <- FALSE
while (cond2 == FALSE) {
    dat2 <- generate_irtree_steps(N = N, J = J)
    cond1 <- suppressWarnings(cor(dat2$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat2$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}
suppressWarnings(rm(cond1, cond2))

test_that("generate_pcm() returns correct output", {
    expect_is(dat1, "list")
    expect_is(dat2, "list")
    expect_equal(ncol(dat1$X), J)
    expect_equal(ncol(dat2$X), J)
    expect_equal(nrow(dat1$X), N)
    expect_equal(nrow(dat2$X), N)
})

# MODEL FITTING -----------------------------------------------------------

M <- 200
warmup <- 100

invisible(capture.output(
    # res1 <- fit_irtree(dat1$X, revItem = dat1$revItem,
    #                    M = M, warmup = warmup, n.chains = 1,
    #                    fitModel = "pcm", fitMethod = "jags"),
    res2 <- fit_irtree(dat1$X, fitModel = "pcm", fitMethod = "stan",
                       revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 1),
    # res3 <- fit_irtree(dat2$X, revItem = dat2$revItem,
    #                    M = M, warmup = warmup, n.chains = 1,
    #                    fitModel = "steps", fitMethod = "jags"),
    res4 <- fit_irtree(dat2$X, fitModel = "steps", fitMethod = "stan",
                       revItem = dat2$revItem,
                       M = M, warmup = warmup, n.chains = 1)
))

test_that("fit_irtree() returns MCMC list", {
    # expect_equal(unique(sapply(res1$samples$mcmc, nrow)), M)
    expect_error(fit_irtree(dat1$X, revItem = dat1$revItem,
                            M = M, warmup = warmup, n.chains = 1,
                            fitModel = "pcm", fitMethod = "jags"))
    expect_equal(unique(sapply(rstan::As.mcmc.list(res2$samples), nrow)), M)
    # expect_equal(unique(sapply(res3$samples$mcmc, nrow)), M)
    expect_error(fit_irtree(dat2$X, revItem = dat2$revItem,
                            M = M, warmup = warmup, n.chains = 1,
                            fitModel = "steps", fitMethod = "jags"))
    expect_equal(unique(sapply(rstan::As.mcmc.list(res4$samples), nrow)), M)
})


# SUMMARIZING MODEL RESULTS -----------------------------------------------

# res1b <- summarize_irtree_fit(res1)
# res1c <- tidyup_irtree_fit(res1b)
# res1d <- suppressMessages(pp_irtree(res1b, iter = iter, N = N))

res2b <- summarize_irtree_fit(res2)
res2c <- tidyup_irtree_fit(res2b)

# res3b <- summarize_irtree_fit(res3)
# res3c <- tidyup_irtree_fit(res3b)
# res3d <- suppressMessages(pp_irtree(res3b, iter = iter, N = N))

res4b <- summarize_irtree_fit(res4)
res4c <- tidyup_irtree_fit(res4b)

test_that("tidyup_irtree_fit() returns correlations", {
    # expect_equal(unique(as.vector(sapply(res1c$Corr, dim))), 1)
    expect_equal(unique(as.vector(sapply(res2c$Corr, dim))), 1)
    # expect_equal(unique(as.vector(sapply(res3c$Corr, dim))), 1)
    expect_equal(unique(as.vector(sapply(res4c$Corr, dim))), 1)
    # expect_equal(unique(as.vector(sapply(res1c$Sigma, dim))), 1)
    expect_equal(unique(as.vector(sapply(res2c$Sigma, dim))), 1)
    # expect_equal(unique(as.vector(sapply(res3c$Sigma, dim))), 1)
    expect_equal(unique(as.vector(sapply(res4c$Sigma, dim))), 1)
})

test_that("tidyup_irtree_fit() returns correct number of parameters", {
    # expect_equal(unique(sapply(res1c$beta, nrow)), J)
    expect_equal(unique(sapply(res2c$beta, nrow)), J)
    # expect_equal(unique(sapply(res3c$beta, nrow)), J)
    expect_equal(unique(sapply(res4c$beta, nrow)), J)
    # expect_equal(unique(sapply(res1c$theta, nrow)), N)
    expect_equal(unique(sapply(res2c$theta, nrow)), N)
    # expect_equal(unique(sapply(res3c$theta, nrow)), N)
    expect_equal(unique(sapply(res4c$theta, nrow)), N)
})

test_that("plot_irtree() returns a valid ggplot", {
    expect_true(ggplot2::is.ggplot(res2c$plot))
    expect_true(ggplot2::is.ggplot(res4c$plot))
})

# RECOVERY ----------------------------------------------------------------

# test_that("Check that true model parameters are correctly recovered", {
#     cor1 <- cor(dat1$theta, res2c$theta$Median)
#     # expect_gt(cor1, .7)
#     cor2 <- cor(dat2$theta, res4c$theta$Median)
#     expect_gt(cor2, .8)
#     # tryCatch(expect_gt(cor2, .8),
#     #          expectation_failure = function(x) {
#     #              message(x)
#     #              message(sprintf(c("Model 'steps' -- theta -- r=%.2f, N=%i, J=%i"),
#     #                              cor2, N, J))
#     #          })
#     cor3 <- cor(as.vector(dat1$thres), as.vector(res2c$beta$Median))
#     # expect_gt(cor3, .6)
#     cor4 <- cor(as.vector(dat2$thres), as.vector(res4c$beta$Median))
#     expect_gt(cor3, .65)
# })

# PPC ---------------------------------------------------------------------

res2d <- post_prob_irtree(res2b, iter = 20)
res2e <- ppc_irtree(prob = res2d, fit = res2b)
invisible(capture.output(res2f <- print(res2e, na.rm = TRUE)))
res2g <- ppc_resp_irtree(res2e)

res4d <- post_prob_irtree(res4b, iter = 20)
res4e <- ppc_irtree(prob = res4d, fit = res4b)
invisible(capture.output(res4f <- print(res4e, na.rm = TRUE)))
res4g <- ppc_resp_irtree(res4e)

test_that("ppc_resp_irtree() returns valid values", {
    expect_is(res2f, "matrix")
    expect_is(res4f, "matrix")
    
    expect_gte(min(res2f), 0)
    expect_gte(min(res4f), 0)
    
    expect_lte(min(res2f), 1)
    expect_lte(min(res4f), 1)
})

test_that("print(ppc_irtree()) returns valid values", {
    expect_is(res2g, "data.frame")
    expect_is(res4g, "data.frame")
    
    expect_equal(as.numeric(unique(res2g$Item)), 1:J)
    expect_equal(as.numeric(unique(res4g$Item)), 1:J)
    
    expect_equal(as.numeric(unique(res2g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res4g$Categ)), 1:5)
    
    expect_equal(unique(res2g$Persons), N)
    expect_equal(unique(res4g$Persons), N)
    
    expect_gte(min(subset(res2g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res4g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    
    expect_lte(min(subset(res2g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res4g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
})
