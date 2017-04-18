library("mpt2irt")
library("magrittr")

# DATA GENERATION ---------------------------------------------------------

context("Data generation")

N <- sample(2:20, 1)
J <- sample(2:20, 1)
# N <- 2
# J <- 2

betas1 <- mpt2irt:::gen_betas(genModel = "ext", J = J)
betas2 <- mpt2irt:::gen_betas(genModel = "2012", J = J)

# dat1 <- generate_irtree_ext(N = N,
#                             J = J,
#                             betas = betas1,
#                             prop.rev = runif(1),
#                             genModel = "ext",
#                             beta_ARS_extreme = rnorm(1))
# dat2 <- generate_irtree_2012(N = N,
#                              J = J,
#                              betas = betas2,
#                              prop.rev = runif(1))

cond2 <- FALSE
while (cond2 == FALSE) {
    dat1 <- generate_irtree_ext(N = N,
                                J = J,
                                betas = betas1,
                                prop.rev = runif(1),
                                genModel = "ext",
                                beta_ARS_extreme = rnorm(1))
    cond1 <- dat1$X %>% cor %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat1$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

cond2 <- FALSE
while (cond2 == FALSE) {
    dat2 <- generate_irtree_2012(N = N,
                                 J = J,
                                 betas = betas2,
                                 prop.rev = runif(1))
    cond1 <- dat2$X %>% cor %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat2$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}
rm(cond1, cond2)

test_that("gen_betas() returns matrix", {
    expect_is(betas1, "matrix")
    expect_is(betas2, "matrix")
})

test_that("generate_irtree() returns correct output", {
    expect_is(dat1, "list")
    expect_is(dat2, "list")
    expect_equal(ncol(dat1$X), J)
    expect_equal(ncol(dat2$X), J)
    expect_equal(nrow(dat1$X), N)
    expect_equal(nrow(dat2$X), N)
})

# MODEL FITTING -----------------------------------------------------------

context("Model fitting")

M <- 100
warmup <- 100

invisible(capture.output(
    res1 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "2012", fitMethod = "jags"),
    res2 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "ext", fitMethod = "stan"),
    res3 <- fit_irtree(dat2$X, revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "ext", fitMethod = "jags"),
    res4 <- fit_irtree(dat2$X, revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "2012", fitMethod = "stan")
))

test_that("fit_irtree() returns MCMC list", {
    expect_equal(unique(sapply(res1$samples$mcmc, nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res2$samples), nrow)), M)
    expect_equal(unique(sapply(res3$samples$mcmc, nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res4$samples), nrow)), M)
})

# res1 <- fit_irtree(dat1$X, revItem = dat1$revItem,
#                    M = 100, warmup = 100, n.chains = 1,
#                    fitModel = "2012", fitMethod = "jags")
# res2 <- fit_irtree(dat1$X, revItem = dat1$revItem,
#                    M = 100, warmup = 100, n.chains = 1,
#                    fitModel = "ext", fitMethod = "stan")
# res3 <- fit_irtree(dat2$X, revItem = dat1$revItem,
#                    M = 100, warmup = 100, n.chains = 1,
#                    fitModel = "ext", fitMethod = "jags")
# res4 <- fit_irtree(dat2$X, revItem = dat1$revItem,
#                    M = 100, warmup = 100, n.chains = 1,
#                    fitModel = "2012", fitMethod = "stan")


# SUMMARIZING MODEL RESULTS -----------------------------------------------

context("Summarizing fitted models")

res1b <- summarize_irtree_fit(res1)
res1c <- tidyup_irtree_fit(res1b, N = N, J = J, revItem = dat1$revItem,
                           traitItem = dat1$traitItem, fitModel = res1$fitModel,
                           fitMethod = res1$fitMethod)
res2b <- summarize_irtree_fit(res2)
res2c <- tidyup_irtree_fit(res2b, N = N, J = J, revItem = dat1$revItem,
                           traitItem = dat1$traitItem, fitModel = res2$fitModel,
                           fitMethod = res2$fitMethod)
res3b <- summarize_irtree_fit(res3)
res3c <- tidyup_irtree_fit(res3b, N = N, J = J, revItem = dat2$revItem,
                           traitItem = dat2$traitItem, fitModel = res3$fitModel,
                           fitMethod = res3$fitMethod)
res4b <- summarize_irtree_fit(res4)
res4c <- tidyup_irtree_fit(res4b, N = N, J = J, revItem = dat2$revItem,
                           traitItem = dat2$traitItem, fitModel = res4$fitModel,
                           fitMethod = res4$fitMethod)

test_that("tidyup_irtree_fit() returns correlations", {
    expect_equal(length(res1c$Corr), 4)
    expect_equal(length(res2c$Corr), 4)
    expect_equal(length(res3c$Corr), 4)
    expect_equal(length(res4c$Corr), 4)
    expect_equal(unique(as.vector(sapply(res1c$Sigma, dim))), 3)
    expect_equal(unique(as.vector(sapply(res2c$Sigma, dim))), 4)
    expect_equal(unique(as.vector(sapply(res3c$Sigma, dim))), 4)
    expect_equal(unique(as.vector(sapply(res4c$Sigma, dim))), 3)
})

test_that("tidyup_irtree_fit() returns correct number of parameters", {
    expect_equal(unique(sapply(res1c$beta, nrow)), J)
    expect_equal(unique(sapply(res2c$beta, nrow)), J)
    expect_equal(unique(sapply(res3c$beta, nrow)), J)
    expect_equal(unique(sapply(res4c$beta, nrow)), J)
    expect_equal(unique(sapply(res1c$theta, nrow)), N)
    expect_equal(unique(sapply(res2c$theta, nrow)), N)
    expect_equal(unique(sapply(res3c$theta, nrow)), N)
    expect_equal(unique(sapply(res4c$theta, nrow)), N)
})


# # posterior predictives for 512 hypothetical persons
# res2d <- pp_irtree(res2b$mcmc, iter = 10, N = 4, revItem = dat1$revItem,
#                    traitItem = dat1$traitItem, fitModel = res2$fitModel)
# 
# plot_irtree(res2b, J = length(dat1$traitItem), revItem = dat1$revItem)
