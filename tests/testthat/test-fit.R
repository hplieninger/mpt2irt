# detach(package:magrittr)
library("mpt2irt")
library("magrittr")

# DATA GENERATION ---------------------------------------------------------

context("Data generation")

N <- sample(10:20, 1)
J <- sample(2:20, 1)
# N <- 10
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
    cond1 <- suppressWarnings(cor(dat1$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
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
    cond1 <- suppressWarnings(cor(dat2$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat2$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}
suppressWarnings(rm(cond1, cond2))

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

M <- 200
warmup <- 200

invisible(capture.output(
    res1 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "2012", fitMethod = "jags"),
    res2 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       M = M, warmup = warmup, n.chains = 2,
                       fitModel = "ext", fitMethod = "stan"),
    res3 <- fit_irtree(dat2$X, revItem = dat2$revItem,
                       M = M, warmup = warmup, n.chains = 2,
                       fitModel = "ext", fitMethod = "jags"),
    res4 <- fit_irtree(dat2$X, revItem = dat2$revItem,
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

iter <- 10

res1b <- summarize_irtree_fit(res1)
res1c <- tidyup_irtree_fit(res1b)
res1d <- suppressMessages(pp_irtree(res1b, iter = iter, N = N))

res2b <- summarize_irtree_fit(res2)
res2c <- tidyup_irtree_fit(res2b)
res2d <- suppressMessages(pp_irtree(res2b, iter = iter, N = N))

res3b <- summarize_irtree_fit(res3)
res3c <- tidyup_irtree_fit(res3b)
res3d <- suppressMessages(pp_irtree(res3b, iter = iter, N = N))

res4b <- summarize_irtree_fit(res4)
res4c <- tidyup_irtree_fit(res4b)
res4d <- suppressMessages(pp_irtree(res4b, iter = iter, N = N))

test_that("tidyup_irtree_fit() returns correlations", {
    expect_equal(unique(as.vector(sapply(res1c$Corr, dim))), 3)
    expect_equal(unique(as.vector(sapply(res2c$Corr, dim))), 4)
    expect_equal(unique(as.vector(sapply(res3c$Corr, dim))), 4)
    expect_equal(unique(as.vector(sapply(res4c$Corr, dim))), 3)
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

test_that("pp_irtree() returns valid values", {
    expect_is(res1d, "data.frame")
    expect_is(res2d, "data.frame")
    expect_is(res3d, "data.frame")
    expect_is(res4d, "data.frame")
    
    expect_equal(as.numeric(levels(res1d$Item)), 1:J)
    expect_equal(as.numeric(levels(res2d$Item)), 1:J)
    expect_equal(as.numeric(levels(res3d$Item)), 1:J)
    expect_equal(as.numeric(levels(res4d$Item)), 1:J)
    
    expect_equal(as.numeric(levels(res1d$Categ)), 1:5)
    expect_equal(as.numeric(levels(res2d$Categ)), 1:5)
    expect_equal(as.numeric(levels(res3d$Categ)), 1:5)
    expect_equal(as.numeric(levels(res4d$Categ)), 1:5)
    
    # expect_equal(ncol(res1d), 9)
    # expect_equal(ncol(res2d), 9)
    # expect_equal(ncol(res3d), 9)
    # expect_equal(ncol(res4d), 9)
    
    expect_equal(unique(res1d$Persons), N)
    expect_equal(unique(res2d$Persons), N)
    expect_equal(unique(res3d$Persons), N)
    expect_equal(unique(res4d$Persons), N)
    
    expect_equal(unique(res1d$Samples), iter*res1b$args$n.chains)
    expect_equal(unique(res2d$Samples), iter*res2b$args$n.chains)
    expect_equal(unique(res3d$Samples), iter*res3b$args$n.chains)
    expect_equal(unique(res4d$Samples), iter*res4b$args$n.chains)
    
    expect_gte(min(subset(res1d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res2d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res3d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res4d, select = -(Item:Samples))), 0)
    
    expect_lte(min(subset(res1d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res2d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res3d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res4d, select = -(Item:Samples))), 1)
})

# RECOVERY ----------------------------------------------------------------

# context("Recovery")
# 
# test_that("Check that true model parameters are correctly recovered", {
#     cor1 <- cor(as.vector(dat1$theta[, c(1, 2, 4)]), as.vector(res1c$theta$Median))
#     # expect_gt(cor1, .7)
#     cor2 <- cor(as.vector(dat1$theta), as.vector(res2c$theta$Median))
#     # expect_gt(cor2, .8)
#     cor3 <- cor(as.vector(dat2$theta), as.vector(res3c$theta$Median[, c(1, 2, 4)]))
#     # expect_gt(cor3, .6)
#     cor4 <- cor(as.vector(dat2$theta), as.vector(res4c$theta$Median))
#     # expect_gt(cor3, .65)
# })
