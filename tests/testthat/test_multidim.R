# detach(package:magrittr)
library("mpt2irt")
library("magrittr")

# DATA GENERATION ---------------------------------------------------------

context("Data generation")

N <- sample(10:20, 1)
J <- sample(6:20, 1)

betas1 <- mpt2irt:::gen_betas(genModel = "ext", J = J)
betas2 <- mpt2irt:::gen_betas(genModel = "2012", J = J)

test_that("gen_betas() returns matrix", {
    expect_is(betas1, "matrix")
    expect_is(betas2, "matrix")
})

# Data 'ext'
cond2 <- FALSE
while (cond2 == FALSE) {
    dat1 <- generate_irtree_ext(N = N,
                                J = J,
                                betas = betas1,
                                traitItem = c(sample(rep(1:2, 3)), sample(1:2, J - 6, T)), 
                                prop.rev = runif(1),
                                genModel = "ext",
                                beta_ARS_extreme = rnorm(1))
    cond1 <- suppressWarnings(cor(dat1$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat1$X[, dat1$revItem == 1] %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

# Data 2012
cond2 <- FALSE
while (cond2 == FALSE) {
    dat2 <- generate_irtree_2012(N = N,
                                 J = J,
                                 betas = betas2,
                                 traitItem = c(sample(rep(1:2, 3)), sample(1:2, J - 6, T)),
                                 prop.rev = runif(1))
    cond1 <- suppressWarnings(cor(dat2$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat2$X[, dat2$revItem == 1] %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

# Data Steps
cond2 <- FALSE
while (cond2 == FALSE) {
    dat3 <- generate_irtree_steps(N = N, J = J,
                                  traitItem = c(sample(rep(1:2, 3)), sample(1:2, J - 6, T)))
    cond1 <- suppressWarnings(cor(dat3$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat3$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

suppressWarnings(rm(cond1, cond2))


test_that("generate_irtree() returns correct output", {
    expect_is(dat1, "list")
    expect_is(dat2, "list")
    expect_is(dat3, "list")
    expect_equal(ncol(dat1$X), J)
    expect_equal(ncol(dat2$X), J)
    expect_equal(ncol(dat3$X), J)
    expect_equal(nrow(dat1$X), N)
    expect_equal(nrow(dat2$X), N)
    expect_equal(nrow(dat3$X), N)
})

# MODEL FITTING -----------------------------------------------------------

context("Model fitting")

M <- 200
warmup <- 100

invisible(capture.output(
    res1 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "2012", fitMethod = "jags"),
    res2 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "ext", fitMethod = "stan"),
    res3 <- fit_irtree(dat2$X, revItem = dat2$revItem,
                       traitItem = dat2$traitItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "ext", fitMethod = "jags"),
    res4 <- fit_irtree(dat2$X, revItem = dat2$revItem,
                       traitItem = dat2$traitItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "2012", fitMethod = "stan"),
    res5 <- fit_irtree(dat3$X, revItem = dat3$revItem,
                       traitItem = dat3$traitItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "steps", fitMethod = "stan"),
    res6 <- fit_irtree(dat1$X, revItem = dat1$revItem,
                       traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1,
                       fitModel = "shift", fitMethod = "stan")
))

dat_over <- data.frame(name = paste0("res", 1:6),
                       dat = c(2012, 2012, "ext", "ext", "shift", 2012),
                       fitModel = NA, fitMethod = NA)

for (iii in 1:6) {
    assign("tmp1", get(paste0("res", iii)))
    dat_over$fitModel[iii] <- tmp1$args$fitModel
    dat_over$fitMethod[iii] <- tmp1$args$fitMethod
}

test_that("fit_irtree() returns MCMC list", {
    expect_equal(unique(sapply(res1$samples$mcmc, nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res2$samples), nrow)), M)
    expect_equal(unique(sapply(res3$samples$mcmc, nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res4$samples), nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res5$samples), nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res6$samples), nrow)), M)
})

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

res5b <- summarize_irtree_fit(res5)
res5c <- tidyup_irtree_fit(res5b)
res5d <- suppressMessages(pp_irtree(res5b, iter = iter, N = N))

res6b <- summarize_irtree_fit(res6)
res6c <- tidyup_irtree_fit(res6b)
# res6d <- suppressMessages(pp_irtree(res6b, iter = iter, N = N))
expect_error(pp_irtree(res6b, iter = iter, N = N))

test_that("tidyup_irtree_fit() returns correlations", {
    expect_equal(unique(as.vector(sapply(res1c$Corr, dim))), 3)
    expect_equal(unique(as.vector(sapply(res2c$Corr, dim))), 4)
    expect_equal(unique(as.vector(sapply(res3c$Corr, dim))), 4)
    expect_equal(unique(as.vector(sapply(res4c$Corr, dim))), 3)
    expect_equal(unique(as.vector(sapply(res5c$Corr, dim))), 2)
    expect_equal(unique(as.vector(sapply(res6c$Corr, dim))), 5)
    expect_equal(unique(as.vector(sapply(res1c$Sigma, dim))), 3)
    expect_equal(unique(as.vector(sapply(res2c$Sigma, dim))), 4)
    expect_equal(unique(as.vector(sapply(res3c$Sigma, dim))), 4)
    expect_equal(unique(as.vector(sapply(res4c$Sigma, dim))), 3)
    expect_equal(unique(as.vector(sapply(res5c$Sigma, dim))), 2)
    expect_equal(unique(as.vector(sapply(res6c$Sigma, dim))), 5)
})

test_that("tidyup_irtree_fit() returns correct number of parameters", {
    expect_equal(unique(sapply(res1c$beta, nrow)), J)
    expect_equal(unique(sapply(res2c$beta, nrow)), J)
    expect_equal(unique(sapply(res3c$beta, nrow)), J)
    expect_equal(unique(sapply(res4c$beta, nrow)), J)
    expect_equal(unique(sapply(res5c$beta, nrow)), J)
    expect_equal(unique(sapply(res6c$beta, nrow)), J)
    expect_equal(unique(sapply(res1c$theta, nrow)), N)
    expect_equal(unique(sapply(res2c$theta, nrow)), N)
    expect_equal(unique(sapply(res3c$theta, nrow)), N)
    expect_equal(unique(sapply(res4c$theta, nrow)), N)
    expect_equal(unique(sapply(res5c$theta, nrow)), N)
    expect_equal(unique(sapply(res6c$theta, nrow)), N)
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
    
    expect_equal(ncol(res1d), 8)
    expect_equal(ncol(res2d), 8)
    expect_equal(ncol(res3d), 8)
    expect_equal(ncol(res4d), 8)
    
    expect_equal(unique(res1d$Persons), N)
    expect_equal(unique(res2d$Persons), N)
    expect_equal(unique(res3d$Persons), N)
    expect_equal(unique(res4d$Persons), N)
    
    expect_equal(unique(res1d$Samples), iter)
    expect_equal(unique(res2d$Samples), iter)
    expect_equal(unique(res3d$Samples), iter)
    expect_equal(unique(res4d$Samples), iter)
    
    expect_gte(min(subset(res1d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res2d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res3d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res4d, select = -(Item:Samples))), 0)
    
    expect_lte(min(subset(res1d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res2d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res3d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res4d, select = -(Item:Samples))), 1)
})
