# detach(package:magrittr)
library("mpt2irt")
library("magrittr")

# DATA GENERATION ---------------------------------------------------------

context("ext2: Data generation")

N <- sample(10:20, 1)
J <- sample(2:20, 1)
# N <- 10
# J <- 2

betas1 <- mpt2irt:::gen_betas(genModel = "ext2", J = J)

cond2 <- FALSE
while (cond2 == FALSE) {
    dat1 <- generate_irtree_ext(N = N,
                                J = J,
                                betas = betas1,
                                prop.rev = runif(1),
                                genModel = "ext2")
    cond1 <- suppressWarnings(cor(dat1$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat1$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

N2 <- sample(10:20, 1)
J2 <- sample(6:20, 1)

betas2 <- mpt2irt:::gen_betas(genModel = "ext2", J = J2)

# Data 'ext'
cond2 <- FALSE
while (cond2 == FALSE) {
    dat2 <- generate_irtree_ext(N = N2,
                                J = J2,
                                betas = betas2,
                                traitItem = c(sample(rep(1:2, 3)), sample(1:2, J2 - 6, T)), 
                                prop.rev = runif(1),
                                genModel = "ext2")
    cond1 <- suppressWarnings(cor(dat2$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat2$X[, dat2$revItem == 1] %>% cor %>% sign %>% magrittr::equals(-1) %>% any
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
    expect_equal(ncol(dat2$X), J2)
    expect_equal(nrow(dat1$X), N)
    expect_equal(nrow(dat2$X), N2)
})

# MODEL FITTING -----------------------------------------------------------

context("ext2: Model fitting")

M <- 200
warmup <- 100

res1 <- fit_irtree(dat1$X, fitModel = "ext2", fitMethod = "stan",
                   revItem = dat1$revItem, traitItem = dat1$traitItem,
                   M = M, warmup = warmup, n.chains = 1)
res2 <- fit_irtree(dat2$X, fitModel = "ext2", fitMethod = "stan",
                   revItem = dat2$revItem, traitItem = dat2$traitItem,
                   M = M, warmup = warmup, n.chains = 1)

test_that("fit_irtree() returns MCMC list", {
    expect_equal(unique(sapply(rstan::As.mcmc.list(res1$samples), nrow)), M)
    expect_equal(unique(sapply(rstan::As.mcmc.list(res2$samples), nrow)), M)
})

# SUMMARIZING MODEL RESULTS -----------------------------------------------

context("ext2: Summarizing fitted models")

iter <- 10

res1b <- summarize_irtree_fit(res1)
res1c <- tidyup_irtree_fit(res1b)
res1d <- suppressMessages(pp_irtree(res1b, iter = iter, N = N))

res2b <- summarize_irtree_fit(res2)
res2c <- tidyup_irtree_fit(res2b)
res2d <- suppressMessages(pp_irtree(res2b, iter = iter, N = N2))

test_that("tidyup_irtree_fit() returns correlations", {
    expect_equal(unique(as.vector(sapply(res1c$Corr, dim))), res1b$args$S)
    expect_equal(unique(as.vector(sapply(res1c$Sigma, dim))), res1b$args$S)
    expect_equal(unique(as.vector(sapply(res2c$Corr, dim))), res2b$args$S)
    expect_equal(unique(as.vector(sapply(res2c$Sigma, dim))), res2b$args$S)
})

test_that("tidyup_irtree_fit() returns correct number of parameters", {
    expect_equal(unique(sapply(res1c$beta, nrow)), J)
    expect_equal(unique(sapply(res1c$theta, nrow)), N)
    expect_equal(unique(sapply(res2c$beta, nrow)), J2)
    expect_equal(unique(sapply(res2c$theta, nrow)), N2)
})

test_that("plot_irtree() returns a valid ggplot", {
    expect_true(ggplot2::is.ggplot(res1c$plot))
    expect_true(ggplot2::is.ggplot(res2c$plot))
})

test_that("pp_irtree() returns valid values", {
    expect_is(res1d, "data.frame")
    expect_is(res2d, "data.frame")
    
    expect_equal(as.numeric(levels(res1d$Item)), 1:J)
    expect_equal(as.numeric(levels(res2d$Item)), 1:J2)
    
    expect_equal(as.numeric(levels(res1d$Categ)), 1:5)
    expect_equal(as.numeric(levels(res2d$Categ)), 1:5)
    
    expect_equal(unique(res1d$Persons), N)
    expect_equal(unique(res2d$Persons), N2)
    
    expect_equal(unique(res1d$Samples), iter*res1b$args$n.chains)
    expect_equal(unique(res2d$Samples), iter*res2b$args$n.chains)
    
    expect_gte(min(subset(res1d, select = -(Item:Samples))), 0)
    expect_gte(min(subset(res2d, select = -(Item:Samples))), 0)
    
    expect_lte(min(subset(res1d, select = -(Item:Samples))), 1)
    expect_lte(min(subset(res2d, select = -(Item:Samples))), 1)
})

# RECOVERY ----------------------------------------------------------------

# context("ext2: Recovery")
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
