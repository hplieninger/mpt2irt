context("Multidimensional Models 'ext', '2012', 'shift', and 'steps'")

# detach(package:magrittr)
library("mpt2irt")
library("magrittr")

# DATA GENERATION ---------------------------------------------------------

context("Multidim: Data generation")

N <- sample(10:20, 1)
J <- sample(6:20, 1)
# N <- 10
# J <- 6

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

context("Multidim: Model fitting")

M <- 200
warmup <- 100

invisible(capture.output(
    res1 <- fit_irtree(dat1$X, fitModel = "2012", fitMethod = "jags",
                       revItem = dat1$revItem, traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1),
    res2 <- fit_irtree(dat1$X, fitModel = "ext", fitMethod = "stan",
                       revItem = dat1$revItem, traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1),
    res3 <- fit_irtree(dat2$X, fitModel = "ext", fitMethod = "jags",
                       revItem = dat2$revItem, traitItem = dat2$traitItem,
                       M = M, warmup = warmup, n.chains = 1),
    res4 <- fit_irtree(dat2$X, fitModel = "2012", fitMethod = "stan",
                       revItem = dat2$revItem, traitItem = dat2$traitItem,
                       M = M, warmup = warmup, n.chains = 1),
    res5 <- fit_irtree(dat3$X, fitModel = "steps", fitMethod = "stan",
                       revItem = dat3$revItem, traitItem = dat3$traitItem,
                       M = M, warmup = warmup, n.chains = 1),
    res6 <- fit_irtree(dat1$X, fitModel = "shift", fitMethod = "stan",
                       revItem = dat1$revItem, traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1)
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

context("Multidim: Summarizing fitted models")

res1b <- summarize_irtree_fit(res1)
res1c <- tidyup_irtree_fit(res1b)

res2b <- summarize_irtree_fit(res2)
res2c <- tidyup_irtree_fit(res2b)

res3b <- summarize_irtree_fit(res3)
res3c <- tidyup_irtree_fit(res3b)

res4b <- summarize_irtree_fit(res4)
res4c <- tidyup_irtree_fit(res4b)

res5b <- summarize_irtree_fit(res5)
res5c <- tidyup_irtree_fit(res5b)

res6b <- summarize_irtree_fit(res6)
res6c <- tidyup_irtree_fit(res6b)

test_that("tidyup_irtree_fit() returns correlations", {
    expect_equal(unique(as.vector(sapply(res1c$Corr, dim))), res1b$args$S)
    expect_equal(unique(as.vector(sapply(res2c$Corr, dim))), res2b$args$S)
    expect_equal(unique(as.vector(sapply(res3c$Corr, dim))), res3b$args$S)
    expect_equal(unique(as.vector(sapply(res4c$Corr, dim))), res4b$args$S)
    expect_equal(unique(as.vector(sapply(res5c$Corr, dim))), res5b$args$S)
    expect_equal(unique(as.vector(sapply(res6c$Corr, dim))), res6b$args$S)
    expect_equal(unique(as.vector(sapply(res1c$Sigma, dim))), res1b$args$S)
    expect_equal(unique(as.vector(sapply(res2c$Sigma, dim))), res2b$args$S)
    expect_equal(unique(as.vector(sapply(res3c$Sigma, dim))), res3b$args$S)
    expect_equal(unique(as.vector(sapply(res4c$Sigma, dim))), res4b$args$S)
    expect_equal(unique(as.vector(sapply(res5c$Sigma, dim))), res5b$args$S)
    expect_equal(unique(as.vector(sapply(res6c$Sigma, dim))), res6b$args$S)
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

test_that("plot_irtree() returns a valid ggplot", {
    expect_true(ggplot2::is.ggplot(res1c$plot))
    expect_true(ggplot2::is.ggplot(res2c$plot))
    expect_true(ggplot2::is.ggplot(res3c$plot))
    expect_true(ggplot2::is.ggplot(res4c$plot))
    expect_true(ggplot2::is.ggplot(res5c$plot))
    expect_true(ggplot2::is.ggplot(res6c$plot))
})

# PPC ---------------------------------------------------------------------

context("PPC")

res1d <- post_prob_irtree(res1b, iter = 20)
res1e <- ppc_irtree(prob = res1d, fit = res1b)
invisible(capture.output(res1f <- print(res1e, na.rm = TRUE)))
res1g <- ppc_resp_irtree(res1e)

res2d <- post_prob_irtree(res2b, iter = 20)
res2e <- ppc_irtree(prob = res2d, fit = res2b)
invisible(capture.output(res2f <- print(res2e, na.rm = TRUE)))
res2g <- ppc_resp_irtree(res2e)

res3d <- post_prob_irtree(res3b, iter = 20)
res3e <- ppc_irtree(prob = res3d, fit = res3b)
invisible(capture.output(res3f <- print(res3e, na.rm = TRUE)))
res3g <- ppc_resp_irtree(res3e)

res4d <- post_prob_irtree(res4b, iter = 20)
res4e <- ppc_irtree(prob = res4d, fit = res4b)
invisible(capture.output(res4f <- print(res4e, na.rm = TRUE)))
res4g <- ppc_resp_irtree(res4e)

res5d <- post_prob_irtree(res5b, iter = 20)
res5e <- ppc_irtree(prob = res5d, fit = res5b)
invisible(capture.output(res5f <- print(res5e, na.rm = TRUE)))
res5g <- ppc_resp_irtree(res5e)

res6d <- post_prob_irtree(res6b, iter = 20)
res6e <- ppc_irtree(prob = res6d, fit = res6b)
invisible(capture.output(res6f <- print(res6e, na.rm = TRUE)))
res6g <- ppc_resp_irtree(res6e)

test_that("ppc_resp_irtree() returns valid values", {
    expect_is(res1f, "matrix")
    expect_is(res2f, "matrix")
    expect_is(res3f, "matrix")
    expect_is(res4f, "matrix")
    expect_is(res5f, "matrix")
    expect_is(res6f, "matrix")
    
    expect_gte(min(res1f), 0)
    expect_gte(min(res2f), 0)
    expect_gte(min(res3f), 0)
    expect_gte(min(res4f), 0)
    expect_gte(min(res5f), 0)
    expect_gte(min(res6f), 0)
    
    expect_lte(min(res1f), 1)
    expect_lte(min(res2f), 1)
    expect_lte(min(res3f), 1)
    expect_lte(min(res4f), 1)
    expect_lte(min(res5f), 1)
    expect_lte(min(res6f), 1)
})

test_that("print(ppc_irtree()) returns valid values", {
    expect_is(res1g, "data.frame")
    expect_is(res2g, "data.frame")
    expect_is(res3g, "data.frame")
    expect_is(res4g, "data.frame")
    expect_is(res5g, "data.frame")
    expect_is(res6g, "data.frame")
    
    expect_equal(as.numeric(unique(res1g$Item)), 1:J)
    expect_equal(as.numeric(unique(res2g$Item)), 1:J)
    expect_equal(as.numeric(unique(res3g$Item)), 1:J)
    expect_equal(as.numeric(unique(res4g$Item)), 1:J)
    expect_equal(as.numeric(unique(res5g$Item)), 1:J)
    expect_equal(as.numeric(unique(res6g$Item)), 1:J)
    
    expect_equal(as.numeric(unique(res1g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res2g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res3g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res4g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res5g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res6g$Categ)), 1:5)
    
    expect_equal(unique(res1g$Persons), N)
    expect_equal(unique(res2g$Persons), N)
    expect_equal(unique(res3g$Persons), N)
    expect_equal(unique(res4g$Persons), N)
    expect_equal(unique(res5g$Persons), N)
    expect_equal(unique(res6g$Persons), N)
    
    expect_gte(min(subset(res1g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res2g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res3g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res4g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res5g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res6g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    
    expect_lte(min(subset(res1g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res2g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res3g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res4g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res5g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res6g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
})
