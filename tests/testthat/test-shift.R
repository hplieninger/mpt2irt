context("Models: 'shift'")

# detach(package:magrittr)
library("mpt2irt")
library("magrittr")

# DATA GENERATION ---------------------------------------------------------

context("shift: Data generation")

N <- sample(10:20, 1)
J <- sample(5:20, 1)

betas1 <- mpt2irt:::gen_betas(genModel = "shift", J = J)

cond2 <- FALSE
while (cond2 == FALSE) {
    dat1 <- generate_irtree_shift(N = N,
                                  J = J,
                                  betas = betas1,
                                  prop.rev = NULL,
                                  revItem = c(sample(c(0,0,1,1)),
                                              sample(0:1, J-4, T)),
                                  genModel = "shift")
    cond1 <- suppressWarnings(cor(dat1$X)) %>% is.na %>% any %>% magrittr::equals(FALSE)
    if (cond1 == TRUE) {
        cond2 <- dat1$X %>% cor %>% sign %>% magrittr::equals(-1) %>% any
    }
}

suppressWarnings(rm(cond1, cond2))

test_that("generate_irtree_shift() returns correct output", {
    expect_is(dat1, "list")
    expect_equal(ncol(dat1$X), J)
    expect_equal(nrow(dat1$X), N)
})

# MODEL FITTING -----------------------------------------------------------

context("shift: Model fitting")

M <- 200
warmup <- 200

invisible(capture.output(
    res1 <- fit_irtree(dat1$X, fitModel = "shift", fitMethod = "stan",
                       revItem = dat1$revItem, traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1),
    res2 <- fit_irtree(dat1$X, fitModel = "shift", fitMethod = "jags",
                       revItem = dat1$revItem, traitItem = dat1$traitItem,
                       M = M, warmup = warmup, n.chains = 1
                       # n.chains = 2,
                       # method = "simple",
                       # add2varlist = c("deviance", "pd", "popt", "dic")
                       )
))

test_that("fit_irtree() returns MCMC list", {
    expect_equal(unique(sapply(rstan::As.mcmc.list(res1$samples), nrow)), M)
    expect_equal(unique(sapply(res2$samples$mcmc, nrow)), M)
})

# SUMMARIZING MODEL RESULTS -----------------------------------------------

context("shift: Summarizing fitted models")

res1b <- summarize_irtree_fit(res1)
res1c <- tidyup_irtree_fit(res1b)

res2b <- summarize_irtree_fit(res2)
res2c <- tidyup_irtree_fit(res2b)

test_that("tidyup_irtree_fit() returns correlations", {
    expect_equal(unique(as.vector(sapply(res1c$Corr, dim))), res1b$args$S)
    expect_equal(unique(as.vector(sapply(res1c$Sigma, dim))), res1b$args$S)
    expect_equal(unique(as.vector(sapply(res2c$Corr, dim))), res2b$args$S)
    expect_equal(unique(as.vector(sapply(res2c$Sigma, dim))), res2b$args$S)
})

test_that("tidyup_irtree_fit() returns correct number of parameters", {
    expect_equal(unique(sapply(res1c$beta, nrow)), J)
    expect_equal(unique(sapply(res1c$theta, nrow)), N)
    expect_equal(unique(sapply(res2c$beta, nrow)), J)
    expect_equal(unique(sapply(res2c$theta, nrow)), N)
})

test_that("plot_irtree() returns a valid ggplot", {
    expect_is(res1c$plot, "ggplot")
    expect_is(res2c$plot, "ggplot")
    # expect_true(ggplot2::is.ggplot(res2c$plot))
    # expect_is(res2c$plot, "ggplot")
    # expect_error(suppressMessages(suppressWarnings(print(res2c$plot))), NA)
})

# RECOVERY ----------------------------------------------------------------

# context("shift: Recovery")
# 
# test_that("Check that true model parameters are correctly recovered", {
#     
#     cor11 <- cor(dat1$theta, res1c$theta$Median)
#     cor12 <- cor(dat1$betas,  res1c$beta$Median)
#     
#     cor21 <- cor(dat1$theta, res2c$theta$Median)
#     cor22 <- cor(dat1$betas,  res2c$beta$Median)
#     
#     cor31 <- cor(res1c$theta$Median, res2c$theta$Median)
#     cor32 <- cor(res1c$beta$Median,  res2c$beta$Median)
#     
#     # expect_true(all(c(1, 5, 9) %in% tail(order(abs(cor2)), 3)),
#     #             label = "Correlations of true and observed betas show expected pattern")
#     
#     expect_gte(min(diag(cor31)), .8,
#               label = "Comparing Stan and JAGS, the minimum correlation of thetas")
#     expect_gte(min(diag(cor32)), .8,
#                label = "Comparing Stan and JAGS, the minimum correlation of beta")
#     
# })

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

test_that("ppc_resp_irtree() returns valid values", {
    expect_is(res1f, "matrix")
    expect_is(res2f, "matrix")
    
    expect_gte(min(res1f), 0)
    expect_gte(min(res2f), 0)
    
    expect_lte(min(res1f), 1)
    expect_lte(min(res2f), 1)
})

test_that("print(ppc_irtree()) returns valid values", {
    expect_is(res1g, "data.frame")
    expect_is(res2g, "data.frame")
    
    expect_equal(as.numeric(unique(res1g$Item)), 1:J)
    expect_equal(as.numeric(unique(res2g$Item)), 1:J)
    
    expect_equal(as.numeric(unique(res1g$Categ)), 1:5)
    expect_equal(as.numeric(unique(res2g$Categ)), 1:5)
    
    expect_equal(unique(res1g$Persons), N)
    expect_equal(unique(res2g$Persons), N)
    
    expect_gte(min(subset(res1g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    expect_gte(min(subset(res2g, select = c(Obs, q025, q975, q16, q84, q50))), 0)
    
    expect_lte(min(subset(res1g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
    expect_lte(min(subset(res2g, select = c(Obs, q025, q975, q16, q84, q50))), 1)
})
