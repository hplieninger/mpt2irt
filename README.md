
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

`mpt2irt` is an `R` package that accompanies a manuscript submitted for publication. Therein, we extend the response style model of Böckenholt (2012) to acquiescence. The model is essentially a hierarchical multinomial processing tree (MPT) model with an item response theory (IRT) structure of its parameters. To estimate the model parameters, we build on Bayesian hierarchical modeling and fit the model in either Stan or JAGS.

Usage
=====

``` r
# This is a minimal working example, where data are generated and subsequently fit.
N <- 20
J <- 10
betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
dat   <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)

# fit model
res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200)
res2 <- summarize_irtree_fit(res1)
res3 <- tidyup_irtree_fit(res2, N = N, J = J, revItem = dat$revItem,
                          traitItem = dat$traitItem, fitModel = res$fitModel)
res3$plot
```
