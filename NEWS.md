## Version 0.2.0

* The Stan part of the package was split off into a new package mpt2irtStan.
This was done (a) to comply with new versions of rstan and rstantools, (b) to be able to have C++ code (`sample_pp()`) in mpt2irt, and (c) for easier maintenance.

## Version 0.1.4

* Fixed bug in estimation of latent variance Sigma for uni-dimensional models
(i.e., Steps and PCM)

## Version 0.1.3

* Added JAGS version of the shift model (introduced only in Stan in 0.1.1)

## Version 0.1.2

* Implementation of more thorough posterior predictive checking including
several discrepancy measures.

## Version 0.1.1

* Implementation of the Steps Model (Verhelst; Tutz) as well as a shift model
(both only in Stan)

* Bug fix with respect to the priors for the item parameters in the PCM

* Core functions now coherently return a list of arguments and setup parameters
in a list called `args`.

## Version 0.1.0

* Basic functionality via `fit_irtree()` to fit the response style model of
Boeckenholt (2012) and an extended version including acquiescence using JAGS
and/or Stan

* Recovery simulations using `recovery_irtree()`

* Functions to plot and summarize a fitted model
