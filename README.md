
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.org/hplieninger/mpt2irt.svg?branch=master)](https://travis-ci.org/hplieninger/mpt2irt)

# mpt2irt

mpt2irt is an R package that accompanies the paper *A new model for
acquiescence at the interface of psychometrics and cognitive psychology*
[(Plieninger &
Heck, 2018)](https://doi.org/10.1080/00273171.2018.1469966). Therein, we
extend the response style model of Böckenholt (2012) to acquiescence.
The model is essentially a hierarchical multinomial processing tree
(MPT) model with an item response theory (IRT) structure of its
parameters. To estimate the model parameters, we build on Bayesian
hierarchical modeling and fit the model in either Stan or JAGS.

## Installation

In order to use the package, you will need either JAGS or RStan.

### Install JAGS

To install JAGS, visit <https://sourceforge.net/projects/mcmc-jags/>.

### Install RStan

To install RStan, visit
[https://github.com/stan-dev](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
and carefully follow the instructions. This may also involve the
following steps:

  - [Installing
    Rtools](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#checking-the-c-toolchain)
  - [Configuration of the C++
    Toolchain](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#configuration-of-the-c-toolchain)
  - [Creating or editing a Makevars file on
    Windows](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows#configuration)

### Install mpt2irt

Actually, the Stan part of mpt2irt was split off for easier maintenance
and is provided in the package
[mpt2irtStan](https://github.com/hplieninger/mpt2irtStan/). User have to
install both packages even though they will interface only with mpt2irt.

The package mpt2irt can be installed directly from GitHub, and this
should automatically also install mpt2irtStan:

``` r
# install.packages("remotes")
remotes::install_github("hplieninger/mpt2irt")
```

However, because compiling the code in mpt2irtStan takes a while and may
need a special setup (see above), users are encouraged to first install
mpt2irtStan via `remotes::install_github("hplieninger/mpt2irtStan")`.
During installation, users may be asked to update or install the
rstantools package (version \>= 2.0.0) and should agree. If this was
successful, the main package mpt2irt can be installed
afterwards.

## Usage

``` r
# This is a minimal working example, where data are generated and subsequently fit.
library("mpt2irt")
N <- 20
J <- 10
betas <- cbind(rnorm(J, .5), rnorm(J, .5), rnorm(J, 1.5), rnorm(J, 0))
dat   <- generate_irtree_ext(N = N, J = J, betas = betas, beta_ARS_extreme = .5)

# fit model
res1 <- fit_irtree(dat$X, revItem = dat$revItem, M = 200)
res2 <- summarize_irtree_fit(res1)
res3 <- tidyup_irtree_fit(res2)
res3$plot
```

## Misc

The proposed Acquiescence Model is a mixture model. Existing approaches
to ARS (e.g., Billiet & McClendon, 2000; Ferrando et al., 2016;
Maydeu-Olivares & Coffman, 2006) view acquiescence as a shift process. A
graphical comparison of the two approaches in terms of the predicted
category probabilities may be found at
<https://hplieninger.shinyapps.io/shift-vs-mixture-ARS>.
