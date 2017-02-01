library(mpt2irt)
context("Data generation")

J <- sample(1:100, 1)
N <- sample(1:100, 1)

betas1 <- mpt2irt:::gen_betas(genModel = "ext", J = J)
betas2 <- mpt2irt:::gen_betas(genModel = "2012", J = J)

dat1 <- generate_irtree_ext(N = N,
                            J = J,
                            betas = betas1,
                            prop.rev = runif(1),
                            genModel = "ext",
                            beta_ARS_extreme = rnorm(1))
dat2 <- generate_irtree_2012(N = N,
                             J = J,
                             betas = betas2,
                             prop.rev = runif(1))

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
