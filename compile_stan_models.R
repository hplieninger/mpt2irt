

# library("rstan")

filename <-  "./inst/models/stan_boeck_2012.stan"
stan.txt <- readChar(filename, file.info(filename)$size)
boeck_stan_2012 = rstan::stan_model( model_code=stan.txt, model_name="Boeckenholt_2012")

filename <- "./inst/models/stan_boeck_ext.stan"
stan.txt <- readChar(filename, file.info(filename)$size)
boeck_stan_ext = rstan::stan_model( model_code=stan.txt, model_name="Boeckenholt_ext")

# filename <-  "inst/models/stan_boeck_ext2.stan"
# stan.txt <- readChar(filename, file.info(filename)$size)
# boeck_stan_ext2 = rstan::stan_model( model_code=stan.txt, model_name="Boeckenholt_ext2")
# 
# filename <-  "inst/models/stan_boeck_ext3.stan"
# stan.txt <- readChar(filename, file.info(filename)$size)
# boeck_stan_ext3 = rstan::stan_model( model_code=stan.txt, model_name="Boeckenholt_ext3")

# filename <-  "inst/models/stan_boeck_ext5.stan"
# stan.txt <- readChar(filename, file.info(filename)$size)
# boeck_stan_ext5 = rstan::stan_model( model_code=stan.txt, model_name="Boeckenholt_ext5")

filename <- "./inst/models/stan_boeck_ext_HH.stan"
stan.txt <- readChar(filename, file.info(filename)$size)
boeck_stan_ext_HH = rstan::stan_model(model_code = stan.txt, model_name = "Boeckenholt_ext_HH")

filename <- "./inst/models/stan_pcm.stan"
stan.txt <- readChar(filename, file.info(filename)$size)
boeck_stan_pcm = rstan::stan_model(model_code = stan.txt, model_name = "Boeckenholt_PCM")


# save(boeck_stan_ext, boeck_stan_ext2, boeck_stan_2012, file=paste0('data/boeck_stan_models.rda'))
# save(boeck_stan_ext, boeck_stan_ext2, boeck_stan_ext3, boeck_stan_2012, file=paste0('data/boeck_stan_models.rda'))
# save(boeck_stan_ext, boeck_stan_2012, file=paste0('data/boeck_stan_models.rda'))
# save(boeck_stan_ext, boeck_stan_2012, boeck_stan_ext_HH, file=paste0('./R/sysdata.rda'))
# devtools::use_data(boeck_stan_ext, boeck_stan_2012, boeck_stan_ext_HH, internal = TRUE, overwrite = TRUE)
devtools::use_data(boeck_stan_ext, boeck_stan_2012, boeck_stan_ext_HH,
                   boeck_stan_pcm,
                   # boeck_stan_ext5,
                   internal = TRUE, overwrite = TRUE)

# save(boeck_stan_ext, boeck_stan_2012, file=paste0('./R/sysdata.rda'))
devtools::use_data(boeck_stan_ext, boeck_stan_2012, boeck_stan_pcm,
                   internal = TRUE, overwrite = TRUE)



# DESCRIPTION file --------------------------------------------------------

# packages_imports <- sort(c("runjags",
#                            # "rjags",
#                            "coda",
#                            "rstan",
#                            "doParallel",
#                            "parallel",
#                            "MASS",
#                            "foreach",
#                            # "ggplot2",
#                            "mail",
#                            "MBESS",
#                            "truncnorm",
#                            "checkmate"))
# for (iii in seq_along(packages_imports)) {
#     devtools::use_package(packages_imports[iii])
# }
# packages_suggests <- sort(c("ggplot2",
#                            "mail",
#                            "reshape2",
#                            "dplyr",
#                            "tidyr"))
# for (iii in seq_along(packages_suggests)) {
#     devtools::use_package(packages_suggests[iii], "Suggests")
# }
