#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4stan_boeck_2012_mod) {


    class_<rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> >("model_stan_boeck_2012")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_stan_boeck_2012_namespace::model_stan_boeck_2012, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4stan_boeck_ext_mod) {


    class_<rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> >("model_stan_boeck_ext")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_stan_boeck_ext_namespace::model_stan_boeck_ext, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4stan_boeck_ext_HH_mod) {


    class_<rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> >("model_stan_boeck_ext_HH")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_stan_boeck_ext_HH_namespace::model_stan_boeck_ext_HH, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4stan_boeck_shift_mod) {


    class_<rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> >("model_stan_boeck_shift")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_stan_boeck_shift_namespace::model_stan_boeck_shift, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4stan_pcm_mod) {


    class_<rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> >("model_stan_pcm")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_stan_pcm_namespace::model_stan_pcm, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4stan_steps_mod) {


    class_<rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> >("model_stan_steps")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_stan_steps_namespace::model_stan_steps, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
