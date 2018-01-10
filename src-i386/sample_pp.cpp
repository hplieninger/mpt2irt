#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector sample_pp(NumericMatrix prob){
    int M = prob.rows();
    NumericVector pp(M);
    NumericVector p;
    for(int i=0; i<M; i++){
        p = prob.row(i);
        pp(i) = sample(5, 1, false, p)(0);
    }
    return pp;
}

/*** R
# prob <- matrix(runif(5*20), 20)
# sample_pp(prob)
*/
