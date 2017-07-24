#include "license.stan" // GPL3+

// ----- Partial Credit Model -----

// 1,...j,...,J: number of items
// 1,...i,...,N: number of participants
// S: number of dimension (== max(traitItem))
// X: N(person) x J(item) x C(categories) array of observed frequencies
// X: N(person) x J(item) matrix with observed response category

// ----- Hyperpriors
// df: degrees of freedom for scaled inverse wishart (>=S, typically df=S+1)
// V: S x S hyperprior for wishart

functions {
    // http://mc-stan.org/documentation/case-studies/pcm_latent_reg.html
    // http://mc-stan.org/documentation/case-studies.html
    vector pcm_probs(real theta1, row_vector beta1, real revx) {
        vector[cols(beta1) + 1] unsummed;
        vector[cols(beta1) + 1] probs;
        vector[cols(beta1) + 1] probsx;
        unsummed = append_row(rep_vector(0.0, 1), theta1 - to_vector(beta1));
        probs = softmax(cumulative_sum(unsummed));
        if (revx == 0) {
            probsx = probs;
        } else {
            for (p in 1:5) {
                probsx[p] = probs[6-p];
            }
        }
        return probsx;
    }
}

data {
	// data and indices
	int<lower=1> N;  					// number of persons
	int<lower=1> J;  					// number of items
	int<lower=1, upper=5> X[N, J];	    // chosen responses of partipants
	int<lower=1> S;						// number of theta-parameters (2012-version: S=3
	int<lower=0, upper=1> revItem[J];   // index for reversed items (=1)
	int<lower=1> traitItem[J];   		// index for trait items (1,...,n.trait)
	int<lower=1> N2;  					// number of persons for whom to draw posterior predictives
	
	// hyperpriors
	vector[S] theta_mu;  // mean of theta
	int<lower=S+1> df;   // df for wishart, often S+1
	cov_matrix[S] V;     // hyperprior for wishart, e.g., diag(S)
}

// ----------------------------------------
parameters {
    // person parameters
    vector[S] theta_raw[N];                      // unscaled latent trait values
    vector<lower=0, upper=100>[S] xi_theta;      // scaling parameters
    cov_matrix[S] Sigma_raw;                     // unscaled covariance matrix of traits
    
    // item parameters
    matrix<lower=-5, upper=5>[J, 4] beta_raw;    // raw item difficulties
    vector<lower=-5, upper=5>[S] mu_beta;        // raw item means
    vector<lower=0>[S] sigma2_beta_raw;          // raw item variance
} 

// ----------------------------------------
transformed parameters {
    matrix[N, S] theta;                  // latent traits
    cov_matrix[S] Sigma;                 // covariance matrix of traits
    
    matrix[J, 4] beta;     			   // item difficulties
    vector<lower=0>[S] sigma_beta_raw;   // raw item variance
    
    simplex[5] p_cat[N, J];              // response category probabilities
  
    // scaling of variance
    Sigma = diag_matrix(xi_theta) * Sigma_raw* diag_matrix(xi_theta);

    // print("log-posterior = ", target());

    // ----- rescaling of item parameters
    for(s in 1:S){
    	sigma_beta_raw[s] = sqrt(sigma2_beta_raw[s]);
    }
    for(j in 1:J){
    	beta[j, 1:4] =  mu_beta[traitItem[j]] + beta_raw[j, 1:4];
    }
    
     for(i in 1:N){	
    	// rescale trait values
    	for(s in 1:S){
    		theta[i,s] = theta_raw[i,s] * xi_theta[s];
    	}
    	
    	// loop across items
    	for(j in 1:J){
    		// response probabilities in five categories
		    p_cat[i, j, 1:5] = pcm_probs(theta[i, traitItem[j]], beta[j, 1:4], revItem[j]);
        }
    }
}

model {
    // ----- independent univariate normals for item difficulties: 
    for(j in 1:J){
    	beta_raw[j, 1:4] ~ normal(0, sigma_beta_raw[traitItem[j]]);
    }
    
    // ----- hyperpriors:
    // implicit uniform on scaling parameters
    mu_beta ~ normal(0, 1);            // raw item mean
    sigma2_beta_raw ~ inv_gamma(1,1);  // raw item variance
    Sigma_raw ~ inv_wishart(df, V);    // person hyperprior
    
    for(i in 1:N){
    	for(j in 1:J){
    		// distribution of observed frequencies
    		X[i,j] ~ categorical(p_cat[i,j]);      // categorical data dim(X)= N x J
    	}
    	
    	// Hierarchical model for participant parameters
    	// fix person mean to zero for weak identification
    	theta_raw[i] ~ multi_normal(theta_mu, Sigma_raw); 
    }
}

// ----- posterior predictive // -----

generated quantities {
    int<lower=1, upper=5> X_pred[N2, J];	    // predicted responses of partipants

	vector<lower=0>[S] sigma_beta;	            // item SD
    sigma_beta = rep_vector(1, S) .* sigma_beta_raw;

    for(i in 1:N2){
    	for(j in 1:J){
    	    X_pred[i,j] = categorical_rng(p_cat[i,j]);
    	}
    }
}
