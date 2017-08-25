#include "license.stan" // GPL3+

// ----- Steps Model aka Sequential Model (Verhelst, Tutz) -----
// 1 2 3 4 5   # raw data
//             #  model as
// 0 1 1 1 1   # step1
// . 0 1 1 1   # step2
// . . 0 1 1   # step3
// . . . 0 1   # step4

// 1,...j,...,J: number of items
// 1,...i,...,N: number of participants
// S: number of parameters (standard Boeckenholt: 3: middle/extremity/trait)
// X: N(person) x J(item) x C(categories) array of observed frequencies
// X: N(person) x J(item) matrix with observed response category
// revItem: vector of length J indicating reversed items (1=reversed)

// ----- Hyperpriors
// df: degrees of freedom for scaled inverse wishart (>=S, typically df=S+1)
// V: S x S hyperprior for wishart

// ----------------------------------------
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
  // person parameters: middle, trait, extreme
  vector[S] theta_raw[N];        			  // unscaled latent trait values
  vector<lower=0, upper=100>[S] xi_theta;     // scaling parameters persons
  cov_matrix[S] Sigma_raw;       			  // unscaled covariance matrix of traits
  
  // item parameters: middle, trait, extreme
  // matrix[J, 4] beta_raw;      				  // raw item difficulties
  matrix<lower=-5, upper=5>[J, 4] beta_raw;   // raw item difficulties
  // vector<lower=0, upper=100>[S] xi_beta;      // scaling parameters items
  // vector[S] mu_beta;      				      // item means
  // vector<lower=-5, upper=5>[S] mu_beta;       // item means
  vector<lower=-5, upper=5>[S*4] mu_beta_vec;       // item means
  // vector<lower=0>[S] sigma2_beta_raw;    	  // raw item variance
  vector<lower=0>[S*4] sigma2_beta_raw;    	  // raw item variance
  
} 

// ----------------------------------------
transformed parameters {
    matrix[N, S] theta;                  // latent traits
    cov_matrix[S] Sigma;                 // covariance matrix of traits
    
    matrix[J, 4] beta;     			     // item difficulties
    // vector<lower=0>[S] sigma_beta_raw;   // raw item variance
    vector<lower=0>[S*4] sigma_beta_raw;   // raw item variance
    
    simplex[5] p_cat[N, J];              // response category probabilities
    real<lower=0, upper=1> node1[N,J];   // item-person probability for 1st node
    real<lower=0, upper=1> node2[N,J];   // item-person probability for 2nd node
    real<lower=0, upper=1> node3[N,J];   // item-person probability for 3rd node
    real<lower=0, upper=1> node4[N,J];   // item-person probability for 4th node

    // scaling of variance
    Sigma = diag_matrix(xi_theta) * Sigma_raw * diag_matrix(xi_theta);
    
    // ----- rescaling of item parameters
    // for(s in 1:S){
    for(s in 1:(S*4)){
    	sigma_beta_raw[s] = sqrt(sigma2_beta_raw[s]);
    }
    for(j in 1:J){
    	// beta[j, 1:4] = mu_beta[traitItem[j]] + beta_raw[j, 1:4];
	    beta[j, 1:4] = to_row_vector(mu_beta_vec[(1+4*(traitItem[j]-1)):(4+4*(traitItem[j]-1))]) +
	        beta_raw[j, 1:4];
    }
    
    for(i in 1:N){	
        // rescale trait values
    	for(s in 1:S){
    		theta[i,s] = theta_raw[i,s] * xi_theta[s];
    	}
    	
    	// loop across items
    	for(j in 1:J){
    		// IRT: 1. additive combination of trait and difficulty
    		//      2. transform to probabilities using Phi (trait: consider reversed items!)
    		// node1[i,j] = Phi_approx((-1)^revItem[j]*(theta[i, traitItem[j]] - beta[j, 1]));
    		// node2[i,j] = Phi_approx((-1)^revItem[j]*(theta[i, traitItem[j]] - beta[j, 2]));
    		// node3[i,j] = Phi_approx((-1)^revItem[j]*(theta[i, traitItem[j]] - beta[j, 3]));
    		// node4[i,j] = Phi_approx((-1)^revItem[j]*(theta[i, traitItem[j]] - beta[j, 4]));
    		
    		node1[i,j] = Phi_approx(theta[i, traitItem[j]] - beta[j, 1]);
            node2[i,j] = Phi_approx(theta[i, traitItem[j]] - beta[j, 2]);
            node3[i,j] = Phi_approx(theta[i, traitItem[j]] - beta[j, 3]);
            node4[i,j] = Phi_approx(theta[i, traitItem[j]] - beta[j, 4]);
    
    		// response probabilities: MPT model for response categories 1, 2, 3, 4, 5
    		if (revItem[j] == 0) {
    		    p_cat[i,j,1] = (1-node1[i, j]);
                p_cat[i,j,2] =    node1[i, j] *(1-node2[i, j]);
                p_cat[i,j,3] =    node1[i, j] *   node2[i, j] *(1-node3[i, j]);
                p_cat[i,j,4] =    node1[i, j] *   node2[i, j] *   node3[i, j] *(1-node4[i, j]);
                p_cat[i,j,5] =    node1[i, j] *   node2[i, j] *   node3[i, j] *   node4[i, j] ;
    		} else {
    		    p_cat[i,j,5] = (1-node1[i, j]);
                p_cat[i,j,4] =    node1[i, j] *(1-node2[i, j]);
                p_cat[i,j,3] =    node1[i, j] *   node2[i, j] *(1-node3[i, j]);
                p_cat[i,j,2] =    node1[i, j] *   node2[i, j] *   node3[i, j] *(1-node4[i, j]);
                p_cat[i,j,1] =    node1[i, j] *   node2[i, j] *   node3[i, j] *   node4[i, j] ;
    		}
        }
    }
}

// ----------------------------------------
model {
    // ----- independent univariate normals for item difficulties: 
    for(j in 1:J){
        // beta_raw[j, 1:4] ~ normal(0, sigma_beta_raw[traitItem[j]]);
        beta_raw[j, 1:4] ~ normal(0, 
                                  sigma_beta_raw[(1+4*(traitItem[j]-1)):(4+4*(traitItem[j]-1))]
                                  );
    }
    
    // ----- hyperpriors:
    // implicit uniform on scaling parameters
    // xi_theta, xi_beta ~ uniform (0, 100);     
    mu_beta_vec ~ normal(0, 1);            // raw item mean
    sigma2_beta_raw ~ inv_gamma(1,1);  // raw item variance
    Sigma_raw ~ inv_wishart(df, V);    // person hyperprior
    
    for(i in 1:N){
    	for(j in 1:J){
    		// distribution of observed frequencies
    		X[i,j] ~ categorical(p_cat[i,j]);      // categorical data dim(X)= N x J
    	}
    	
    	// Hierarchical model for participant parameters
    	// fix person mean to zero for weak identification
    	if (S > 1) {
    	    theta_raw[i] ~ multi_normal(theta_mu, Sigma_raw); 
    	} else {
    	    theta_raw[i] ~ normal(theta_mu, Sigma_raw[1, 1]); 
    	}
    }
}

// ----- posterior predictive
generated quantities {
    cov_matrix[S] Corr;                     // Correlation matrix
    // vector<lower=0>[S] sigma_beta;	        // item SD
    // vector<lower=0>[S*4] sigma_beta;	        // item SD
    matrix<lower=0>[S,4] sigma_beta;	        // item SD
    matrix[S,4] mu_beta;	        // item SD
    int<lower=1, upper=5> X_pred[N2, J];    // predicted responses of partipants
    
    Corr = diag_matrix(inv_sqrt(diagonal(Sigma))) * Sigma * diag_matrix(inv_sqrt(diagonal(Sigma)));

    // sigma_beta = rep_vector(1, S) .* sigma_beta_raw;
    // sigma_beta = rep_vector(1, S*4) .* sigma_beta_raw;
    for(s in 1:S){
        sigma_beta[s, 1:4] = to_row_vector(rep_vector(1, 4)) .*
            to_row_vector(sigma_beta_raw[(1+(s-1)*4):(4+(s-1)*4)]);
        mu_beta[s, 1:4] = to_row_vector(mu_beta_vec[(1+(s-1)*4):(4+(s-1)*4)]);
    }
    
    for(i in 1:N2){
        for(j in 1:J){
            X_pred[i,j] = categorical_rng(p_cat[i,j]);
     	}
    }
}
