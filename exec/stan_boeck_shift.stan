#include "license.stan" // GPL3+

// ----- Boeckenholt 2012 plus ARS shift -----

// Target trait parameter: t_ij = Phi(theta_ti + theta_ai - beta_tj) (i.e., trait + ars)

// 1,...j,...,J: number of items
// 1,...i,...,N: number of participants
// S: number of parameters (extended Boeckenholt: 4: middle/trait/extremity/acquies)
// X: N(person) x J(item) x C(categories) array of observed frequencies
// X: N(person) x J(item) matrix with observed response category

// ----- Hyperpriors
// df: degrees of freedom for scaled inverse wishart (>=S, typically df=S+1)
// V: S x S hyperprior for wishart

// ----------------------------------------
data {
	// data and indices
	int<lower=1> N;  					 // number of persons
	int<lower=1> J;  					 // number of items
	int<lower=1, upper=5> X[N, J];	     // chosen responses of partipants
	int<lower=1> S;						 // number of theta-parameters (Shift: S=4)
	int<lower=0, upper=1> revItem[J];    // index for reversed items (=1)
	int<lower=1> traitItem[J];   		 // index for trait items (1,...,n.trait)
	int<lower=1> N2;  					 // number of persons for whom to draw posterior predictives
	
	// hyperpriors
	vector[S] theta_mu;  // mean of theta
	int<lower=S+1> df;   // df for wishart, often S+1
	cov_matrix[S] V;     // hyperprior for wishart, e.g., diag(S)
}

// ----------------------------------------
parameters {
    // person parameters: middle, trait, extreme
    vector[S] theta_raw[N];                      // unscaled latent trait values
    vector<lower=0, upper=100>[S] xi_theta;      // scaling parameters
    cov_matrix[S] Sigma_raw;                     // unscaled covariance matrix of traits
    
    // item parameters: middle, extreme, trait
    matrix<lower=-5, upper=5>[J, 3] beta_raw;    // raw item difficulties
    // vector<lower=0, upper=100>[S] xi_beta;       // scaling parameters items
    vector<lower=-5, upper=5>[S-1] mu_beta;      // raw item means
    vector<lower=0>[S-1] sigma2_beta_raw;        // raw item variance
} 

// ----------------------------------------
transformed parameters {
    matrix[N, S] theta;                  // latent traits
    cov_matrix[S] Sigma;                 // covariance matrix of traits
    
    matrix[J, 3] beta;     			     // item difficulties
    vector<lower=0>[S-1] sigma_beta_raw; // raw item variance
    
    simplex[5] p_cat[N, J];              // response category probabilities
    real<lower=0, upper=1> middle[N,J];  // item-person probability to select middle category
    real<lower=0, upper=1> extreme[N,J]; // item-person probability to respond extremely
    real<lower=0, upper=1> trait[N,J];   // item-person probability to respond positive/negative
  
    // scaling of variance
    Sigma = diag_matrix(xi_theta) * Sigma_raw * diag_matrix(xi_theta);
    
    // ----- rescaling of item parameters
    for(s in 1:(S-1)){
    	sigma_beta_raw[s] = sqrt(sigma2_beta_raw[s]);
    }
    for(j in 1:J){
    	// Response styles: MRS and ERS
    	for(s in 1:2){
    		// independent univariate normal 
    		// beta[j,s] = mu_beta[s] + xi_beta[s] * beta_raw[j,s];
    		beta[j,s] = mu_beta[s] + beta_raw[j,s];
    	}
    	// Trait(s):
    	// beta[j,3] =  mu_beta[3+traitItem[j]]+ xi_beta[2+traitItem[j]] * beta_raw[j,3];
    	beta[j,3] =  mu_beta[2 + traitItem[j]] + beta_raw[j, 3];
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
    		middle[i, j]  = Phi_approx(theta[i,1] - beta[j,1]);   // Phi_approx()
    		extreme[i, j] = Phi_approx(theta[i,2] - beta[j,2]);
    		// standard items: theta-beta  // reversed items: beta-theta
            trait[i, j]   = Phi_approx(theta[i,3] + (-1)^revItem[j] *(theta[i,3+traitItem[j]]-beta[j,3]) );
    		
    		// response probabilities: MPT model for response categories 1, 2, 3, 4, 5
    		p_cat[i,j,1] = (1-middle[i,j])*(1-trait[i,j])*   extreme[i,j];
    		p_cat[i,j,2] = (1-middle[i,j])*(1-trait[i,j])*(1-extreme[i,j]);
    		p_cat[i,j,3] =    middle[i,j];
    		p_cat[i,j,4] = (1-middle[i,j])*    trait[i,j]*(1-extreme[i,j]);
    		p_cat[i,j,5] = (1-middle[i,j])*    trait[i,j]*   extreme[i,j];
    	}
    }
}
		

// ----------------------------------------
model {
    // ----- independent univariate normals for item difficulties: 
    for(j in 1:J){
    	// MRS, ERS, ARS:
    	for(s in 1:2){
            beta_raw[j,s] ~ normal(0, sigma_beta_raw[s]);
    	}
    	// Trait(s):
    	// if (traitItem[j] > 0) {
    	    beta_raw[j,3] ~ normal(0, sigma_beta_raw[2 + traitItem[j]]);
    	// }
	}
	
    // ----- hyperpriors:
    // implicit uniform on scaling parameters
    // xi_theta, xi_beta ~ uniform (0, 100);     
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

generated quantities {
    int<lower=1, upper=5> X_pred[N2, J];	  // predicted responses of partipants
    
	vector<lower=0>[S-1] sigma_beta;	      // item SD
	
	cov_matrix[S] Corr;

    sigma_beta = rep_vector(1, (S-1)) .* sigma_beta_raw;
    
    Corr = diag_matrix(inv_sqrt(diagonal(Sigma))) * Sigma * diag_matrix(inv_sqrt(diagonal(Sigma)));

    // ----- POSTERIOR PREDICTIVE
    for(i in 1:N2){
    	for(j in 1:J){
    	    X_pred[i,j] = categorical_rng(p_cat[i,j]);
    	}
    }
}
