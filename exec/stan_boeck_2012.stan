#include "license.stan" // GPL3+

############################### Boeckenholt (2012): Response styles ##########################

# 1,...j,...,J: number of items
# 1,...i,...,N: number of participants
# S: number of parameters (standard Boeckenholt: 3: middle/extremity/trait)
# X: N(person) x J(item) x C(categories) array of observed frequencies
# X: N(person) x J(item) matrix with observed response category
# revItem: vector of length J indicating reversed items (1=reversed)

#### Hyperpriors
# df: degrees of freedom for scaled inverse wishart (>=S, typically df=S+1)
# V: S x S hyperprior for wishart

#####################################################
data {
	# data and indices
	int<lower=1> N;  					# number of persons
	int<lower=1> J;  					# number of items
	int<lower=1, upper=5> X[N, J];	    # chosen responses of partipants
	int<lower=1> S;						# number of theta-parameters (2012-version: S=3
	int<lower=0, upper=1> revItem[J];   # index for reversed items (=1)
	int<lower=1> traitItem[J];   		# index for trait items (1,...,n.trait)
	# real<lower=0> T1_CONST;             # constant added to mean-expected frequencies of zero
	int<lower=1> N2;  					 # number of persons for whom to draw posterior predictives
	
	# hyperpriors
	vector[S] theta_mu;  # mean of theta
	int<lower=S+1> df;   # df for wishart, often S+1
	cov_matrix[S] V;     # hyperprior for wishart, e.g., diag(S)
}

#####################################################
parameters {
  # person parameters: middle, trait, extreme
  vector[S] theta_raw[N];        			 # unscaled latent trait values
  vector<lower=0, upper=100>[S] xi_theta;    # scaling parameters persons
  cov_matrix[S] Sigma_raw;       			 # unscaled covariance matrix of traits
  
  # item parameters: middle, trait, extreme
  // matrix[J, 3] beta_raw;      				# raw item difficulties
  matrix<lower=-5, upper=5>[J, 3] beta_raw;      				# raw item difficulties
  // vector<lower=0, upper=100>[S] xi_beta;    # scaling parameters items
  // vector[S] mu_beta;      				# item means
  vector<lower=-5, upper=5>[S] mu_beta;      				# item means
  vector<lower=0>[S] sigma2_beta_raw;    		    # raw item variance
  
} 

#####################################################
transformed parameters {
  matrix[N, S] theta;                  # latent traits
  cov_matrix[S] Sigma;                 # covariance matrix of traits
  
  matrix[J, 3] beta;     			   # item difficulties
  vector<lower=0>[S] sigma_beta_raw;   # raw item variance

  simplex[5] p_cat[N, J];              # response category probabilities
  real<lower=0, upper=1> middle[N,J];  # item-person probability to select middle category
  real<lower=0, upper=1> extreme[N,J]; # item-person probability to respond extremely
  real<lower=0, upper=1> trait[N,J];   # item-person probability to respond positive/negative
  

# scaling of variance
Sigma = diag_matrix(xi_theta) * Sigma_raw* diag_matrix(xi_theta);

#  print("lp before =",get_lp());

### rescaling of item parameters
for(s in 1:S){
	sigma_beta_raw[s] = sqrt(sigma2_beta_raw[s]);
}
for(j in 1:J){
	## Response styles: MRS and ERS
	for(s in 1:2){
		 # independent univariate normal 
		// beta[j,s] = mu_beta[s] + xi_beta[s] * beta_raw[j,s];
		beta[j,s] = mu_beta[s] + beta_raw[j,s];
	}
	# Trait(s):
	beta[j,3] = mu_beta[2+traitItem[j]] + beta_raw[j,3];
}


 for(i in 1:N){	
	# rescale trait values
	for(s in 1:S){
		theta[i,s] = theta_raw[i,s] * xi_theta[s];
	}
	
	# loop across items
	for(j in 1:J){
		# IRT: 1. additive combination of trait and difficulty
		#      2. transform to probabilities using Phi (trait: consider reversed items!)
		middle[i,j] = Phi_approx(theta[i,1]-beta[j,1]);   #Phi_approx()
		extreme[i,j] =  Phi_approx(theta[i,2]-beta[j,2]);
		# standard items: theta-beta  // reversed items: beta-theta
		trait[i,j] =  Phi_approx( (-1)^revItem[j] *(theta[i,2+traitItem[j]]-beta[j,3]) );
				
		# response probabilities: MPT model for response categories 1, 2, 3, 4, 5
		p_cat[i,j,1] = (1-middle[i,j])*(1-trait[i,j])*extreme[i,j];
		p_cat[i,j,2] = (1-middle[i,j])*(1-trait[i,j])*(1-extreme[i,j]);
		p_cat[i,j,3] =    middle[i,j];
		p_cat[i,j,4] = (1-middle[i,j])*trait[i,j]*(1-extreme[i,j]);
		p_cat[i,j,5] = (1-middle[i,j])*trait[i,j]*extreme[i,j];
  }
}

}

#####################################################
model {


######## independent univariate normals for item difficulties: 
for(j in 1:J){
	## MRS, ERS:
	for(s in 1:2){
		beta_raw[j,s] ~ normal(0, sigma_beta_raw[s]);
	}
	## Trait(s):
	beta_raw[j,3] ~ normal(0, sigma_beta_raw[2+traitItem[j]]);
}

#######  hyperpriors:
# implicit uniform on scaling parameters
# xi_theta, xi_beta ~ uniform (0, 100);     
mu_beta ~ normal(0, 1);      # raw item mean
sigma2_beta_raw ~ inv_gamma(1,1);  # raw item variance
Sigma_raw ~ inv_wishart(df, V);    # person hyperprior

for(i in 1:N){
	for(j in 1:J){
		# distribution of observed frequencies
		X[i,j] ~ categorical(p_cat[i,j]);      # categorical data dim(X)= N x J
	}
	
	# Hierarchical model for participant parameters
	# fix person mean to zero for weak identification
	theta_raw[i] ~ multi_normal(theta_mu, Sigma_raw); 
}

}

##################################################### posterior predictive
generated quantities {
    int<lower=1, upper=5> X_pred[N2, J];	    # predicted responses of partipants
    
    // matrix[N,J] dev_obs;
    // matrix[N,J] dev_pred;
    // real<lower=0> T_obs;
    // real<lower=0> T_pred;
    // int<lower=0, upper=1> post_p;
    vector<lower=0>[S] sigma_beta;	 # item SD
    // vector[S] mu_beta;	             #  item mean

    // # rescaled item mean and SD
    // mu_beta = xi_beta .* mu_beta_raw;
    // sigma_beta = xi_beta .* sigma_beta_raw;
    sigma_beta = rep_vector(1, S) .* sigma_beta_raw;

    for(i in 1:N2){
        for(j in 1:J){
            X_pred[i,j] = categorical_rng(p_cat[i,j]);
     		
     		// dev_obs[i,j] = ((p_cat[i,j, X[i,j]]-1) / p_cat[i,j, X[i,j]])^2;
     		// dev_pred[i,j] = ((p_cat[i,j, X_pred[i,j]]-1) / p_cat[i,j, X_pred[i,j]])^2;
     	}
    }
    
    // # posterior predictive checks
    // T_obs = sum(dev_obs);
    // T_pred = sum(dev_pred);
    // post_p = T_pred > T_obs;
}
