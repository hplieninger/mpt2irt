############################### Boeckenholt (2012): Response styles ##########################

# 1,...j,...,J: number of items
# 1,...i,...,N: number of participants
# S: number of parameters (standard Böckenholt: 3: middle/trait/extremity)
# X: N(person) x J(item) x C(categories) array of observed frequencies
# X: N(person) x J(item) matrix with observed response category
# revItem: vector of length J indicating reversed items (1=reversed)

#### Hyperpriors
# df: degrees of freedom for scaled inverse wishart (>=S, typically df=S+1)
# V: S x S hyperprior for wishart

data {
	# data and indices
	int<lower=1> N;  					# number of persons
	int<lower=1> J;  					# number of items
	int<lower=1, upper=5> X[N, J];	    # chosen responses of partipants
	int<lower=1> S;						# number of theta-parameters (2012-version: S=3
	int<lower=0, upper=1> revItem[J];   # index for reversed items (=1)
	
	# hyperpriors
	vector[S] theta_mu;  # mean of theta
	int<lower=S+1> df;   # df for wishart, often S+1
	cov_matrix[S] V;     # hyperprior for wishart, e.g., diag(S)
}

parameters {
  # person parameters: middle, trait, extreme
  vector[S] theta_raw[N];         # unscaled latent trait values
  vector[S] xi_theta;             # scaling parameters
  cov_matrix[S] Sigma_raw;  # unscaled covariance matrix of traits
  
  # item parameters: middle, trait, extreme
  matrix[J, S] beta;                     # beta estimated freely
} 

transformed parameters {
  simplex[5] p_cat[N, J];              # response category probabilities
  matrix[N, S] theta;                  # latent traits
  real<lower=0, upper=1> middle[N,J];  # item-person probability to select middle category
  real<lower=0, upper=1> trait[N,J];   # item-person probability to respond positive/negative
  real<lower=0, upper=1> extreme[N,J]; # item-person probability to respond extremely
  cov_matrix[S] Sigma;           # covariance matrix of traits
  
# scaling of variance
Sigma <- diag_matrix(xi_theta) * Sigma_raw* diag_matrix(xi_theta);

 for(i in 1:N){	
	# rescale trait values
	for(s in 1:S){
		theta[i,s] <- theta_raw[i,s] * xi_theta[s];
	}
	
	# loop across items
	for(j in 1:J){
		# IRT: 1. additive combination of trait and difficulty
		#      2. transform to probabilities using Phi (trait: consider reversed items!)
		middle[i,j] <- Phi(theta[i,1]-beta[j,1]);   #Phi_approx()
		if(revItem[j] == 1)
			trait[i,j] <-  Phi(beta[j,2]-theta[i,2]);     # reversed item: beta-theta
		else
			trait[i,j] <-  Phi(theta[i,2]-beta[j,2]);     # standard items: theta-beta
		extreme[i,j] <-  Phi(theta[i,3]-beta[j,3]);
		
		# response probabilities: MPT model for response categories 1, 2, 3, 4, 5
		p_cat[i,j,1] <- (1-middle[i,j])*(1-trait[i,j])*extreme[i,j];
		p_cat[i,j,2] <- (1-middle[i,j])*(1-trait[i,j])*(1-extreme[i,j]);
		p_cat[i,j,3] <-    middle[i,j];
		p_cat[i,j,4] <- (1-middle[i,j])*trait[i,j]*(1-extreme[i,j]);
		p_cat[i,j,5] <- (1-middle[i,j])*trait[i,j]*extreme[i,j];
  }
}

}

model {

xi_theta ~ normal (1, 100);           # scaling parameters
Sigma_raw ~ inv_wishart(df, V);    # raw covariance matrix of traits

# estimate item difficulties freely: 
for(j in 1:J){
	for(s in 1:S){
		beta[j,s] ~ normal(0, 100) T[-10,10];
	}
}

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

################## posterior predictive

generated quantities {
	int<lower=1, upper=5> X_pred[N, J];	    # predicted responses of partipants
	matrix[N,J] dev_obs;
	matrix[N,J] dev_pred;
	real<lower=0> T_obs;
	real<lower=0> T_pred;
	int<lower=0, upper=1> post_p;
	
for(i in 1:N){
	for(j in 1:J){
		X_pred[i,j] <- categorical_rng(p_cat[i,j]);
		
		dev_obs[i,j] <- ((p_cat[i,j, X[i,j]]-1) / p_cat[i,j, X[i,j]])^2;
		dev_pred[i,j] <- ((p_cat[i,j, X_pred[i,j]]-1) / p_cat[i,j, X_pred[i,j]])^2;
	}
}


# posterior predictive checks
T_obs <- sum(dev_obs);
T_pred <- sum(dev_pred);
post_p <- T_pred > T_obs;

}
