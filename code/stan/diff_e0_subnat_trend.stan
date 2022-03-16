data {
  
  int p;      // 29 number of time points (years)  
  int d;      // 43 number of district  

  vector<lower=0>[p] t;      // time index  from 0 to 28
  matrix[p, d] D_e;      // life expectancy at birth over years
 
  
}

parameters {
  
  vector[d] alpha;
  vector[d] beta;
  vector<lower=0>[d] sigma;
  
  real mu_alpha;
  real mu_beta;
  real mu_sigma;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_sigma;
  
}

transformed parameters {
  matrix[p, d] mu;
  
  for (i in 1:d){
       mu[, i] = alpha[i] + beta[i] * t;
  }
  
}

model {
  
  // Priors
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
  sigma ~ normal(mu_sigma, sigma_sigma);
  
  mu_alpha ~ normal(0, 1);
  mu_beta ~ normal(0, 1);
  mu_sigma ~ normal(0, 1);
  sigma_alpha ~ normal(0, 1);
  sigma_beta ~ normal(0, 1);
  sigma_sigma ~ normal(0, 1);

  
  
  // Likelihood
  for (i in 1:d) {
      D_e[, i] ~ normal(mu[, i], sigma[i]);

  }
  
  
}