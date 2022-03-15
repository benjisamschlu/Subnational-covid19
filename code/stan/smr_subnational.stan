data {
  
  int A;      // number of age gps
  int T;      // number of years 
  int R;      // number of districts
  
  int D[R, T];
  matrix[R, T] EXP;
  
}

parameters {
  
  matrix<lower=0>[R, T] theta;
  real<lower=0> mu_theta[T];
  real<lower=0> sigma_theta[T];
  real<lower=0> mu_THETA;
  
}

transformed parameters {
  
  matrix[R, T] log_D_hat;

  
  log_D_hat = log(EXP .* theta);
  
}

model {
  
  // Priors imposing hierarchy
  for (t in 1:T) {
    theta[, t] ~ normal(mu_theta[t], sigma_theta[t]);
    mu_theta[t] ~ normal(mu_THETA, 1);

  }
  sigma_theta ~ normal(1, 1);
  mu_THETA ~ normal(1, 1);
  
  
  // Likelihood
  to_array_1d(D) ~ poisson_log(to_array_1d(log_D_hat'));
  
  
}