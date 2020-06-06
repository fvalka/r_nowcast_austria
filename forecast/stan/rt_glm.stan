data {
  int<lower=0> N;
  int<lower=0> N0;
  int<lower=0> M;
  int<lower=0> K;
  vector<lower=0>[N] walking;
  vector<lower=0>[N] driving;
  vector<lower=0>[M] Rt;
}

transformed data {
  vector[N] interaction;
  interaction = walking .* driving;
}

parameters {
  vector<lower=0>[K] beta_walking;
  vector<lower=0>[K] beta_driving; 
  vector<lower=0>[K] beta_interaction;
  real alpha;
  real<lower=0> sigma;
  real<lower=0> R0;
  real<lower=0> kappa;
}

model {
  sigma ~ cauchy(0., 2.5);
  alpha ~ normal(2.3, 1.0);
  beta_walking ~ exponential(1);
  beta_driving ~ exponential(1);
  beta_interaction ~ exponential(1);
  kappa ~ normal(0., 0.5);
  R0 ~ normal(2.3, kappa);
  
  for (m in 1:M) {
    real eta = alpha;
    eta += dot_product(beta_walking, walking[(m+N0-K+1):(m+N0)]);
    eta += dot_product(beta_driving, driving[(m+N0-K+1):(m+N0)]);
    eta += dot_product(beta_interaction, interaction[(m+N0-K+1):(m+N0)]);
    eta *= R0;
    Rt[m] ~ normal(eta, sigma);
  }
}

generated quantities {
  vector[N-K+1] Rt_rep;
  for (m in K:N) {
    real eta = alpha;
    eta += dot_product(beta_walking, walking[(m-K+1):(m)]);
    eta += dot_product(beta_driving, driving[(m-K+1):(m)]);
    eta += dot_product(beta_interaction, interaction[(m-K+1):(m)]);
    eta *= R0;
    Rt_rep[m-K+1] = normal_rng(eta, sigma);
  }
}
