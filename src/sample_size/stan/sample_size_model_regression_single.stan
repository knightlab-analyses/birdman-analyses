data {
  int<lower=1> N;
  vector[N] log_depths;
  array[N] int counts;
  matrix[N, 2] x;
}

parameters {
  real<offset=-8, multiplier=2> beta_0;
  real<multiplier=5> beta_1;
  real<lower=0> inv_disp;
}

transformed parameters {
  vector[2] beta_var = [beta_0, beta_1]';
  vector[N] lam = x*beta_var + log_depths;
}

model {
  beta_0 ~ normal(-8, 2);
  beta_1 ~ normal(0, 5);
  inv_disp ~ lognormal(0, 1);

  counts ~ neg_binomial_2_log(lam, inv(inv_disp));
}
