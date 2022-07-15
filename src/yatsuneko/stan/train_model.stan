data {
  int<lower=1> N;
  int<lower=1> p;
  vector[N] log_depth;
  real A;
  array[N] int y;
  real<lower=0> B_p;
  real<lower=0> disp_scale;
  matrix[N, p] x;
}
parameters {
  real<offset=A, multiplier=2> beta_0;
  real beta_1;
  real<lower=0> phi;
}
transformed parameters {
  vector[p] beta_var = [beta_0, beta_1]';
  vector[N] lam = x*beta_var + log_depth;
}
model {
  beta_0 ~ normal(A, 2);
  beta_1 ~ normal(0, B_p);
  phi ~ lognormal(0, disp_scale);

  y ~ neg_binomial_2_log(lam, inv(phi));
}
generated quantities {
  vector[N] y_predict;
  vector[N] log_lhood;

  for (n in 1:N) {
    y_predict[n] = neg_binomial_2_log_rng(lam[n], inv(phi));
    log_lhood[n] = neg_binomial_2_log_lpmf(y[n] | lam[n], inv(phi));
  }
}
