data {
  int N;
  array[N] int y;
  int depth;
  vector[N] x;
  real beta_0;
}
parameters {
  real beta_1;
}
transformed parameters {
  vector[N] lam;
  lam = beta_0 + x*beta_1 + log(depth);
}
model {
  beta_1 ~ normal(0, 3);
  y ~ poisson_log(lam);
}
