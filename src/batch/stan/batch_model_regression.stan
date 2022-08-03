data {
  int N;
  vector[N] log_depths;
  array[N] int y;
  matrix[N, 2] x;

  int num_batches;
  array[N] int<lower=1, upper=num_batches> batch_map;
}
parameters {
  real beta_0;
  real beta_1;
  real<lower=0> inv_disp;

  array[num_batches] real<lower=0> batch_disps;
  array[num_batches] real batch_offsets;
}
transformed parameters {
  vector[2] beta_var = [beta_0, beta_1]';
  vector[N] lam = x * beta_var + log_depths;
  vector[N] alpha;

  for (n in 1:N) {
    lam[n] += batch_offsets[batch_map[n]];
    alpha[n] = inv_disp + batch_disps[batch_map[n]];
  }
}
model {
  beta_0 ~ normal(-6, 3);
  beta_1 ~ normal(0, 5);

  inv_disp ~ lognormal(0, 0.5);
  batch_disps ~ lognormal(0, 0.5);

  y ~ neg_binomial_2_log(lam, inv(alpha));
}
