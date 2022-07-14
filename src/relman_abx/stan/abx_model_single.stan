data {
  int<lower=1> N;
  real A;
  int<lower=1> p;
  array[N] real depth;
  array[N] int y;
  matrix[N, p] x;
  real<lower=0> B_p;
  real<lower=0> phi_s;

  // Random Effects
  int<lower=1> num_subjs;
  array[N] int<lower=1, upper=num_subjs> subj_map;
  real<lower=0> re_p;
}

parameters {
  real<offset=A, multiplier=2> beta_0;
  vector[p-1] beta_x;
  real<lower=0> reciprocal_phi;
  matrix[num_subjs, p] subj_re;  // Subject effect on intercept + slope of each treatment
}

transformed parameters {
  vector[p] beta_var = append_row(beta_0, beta_x);
  real phi = inv(reciprocal_phi);
  vector[N] lam;

  for (n in 1:N) {
    lam[n] = depth[n];
    for (i in 1:p) {
      // when i = 1 -> Intercept (1)
      // Add subj effect
      lam[n] = lam[n] + (beta_var[i] + subj_re[subj_map[n], i]) * x[n, i];
    }
  }
}

model {
  beta_0 ~ normal(A, 2);
  beta_x ~ normal(0, B_p);

  for (i in 1:num_subjs) {
    for (j in 1:p) {
      subj_re[i, j] ~ normal(0, re_p);
    }
  }
  reciprocal_phi ~ cauchy(0, phi_s);

  y ~ neg_binomial_2_log(lam, phi);
}

generated quantities {
  vector[N] y_predict;
  vector[N] log_lhood;

  for (n in 1:N) {
    y_predict[n] = neg_binomial_2_log_rng(lam[n], phi);
    log_lhood[n] = neg_binomial_2_log_lpmf(y[n] | lam[n], phi);
  }
}
