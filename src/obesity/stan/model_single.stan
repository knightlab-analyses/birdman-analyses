data {
  int<lower=1> N;
  real A;
  int<lower=1> p;
  array[N] real depth;
  array[N] int y;
  matrix[N, p] x;
  real<lower=0> B_p;
  real<lower=0> disp_scale;

  // Random Effects
  int<lower=1> num_studies;
  array[N] int<lower=1, upper=num_studies> study_map;
  real<lower=0> re_p;
}

parameters {
  real<offset=A, multiplier=2> beta_0;
  vector[p-1] beta_x;
  real<lower=0> base_phi;
  vector<lower=0>[num_studies] study_disp;
  vector[num_studies] study_re;
}

transformed parameters {
  vector[p] beta_var = append_row(beta_0, beta_x);
  vector[N] transformed_phi;
  vector[N] lam = x * beta_var;

  for (n in 1:N) {
    lam[n] += depth[n] + study_re[study_map[n]];
    transformed_phi[n] = base_phi + study_disp[study_map[n]];
  }
}

model {
  beta_0 ~ normal(A, 2);
  beta_x ~ normal(0, B_p);
  base_phi ~ lognormal(log(10), disp_scale);

  for (i in 1:num_studies) {
    study_disp[i] ~ lognormal(0, 1);
    study_re[i] ~ normal(0, re_p);
  }

  y ~ neg_binomial_2_log(lam, inv(transformed_phi));
}

generated quantities {
  vector[N] y_predict;
  vector[N] log_lhood;

  for (n in 1:N) {
    y_predict[n] = neg_binomial_2_log_rng(lam[n], inv(transformed_phi[n]));
    log_lhood[n] = neg_binomial_2_log_lpmf(y[n] | lam[n], inv(transformed_phi[n]));
  }
}
