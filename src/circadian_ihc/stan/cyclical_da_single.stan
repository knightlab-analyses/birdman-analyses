data {
  int<lower=0> N;                      // Number of samples
  real A;
  int<lower=0> p;                      // Number of fixed effects (incl Intercept)
  matrix[N, p] x;                      // Fixed effect design matrix
  int<lower=0> S;                      // Number of subjects
  array[N] real depth;
  array[N] real time;
  array[N] int y;
  array[N] int<lower=1, upper=S> subj_ids;

  real<lower=0> B_p;
  real<lower=0> inv_disp_sd;
  real<lower=0> subj_prior;
}

parameters {
  real<offset=A, multiplier=2> beta_0;
  real<multiplier=B_p> beta_1;
  vector[p-1] beta_x;
  real<lower=0> inv_disp;
  matrix[2, S] subj_re; // Subject intercept, slope
  vector[p] cos_coefs;
  real phase_shift;
  real period_offset_base;
  real period_offset_ihc;
}

transformed parameters {
  vector[p] beta_var = [beta_0, beta_1]';
  vector[N] lam;
  vector[N] amp_terms = x * cos_coefs;

  for (n in 1:N) {
    real period = 24 + period_offset_base + period_offset_ihc*x[n, 2];
    real cos_term = amp_terms[n] * cos(2 * pi() * (time[n] - phase_shift) / period);
    lam[n] = depth[n] + cos_term
      + beta_0 + subj_re[1, subj_ids[n]]  // Global intercept + subj intercept
      + (beta_1 + subj_re[2, subj_ids[n]])*x[n, 2];  // Global IHC slope + subj slope
  }
}

model {
  // Setting priors
  beta_0 ~ normal(A, 2);
  beta_1 ~ normal(0, B_p);

  inv_disp ~ lognormal(0, inv_disp_sd);

  cos_coefs ~ normal(0, 2);
  phase_shift ~ normal(0, 2);

  period_offset_base ~ std_normal();
  period_offset_ihc ~ std_normal();

  for (s in 1:S) {
    subj_re[1, s] ~ normal(0, subj_prior);
    subj_re[2, s] ~ normal(0, subj_prior);
  }

  y ~ neg_binomial_2_log(lam, inv(inv_disp));
}

generated quantities {
  vector[N] y_predict;
  vector[N] log_lhood;

  for (n in 1:N) {
    y_predict[n] = neg_binomial_2_log_rng(lam[n], inv(inv_disp));
    log_lhood[n] = neg_binomial_2_log_lpmf(y[n] | lam[n], inv(inv_disp));
  }
}
