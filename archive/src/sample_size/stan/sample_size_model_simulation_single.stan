data {
  int N;
  vector[N] log_depths;
  real beta_0;
  real beta_1;
  real<lower=0> inv_disp;
}
parameters {
}
model {
}
generated quantities {
  vector[2] beta_var = [beta_0, beta_1]';

  vector[N %/% 2] ctrls = to_vector(rep_array(0, N %/% 2));
  vector[N %/% 2] cases = to_vector(rep_array(1, N %/% 2));
  vector[N] case_ctrl = append_row(ctrls, cases);
  vector[N] intercept = to_vector(rep_array(1, N));
  matrix[N, 2] x = append_col(intercept, case_ctrl);

  vector[N] lam = x*beta_var + log_depths;
  array[N] int sim_counts = neg_binomial_2_log_rng(lam, inv(inv_disp));
}
