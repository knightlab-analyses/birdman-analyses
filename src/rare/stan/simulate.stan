data {
  int N;
  int D;
  int depth;
  vector[D] beta_0;
  vector[D] beta_1;
}
parameters {
}
model {
}
generated quantities {
  vector[N %/% 2] ctrls = to_vector(rep_array(0, N %/% 2));
  vector[N %/% 2] cases = to_vector(rep_array(1, N %/% 2));
  vector[N] case_ctrl = append_row(ctrls, cases);
  vector[N] intercept = to_vector(rep_array(1, N));
  matrix[N, 2] x = append_col(intercept, case_ctrl);
  matrix[D, 2] beta_var = append_col(beta_0, beta_1);

  array[N, D] int sim_counts;

  for (n in 1:N) {
    vector[D] gamma_var = to_vector(x[n] * beta_var');
    sim_counts[n] = multinomial_logit_rng(gamma_var, depth);
  }
}
