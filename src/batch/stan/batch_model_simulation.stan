data {
  int N;
  int D;
  vector[N] log_depths;
  vector[D] beta_0;
  vector[D] beta_1;
  vector[D] inv_disp;

  int num_batches;
  array[N] int<lower=1, upper=num_batches> batch_map;
  matrix[num_batches, D] batch_offsets;
  matrix[num_batches, D] batch_disps;
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

  array[D] vector[N] big_alpha;
  array[D] vector[N] big_lam;

  array[D, N] int sim_counts;

  for (d in 1:D) {
    vector[2] beta_var = [beta_0[d], beta_1[d]]';
    vector[N] lam = x*beta_var + log_depths;
    array[N] real alpha;

    for (n in 1:N) {
      lam[n] += batch_offsets[batch_map[n], d];
      alpha[n] = inv_disp[d] + batch_disps[batch_map[n], d];
    }
    sim_counts[d] = neg_binomial_2_log_rng(lam, alpha);
    big_lam[d] = lam;
    big_alpha[d] = to_vector(alpha);
  }
}
