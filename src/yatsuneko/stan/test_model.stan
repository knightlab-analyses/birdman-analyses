data {
  int<lower=1> D;                           // Number of microbes
  int<lower=1> N;                           // Number of samples
  int<lower=1> draws;                       // Number of draws

  array[N] real log_depths;                 // Log sequencing depths
  array[N, D] int<lower=0> y;               // Count data
  matrix[N, 2] x;                           // Design matrix (all ones)
  array[draws] matrix[2, D] post_beta_var;  // Posterior draws
  array[draws] vector[D] post_phi;
}
parameters {

}
model {

}
generated quantities {
  array[draws] matrix[N, 2] all_log_lhood;  // Log-likelihood of each class

  for (i in 1:draws) {
    matrix[N, 2] log_lhood = rep_matrix(rep_row_vector(0, 2), N);
    matrix[2, D] beta_var = post_beta_var[i];
    matrix[N, D] lam1 = col(x, 1) * row(beta_var, 1);  // Intercept only
    matrix[N, D] lam2 = x * beta_var;  // Intercept + beta

    vector[D] phi = post_phi[i];

    for (n in 1:N) {
      for (d in 1:D) {
        log_lhood[n, 1] += neg_binomial_2_log_lpmf(y[n, d] | lam1[n, d] + log_depths[n], inv(phi[d]));
        log_lhood[n, 2] += neg_binomial_2_log_lpmf(y[n, d] | lam2[n, d] + log_depths[n], inv(phi[d]));
      }
    }
    all_log_lhood[i] = log_lhood;
  }
}
