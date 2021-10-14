data {
  int<lower=0> N1;  // number of observations (group 1)
  int<lower=0> N2;  // number of observations (group 2)
  vector<lower=0>[N1] y1;  // response  (group 1);
  vector<lower=0>[N2] y2;  // response  (group 2);
}
parameters {
  real mu_1;  // mean of group 1 for the Gaussian component
  real beta_mu;  // difference in means for the Gaussian component
  real<lower=0> sigma;  // pooled standard deviation
  real lndt_1;  // log of non-decision time for group 1
  real beta_lndt;  // difference in ndt
}
transformed parameters {
  real mu_2 = mu_1 + beta_mu; 
  real<lower=0> ndt_1 = exp(lndt_1);
  real<lower=0> ndt_2 = exp(lndt_1 + beta_lndt);
}
model {
  target += lognormal_lpdf(y1 - ndt_1 | mu_1, sigma);
  target += lognormal_lpdf(y2 - ndt_2 | mu_2, sigma);
  target += std_normal_lpdf(mu_1);
  target += std_normal_lpdf(beta_mu);
  target += student_t_lpdf(sigma | 4, 0, 1) - 
            student_t_lccdf(0 | 4, 0, 1);
  target += std_normal_lpdf(lndt_1);
  target += std_normal_lpdf(beta_lndt);
}
generated quantities {
  real<lower=0> y1rep[N1];
  real<lower=0> y2rep[N2];
  real MG = ndt_1 + exp(mu_1 + 0.5*sigma^2);
  real SOL = ndt_2 + exp(mu_2 + 0.5*sigma^2);

  for (i in 1:N1) {
    y1rep[i] = lognormal_rng(mu_1, sigma) + ndt_1;
  }
  for (i in 1:N2) {
    y2rep[i] = lognormal_rng(mu_2, sigma) + ndt_2;
  }
  
}


