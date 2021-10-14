//A robust Bayesian version of the t   test can accommodate unequal standard deviations as well as outliers, by replacing the normal likelihood with a Student’s t  likelihood. See homework instruction. We will talk more about this in a later week after we discuss regression.
//
  

data {
  int<lower=0> N1;  // number of observations (group 1)
  int<lower=0> N2;  // number of observations (group 2)
  vector[N1] y1;  // response time (group 1);
  vector[N2] y2;  // response time (group 2);
}
parameters {
  real<lower=0> mu_1;  // mean of group 1
  real beta;  // difference in means
  real lsigma_1;  // log of scale parameter for group 1
  real beta_lsigma;  // difference in log standard deviation
  real<lower=1> nu;  // degrees of freedom of student's t distribution
}
transformed parameters {
  real<lower=0> mu_2 = mu_1 + beta; 
  real<lower=0> sigma_1 = exp(lsigma_1);
  real<lower=0> sigma_2 = exp(lsigma_1 + beta_lsigma);
}
model {
  y1 ~ student_t(nu, mu_1, sigma_1);
  y2 ~ student_t(nu, mu_2, sigma_2);
  // prior
  mu_1 ~ normal(0.5, 2.5);
  beta ~ normal(0, 2.5);
  lsigma_1 ~ student_t(4, 0, 1);
  beta_lsigma ~ std_normal();
  nu ~ gamma(2, 0.1);
}
generated quantities {
  real sigma_ratio = sigma_2 / sigma_1;
  real y1rep[N1];
  real y2rep[N2];
  for (i in 1:N1) {
    y1rep[i] = student_t_rng(nu, mu_1, sigma_1);
  }
  for (i in 1:N2) {
    y2rep[i] = student_t_rng(nu, mu_2, sigma_2);
  }
}

