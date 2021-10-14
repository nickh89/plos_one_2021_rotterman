data {
  int<lower=1> n;
  int<lower=1> n_Muscle;
  int<lower=1, upper=n_Muscle> Muscle[n];
  real Dis_2_Soma[n];
}
parameters {
  real overall_mean;
  vector[n_Muscle] Muscle_zoffset;
  real<lower=0> Dis_2_Soma_sd;
  real<lower=0> Muscle_mean_sd;
}
transformed parameters {
  vector[n_Muscle] Muscle_mean;
  Muscle_mean = overall_mean + Muscle_zoffset * Muscle_mean_sd;
}
model {
  Dis_2_Soma_sd ~ cauchy(0, 1);         // => half-cauchy(0, 1)
  Muscle_mean_sd ~ cauchy(0, 1);   // => half-cauchy(0, 1)
  overall_mean ~ normal(0, 5);
  Muscle_zoffset ~ normal(0, 1);   // => condition_mean ~ normal(overall_mean, condition_mean_sd)
  for (i in 1:n) {
    Dis_2_Soma[i] ~ normal(Muscle_mean[Muscle[i]], Dis_2_Soma_sd);
  }
}
