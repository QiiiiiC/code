functions {
  real frac_function(real N, real u, real v, real t1, real t2) {
    real k1 = (2 * N * v + 50) * 50 / pow((N * v + 50), 2) * (exp(-(N * v + 50) / (N * 50) * t2) - exp(-(N * v + 50) / (N * 50) * t1));
    real k2 = -(2 * N * u + 50) * 50 / pow((N * u + 50), 2) * (exp(-(N * u + 50) / (N * 50) * t2) - exp(-(N * u + 50) / (N * 50) * t1));
    real k3 = v / (N * v + 50) * (t2 * exp(-(N * v + 50) / (50 * N) * t2) - t1 * exp(-(N * v + 50) / (50 * N) * t1));
    real k4 = -u / (N * u + 50) * (t2 * exp(-(N * u + 50) / (50 * N) * t2) - t1 * exp(-(N * u + 50) / (50 * N) * t1));
    real r = exp(t1 / N);
    return r * (k1 + k2 + k3 + k4);
  }
}
data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] u; // starting interval
  vector<lower=0>[N_obs] v; // ending interval
  vector<lower=0>[N_obs] y;  // observed data
}
parameters {
  real<lower=0> N;                // effective population sizes

}

model {
    N ~ gamma(6.25,0.00125);
    for (i in 1:N_obs) {
        real mu = frac_function(N, u[i], v[i], 0, 100000000);
        y[i] ~ normal(mu, 0.25);
  }
}