functions{
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
  vector<lower=0>[N_obs] u;     // starting interval
  vector<lower=0>[N_obs] v;     // ending interval
  vector<lower=0>[N_obs] y;     // observed data
  array[N_obs] int group;     //number of different groups
}

parameters {
    vector<lower=0, upper=10000>[3] N;     // effective population sizes
    real<lower=0> T;     //merge times     
}

model {
    for (i in 1:3) {
      N[i] ~ gamma(6.25,0.00125);     //prior such that mean is 5000 and variance is 2000
    }
    T ~ exponential(1.0/25);

    for (i in 1:N_obs) {
      if (group[i]==1){
        y[i] ~ normal(frac_function(N[1],u[i],v[i],0,T) + frac_function(N[3],u[i],v[i],T, 1000000),0.25);
      }
      else if (group[i]==2){
        y[i] ~ normal(frac_function(N[3],u[i],v[i],T, 1000000),0.25);
      }
      else{
        y[i] ~ normal(frac_function(N[2],u[i],v[i],0,T) + frac_function(N[3],u[i],v[i],T, 1000000),0.25);
      }
  }
}