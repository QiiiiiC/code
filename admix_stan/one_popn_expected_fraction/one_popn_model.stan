functions {
  real frac_function(real N, real u, real v, real t1, real t2) {
    real k1 = (2 * N * v + 50) * 50 / pow((N * v + 50), 2) * (exp(-(N * v + 50) / (N * 50) * t2) - exp(-(N * v + 50) / (N * 50) * t1));
    real k2 = -(2 * N * u + 50) * 50 / pow((N * u + 50), 2) * (exp(-(N * u + 50) / (N * 50) * t2) - exp(-(N * u + 50) / (N * 50) * t1));
    real k3 = v / (N * v + 50) * (t2 * exp(-(N * v + 50) / (50 * N) * t2) - t1 * exp(-(N * v + 50) / (50 * N) * t1));
    real k4 = -u / (N * u + 50) * (t2 * exp(-(N * u + 50) / (50 * N) * t2) - t1 * exp(-(N * u + 50) / (50 * N) * t1));
    real r = exp(t1 / N);
    return r * (k1 + k2 + k3 + k4);
  }
  real theo_mean(real N, real u, real v) {
    real o = 100*pow(N,2)*(v-u)*(25*(v+u) + (u*v*N))/(pow(50+u*N,2) * pow(50+v*N,2));
    return o;
  }
  real theo_sigma(real N, real u, real v,real l){
    real a1 = 10*N*(25*v + u*(25+N*v));
    real a2 = (50+N*v)*(50+N*u);
    real a3 = sqrt(2*(v-u)/(l*(100+N*(u+v))));
    return a1/a2*a3;
  }
}
data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] length;
  vector<lower=0>[N_obs] u; // starting interval
  vector<lower=0>[N_obs] v; // ending interval
  vector<lower=0>[N_obs] y;  // observed data
  vector<lower=0>[N_obs] number; //number of observations within bin [u,v]
}
parameters {
  real<lower=0> N;                // effective population sizes

}

model {
    N ~ gamma(6.5,0.00125);
    for (i in 1:N_obs) {
        real mu = theo_mean(N,u[i],v[i]);
        real sigma = theo_sigma(N,u[i],v[i],length[i]);
        y[i] ~ normal(mu, sigma/sqrt(number));
  }
}