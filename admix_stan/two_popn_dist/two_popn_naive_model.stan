// This model assumes identical effective population sizes. 

functions{
    real frac_function(real N, real u, real v, real t1, real t2) {
        real k1 = (2 * N * v + 50) * 50 / pow((N * v + 50), 2) * (exp(-(N * v + 50) / (N * 50) * t2) - exp(-(N * v + 50) / (N * 50) * t1));
        real k2 = -(2 * N * u + 50) * 50 / pow((N * u + 50), 2) * (exp(-(N * u + 50) / (N * 50) * t2) - exp(-(N * u + 50) / (N * 50) * t1));
        real k3 = v / (N * v + 50) * (t2 * exp(-(N * v + 50) / (50 * N) * t2) - t1 * exp(-(N * v + 50) / (50 * N) * t1));
        real k4 = -u / (N * u + 50) * (t2 * exp(-(N * u + 50) / (50 * N) * t2) - t1 * exp(-(N * u + 50) / (50 * N) * t1));
        real r = exp(t1 / N);
        return r * (k1 + k2 + k3 + k4);
    }
    real p_within(real l, real N, real u){
        real up = (100 * pow(N,2))/(pow((l * N + 50),3));
        real denom = 50 * N / (pow(u * N + 50,2));
        return up/denom;

    }
    real p_between(real l, real N, real T, real u){
        real c = (l*N + 50)/(50*N);
        real k1 = exp(-c*T) * (pow(c*T,2) + 2*c*T + 2) / pow(c,3);
        real up = l / (2500*N) * k1;
        real denom_c = (N*u + 50)/(50*N);
        real denom = exp(-denom_c*T)*(denom_c*T + 1)/(pow(denom_c,2) * 50 * N);
        return up/denom;
    }
}

data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] u;     // starting interval
  vector<lower=0>[N_obs] v;     // ending interval
  vector<lower=0>[N_obs] y;     // observed data
  array[N_obs] int group;     //number of different groups
  vector[3] N;
}

parameters {
    //vector<lower=0, upper=10000>[3] N;     // effective population sizes
    real<lower=0> T;     //merge times     
}

model {
    //for (i in 1:3) {
      //N[i] ~ gamma(6.25,0.00125);     //prior such that mean is 5000 and variance is 2000
    //}
    T ~ exponential(1.0/20);

    for (i in 1:N_obs) {
      if (group[i]==1){
        target += log(p_within(y[i],N[1],0.35));
      }
      else if (group[i]==2){
        target += log(p_between(y[i],N[1],T,0.35));
      }
      else{
        target += log(p_within(y[i],N[1],0.35));
      }
  }
}