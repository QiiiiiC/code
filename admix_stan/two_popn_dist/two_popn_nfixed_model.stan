// This model assumes identical effective population sizes. 

functions{
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
      if (group[i]==2){
        target += log(p_between(y[i],N[1],T,0.35));
      }
  }
}