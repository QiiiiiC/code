functions{
 
    real p_within(real l, real N, real u){
        real up = (100 * pow(N,2))/(pow((l * N + 50),3));
        real denom = 50 * N / (pow(u * N + 50,2));
        return up/denom;
    }
}
data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] y;     // observed data
}

parameters {
    real N;     
}

model {

    N ~ gamma(6.25,0.00125);     //prior such that mean is 5000 and variance is 2000


    for (i in 1:N_obs) {
        target += log(p_within(y[i],N,0.4));
      }
}