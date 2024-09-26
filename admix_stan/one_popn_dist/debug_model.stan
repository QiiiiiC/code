functions{
 
    real number_within(real N, real u,real v, real m){
        real up = 50*pow(N,2)*(v-u)*(100+N*(v+u));
        real denom =(pow(v * N + 50,2))*(pow(u * N + 50,2));
        return m*up/denom;
    }
}
data {
  int<lower=0> N_obs;      // number of observations
  array[N_obs] int y;     // observed data
}

parameters {
    real<lower = 0> N;     
}
transformed parameters {
   real lambda = number_within(N,0.35,200000,200);
}

model {

    N ~ gamma(6.25,0.00125);     //prior such that mean is 5000 and variance is 2000

    for (i in 1:N_obs) {
        y[i] ~ poisson(lambda);

      }
}