functions{
    real frac_function(real N, real u, real v, real t1, real t2) {
        real k1 = (2 * N * v + 50) * 50 / pow((N * v + 50), 2) * (exp(-(N * v + 50) / (N * 50) * t2) - exp(-(N * v + 50) / (N * 50) * t1));
        real k2 = -(2 * N * u + 50) * 50 / pow((N * u + 50), 2) * (exp(-(N * u + 50) / (N * 50) * t2) - exp(-(N * u + 50) / (N * 50) * t1));
        real k3 = v / (N * v + 50) * (t2 * exp(-(N * v + 50) / (50 * N) * t2) - t1 * exp(-(N * v + 50) / (50 * N) * t1));
        real k4 = -u / (N * u + 50) * (t2 * exp(-(N * u + 50) / (50 * N) * t2) - t1 * exp(-(N * u + 50) / (50 * N) * t1));
        real r = exp(t1 / N);
        return r * (k1 + k2 + k3 + k4);
  }
    real dist_fun(real N, real t1, real t2, real l) {
        real c = -(l / 50 + 1 / N);
        real output = exp(c * t2) * (pow(c * t2 , 2) - 2*c*t2 + 2)/(pow(c,3)) - exp(c * t1) * (pow(c * t1, 2) - 2*c*t1 + 2)/(pow(c,3));
        real r = exp(t1/N)/(2500 * N);
        return r*output;
    }
    real dist_denom(real N, real t1, real t2, real u, real v) {
        real c1 = -(v/50 + 1/N);
        real c2 = -(u/50 + 1/N);
        real k1 = (-1/(50*N)) * ((exp(c1*t2)*(c1*t2-1)/pow(c1,2)) - (exp(c1*t1)*(c1*t1-1)/pow(c1,2)));
        real k2 = (1/(50*N)) * ((exp(c2*t2)*(c2*t2-1)/pow(c2,2)) - (exp(c2*t1)*(c2*t1-1)/pow(c2,2)));
        real r = exp(t1/N);
        return r * (k1 + k2);
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
        target += log((dist_fun(N[1],0,T,y[i]) + dist_fun(N[3],T,1000000,y[i]))/(dist_denom(N[1],0,T,u[i],v[i]) + dist_denom(N[3],T,10000000,u[i],v[i])));
      }
      else if (group[i]==2){
        target += log(dist_fun(N[3],T,100000000,y[i])/dist_denom(N[3],T,10000000,u[i],v[i]));
      }
      else{
        target += log((dist_fun(N[2],0,T,y[i]) + dist_fun(N[3],T,1000000,y[i]))/(dist_denom(N[2],0,T,u[i],v[i]) + dist_denom(N[3],T,10000000,u[i],v[i])));
      }
  }
}