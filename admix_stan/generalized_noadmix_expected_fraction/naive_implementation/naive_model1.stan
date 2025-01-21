functions{
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
        return r*(k1 + k2);
    }
}
data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] y;     // observed data
  vector<lower=0>[N_obs] u;
  vector<lower=0>[N_obs] v;
  vector<lower=0>[N_obs] group;
  vector<lower=0>[2] N_given;
  vector<lower=0>[3] N;

}

parameters {
    vector<lower=0>[2] T;     
}

model {
    //for (i in 1:3){
        //N[i] ~ gamma(6.25,0.00125);     //prior such that mean is 5000 and variance is 2000
    //}

    for (i in 1:2){
        T[i] ~ exponential(1.0/20);
    }



    for (i in 1:N_obs) {
        if (group[i]==1){
            target += log((dist_fun(N[1],0,T[1]+T[2],y[i]) + exp(-(T[1]+T[2])/N[1])*dist_fun(N_given[2],T[1]+T[2],100000000,y[i]))/(dist_denom(N[1],0,T[1]+T[2],u[i],v[i])+exp(-(T[1]+T[2])/N[1])*dist_denom(N_given[2],T[1]+T[2],100000000,u[i],v[i])));
        }
        if (group[i]==2){
            target += log(dist_fun(N_given[2],T[1]+T[2],100000000,y[i])/(dist_denom(N_given[2],T[1]+T[2],100000000,u[i],v[i])));
        }
        if (group[i]==3){
            target += log(dist_fun(N_given[2],T[1]+T[2],100000000,y[i])/(dist_denom(N_given[2],T[1]+T[2],100000000,u[i],v[i])));
        }
        if (group[i]==4){
            target += log((dist_fun(N[2],0,T[1],y[i])+exp(-T[1]/N[2])*dist_fun(N_given[1],T[1],T[2],y[i])+exp(-T[1]/N[2])*exp(-T[2]/N_given[1])*dist_fun(N_given[2],T[1]+T[2],100000000,y[i]))/
            (dist_denom(N[2],0,T[1],u[i],v[i])+exp(-T[1]/N[2])*dist_denom(N_given[1],T[1],T[2],u[i],v[i])+exp(-T[1]/N[2])*exp(-T[2]/N_given[1])*dist_denom(N_given[2],T[1]+T[2],100000000,u[i],v[i]))
            );
        }
        if (group[i]==5){
            target += log((dist_fun(N_given[1],T[1],T[1]+T[2],y[i])+exp(-(T[2]/N_given[1]))*dist_fun(N_given[2],T[1]+T[2],100000000,y[i]))/
            (dist_denom(N_given[1],T[1],T[1]+T[2],u[i],v[i])+exp(-(T[2]/N_given[1]))*dist_denom(N_given[2],T[1]+T[2],100000000,u[i],v[i]))
            );
        }
        if (group[i]==6){
            target += log((dist_fun(N[3],0,T[1],y[i])+exp(-T[1]/N[3])*dist_fun(N_given[1],T[1],T[2],y[i])+exp(-T[1]/N[3])*exp(-T[2]/N_given[1])*dist_fun(N_given[2],T[1]+T[2],100000000,y[i]))/
            (dist_denom(N[3],0,T[1],u[i],v[i])+exp(-T[1]/N[3])*dist_denom(N_given[1],T[1],T[2],u[i],v[i])+exp(-T[1]/N[3])*exp(-T[2]/N_given[1])*dist_denom(N_given[2],T[1]+T[2],100000000,u[i],v[i]))
            );
        }
      }
}