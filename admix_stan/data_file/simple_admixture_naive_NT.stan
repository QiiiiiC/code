functions {
    real p_L_divide_l(real N, real t1, real t2, real l) {
        real c = -(l / 50 + 1 / N);
        real output = exp(c * t2) * (pow(c * t2 , 2) - 2*c*t2 + 2)/(pow(c,3)) - exp(c * t1) * (pow(c * t1, 2) - 2*c*t1 + 2)/(pow(c,3));
        real r = exp(t1/N)/(2500 * N);
        return r*output;
    }
    real int_p_L_divide_l(real N, real t1, real t2, real u, real v) {
        real c1 = -(v/50 + 1/N);
        real c2 = -(u/50 + 1/N);
        real k1 = (-1/(50*N)) * ((exp(c1*t2)*(c1*t2-1)/pow(c1,2)) - (exp(c1*t1)*(c1*t1-1)/pow(c1,2)));
        real k2 = (1/(50*N)) * ((exp(c2*t2)*(c2*t2-1)/pow(c2,2)) - (exp(c2*t1)*(c2*t1-1)/pow(c2,2)));
        real r = exp(t1/N);
        return r*(k1 + k2);
    }
    real int_p_L_divide_l_limit(real N, real t, real u, real v) {
        real cu = -(1/N + u/50);
        real cv = -(1/N + v/50);
        return -1/(50*N)*exp(t/N) * (exp(cu*t)*(cu*t-1)/pow(cu,2) - exp(cv*t)*(cv*t-1)/pow(cv,2));
    }
    real p_L_in_R(real N, real t1, real t2, real u, real v){
      real cu = -(1/N + u/50);
      real cv = -(1/N + v/50);
      real a1 = exp(cu*t2)*(cu*t2-1) - exp(cu*t1)*(cu*t1-1)/(pow(cu,2));
      real a2 = -exp(cv*t2)*(cv*t2-1) - exp(cv*t1)*(cv*t1-1)/(pow(cv,2));
      return a1;
    }
}

data {
  int<lower=0> N_obs;      // number of observations
  vector<lower=0>[N_obs] u;     // starting interval
  vector<lower=0>[N_obs] v;     // ending interval
  array[N_obs] int number;     // observed data
  array[N_obs,2] int group;     //number of different groups
  real length;
  real<lower=0> N;
  real<lower=1,upper = 1000>T2;         //merge times
}

parameters {
    real<lower=1, upper = T2> T1;     //admixture times 
    real<lower=0,upper = 1> fraction;

}


model {
    T1 ~ uniform(1,200);
    fraction ~ uniform(0,1);

    for (i in 1:N_obs) {
      // This is the sharing number between population A and ADMIX
      real lambda;
      lambda = 1;
      if (group[i][1]==0 && group[i][2]==1){
        lambda = (fraction * int_p_L_divide_l(N,T1,T2,u[i],v[i]) + int_p_L_divide_l_limit(N, T2, u[i], v[i]))*length;
        target += poisson_lpmf(number[i] | lambda);
      }
      // This is the sharing number between population Admix and Admix
      else if (group[i][1]==1 && group[i][2]==1){
        lambda = (int_p_L_divide_l(N,0,T1,u[i],v[i]) + pow(fraction,2) * int_p_L_divide_l(N,T1,T2,u[i],v[i]) + pow(1-fraction,2)*int_p_L_divide_l(N,T1,T2,u[i],v[i]) + int_p_L_divide_l_limit(N, T2, u[i], v[i]))*length;
        target += poisson_lpmf(number[i] | lambda);
      }
      // This is the sharing number between population Admix and C
      else if (group[i][1]==1 && group[i][2]==2){
        lambda = ((1-fraction) * int_p_L_divide_l(N,T1,T2,u[i],v[i]) + int_p_L_divide_l_limit(N, T2, u[i], v[i]))*length;
        target += poisson_lpmf(number[i] | lambda);
      }
      
    }

}