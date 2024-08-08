functions {
  vector append_value_to_vector(vector v, real new_value) {
    int n = num_elements(v);
    vector[n + 1] new_vector;
    for (i in 1:n) {
      new_vector[i] = v[i];
    }
    new_vector[n + 1] = new_value;
    return new_vector;
  }

  vector prepend_value_to_vector(vector v, real new_value) {
    int n = num_elements(v);
    vector[n + 1] new_vector;
    new_vector[1] = new_value; // Prepend the new value
    for (i in 1:n) {
      new_vector[i + 1] = v[i]; // Copy original elements
    }
    return new_vector;
  }

  vector indicator_vector(int d, int i) {
    vector[d] v = rep_vector(0.0, d); 
    v[i] = 1.0; 
    return v;
  }


  real frac_function(real N, real u, real v, real t1, real t2) {
    real k1 = (2 * N * v + 50) * 50 / pow((N * v + 50), 2) * (exp(-(N * v + 50) / (N * 50) * t2) - exp(-(N * v + 50) / (N * 50) * t1));
    real k2 = -(2 * N * u + 50) * 50 / pow((N * u + 50), 2) * (exp(-(N * u + 50) / (N * 50) * t2) - exp(-(N * u + 50) / (N * 50) * t1));
    real k3 = v / (N * v + 50) * (t2 * exp(-(N * v + 50) / (50 * N) * t2) - t1 * exp(-(N * v + 50) / (50 * N) * t1));
    real k4 = -u / (N * u + 50) * (t2 * exp(-(N * u + 50) / (50 * N) * t2) - t1 * exp(-(N * u + 50) / (50 * N) * t1));
    real r = exp(t1 / N);
    return r * (k1 + k2 + k3 + k4);
  }

  real expected_ratio(vector N, vector T, matrix[] A, real u, real v,int j, int k) {
    int N_events = num_elements(T);     // Number of events
    int N_popn = num_elements(N);
    real out = 0;     // Initialize output matrix
    real weight = 1;     // Initialize weight matrix
    vector[N_popn] dist_j = indicator_vector(N_popn, j);     // Initial probability distribution of population j
    vector[N_popn] dist_k = indicator_vector(N_popn, k);     // Initial probability distribution of population j

    vector[N_events+1] T_new = append_value_to_vector(T,100000000);     // Append large value to T
    vector[N_events+2] T_final = prepend_value_to_vector(T_new,0);     //get the final T as [0,T,100000000]

    for (i in 1:(N_events+1)) {
      dist_j = A[i]*dist_j;
      dist_k = A[i]*dist_k;
      vector[N_popn] dist_together = dist_j .* dist_k;
      for (l in 1:N_popn){
        out += weight * frac_function(N[l],u,v,T[i],T[i+1]) * dist_together[l];
      }
      vector[N_popn] ff = rep_vector(0.0, N_popn);
      for (l in 1:N_popn) {
        ff[l] = 1 - exp(-(T[i+1] - T[i]) / N[l]);
      }
      real gj = dot_product(dist_together, ff);
      weight -= gj * weight;
    }
    return out;
  }
}


data {
  int<lower=0> N_obs;      // number of observations
  int<lower=0> N_starting;     //number of starting populations
  int<lower=0> N_popn;     //number of total populations
  int<lower=0> N_events;     //number of merge events
  int<lower=0> N_bins;     //number of different bins considered
  vector<lower=0>[N_bins] L;     //bins list
  array[N_events+1] matrix[N_popn, N_popn] A;     //membership matrices     
  vector<lower=0>[N_obs] u;     // starting interval
  vector<lower=0>[N_obs] v;     // ending interval
  vector<lower=0>[N_obs] y;     // observed data
  array[N_obs,2] int group;     //number of different groups
}

parameters {
    vector<lower=0>[N_popn] N;     // effective population sizes
    vector<lower=0>[N_events] T;     //merge times     
}

model {
    for (i in 1:N_popn) {
      N[i] ~ gamma(6.25,0.00125);     //prior such that mean is 5000 and variance is 2000
    }
    for (i in 1:N_events) {
      T[i] ~ exponential(1.0/25);     //prior such that mean is 25
    }

    for (i in 1:N_obs) {
      real true_mean = expected_ratio(N, T, A, u[i], v[i], group[i,1], group[i,2]);
      y[i] ~ normal(true_mean,0.25);
  }
}