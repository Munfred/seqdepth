// You need to specify the kind of input data, incl. number of observations.
data { 
  int<lower=1> N;  // total number of observations (integer); at least 1
  // for a given size of numebr of cells is the number of sampled depths * 30 
  // (or howver many replicates)
  real validation_error[N];     // VE == scvi validation_error: outcome variable with N elements (real-valued)
  real umis_per_cell[N];  // UPC ==umis per cell: predictor variable with N elements (real-valued)
}

// the parameters to be estimated from the data
parameters { 
  real<lower = 0> intercept;                 // = predicted outcome at breakpoint
  real <upper = 0>  slope_before;              // slope before the breakpoint
  real slope_after;               // slope after the breakpoint
  real<lower = 0, upper = 10000> bp; // the breakpoint: the number of UMIs at which saturation begins
  real<lower = 0> before_variance;          // standard deviation of residuals error before the breakpoint
  real<lower = 0> after_variance;          // standard deviation of residuals  error after the breakpoint

                                  //  (always positive, hence <lower = 0>)
} 

// Functions of estimated parameters.
transformed parameters{
  vector[N] conditional_mean; // the estimated average validation_error for each observation
  real slope_difference;      // the difference between slope_after and slope_before
  real slope_ratio;      // the ratio between slope_after / slope_before

  slope_difference = slope_after - slope_before;  
  slope_ratio = slope_after / slope_before;  

  // conditional_mean depends on whether UPC is before or after bp 
  for (i in 1:N) {
    if (umis_per_cell[i] < bp) {
      conditional_mean[i] = intercept + slope_before * (umis_per_cell[i] - bp);
    } else {
      conditional_mean[i] = intercept + slope_after * (umis_per_cell[i] - bp);
    }
  }
}

// The model itself specifies how the data are expected to have
// been generated and what the prior expectations for the model parameters are.
model {
  // Set priors
  intercept ~ normal(8, 1);  // Average validation_error at breakpoint
  slope_before ~ normal(0, 0.2);  // Slope before breakpoint
  slope_after ~ normal(0, 0.2);   // Slope after breakpoint
  bp ~ normal(8, 3);           // Breakpoint at whih saturation befins, pretty wide, but somewhere in childhood/puberty
  before_variance ~ normal(0, 0.2);        // Residual error before the breakpoint, 
  // stdev in validation error for 30 replicates across datasets ranges between 400 to 8, 300 sounds reasonable here
  after_variance ~ normal(0, 0.2);        // Residual error after the breakpoint
  
  // How the data are expected to have been generated:
  // normal distribution with mu = conditional_mean and 
  // std = error, estimated from data.
  for (i in 1:N) {
    //validation_error[i] ~ normal(conditional_mean[i], error);
    if (umis_per_cell[i] < bp) {
      validation_error[i] ~ normal(conditional_mean[i], before_variance);
    } else {
      validation_error[i] ~ normal(conditional_mean[i], after_variance);
    }
  }
}
/*
generated quantities {
  vector[N] sim_validation_error;               // Simulate new data using estimated parameters.
  vector[N] log_lik;               // Useful for model comparisons; not done here.
  vector[15000] sim_conditional_mean; // Useful for plotting.

  // Compute conditional means for AOAs between 1 and 50.
  for (i in 1:15000) {
    if (i < bp) {
      sim_conditional_mean[i] = intercept + slope_before * (i - bp);
    } else {
      sim_conditional_mean[i] = intercept + slope_after * (i - bp);
    }
  }

  for (i in 1:N) {
    //sim_validation_error[i] = normal_rng(conditional_mean[i], error);
    //log_lik[i] = normal_lpdf(GJT[i] | conditional_mean[i], error);
    if (umis_per_cell[i] < bp) {
      sim_validation_error[i] = normal_rng(conditional_mean[i], before_variance);
      log_lik[i] = normal_lpdf(validation_error[i] | conditional_mean[i], before_variance);
    } else {
      sim_validation_error[i] = normal_rng(conditional_mean[i], after_variance);
      log_lik[i] = normal_lpdf(validation_error[i] | conditional_mean[i], before_variance);

    }

  }
}

*/