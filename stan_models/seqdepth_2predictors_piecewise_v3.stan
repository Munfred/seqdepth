// You need to specify the kind of input data, incl. number of observations.
data { 
  int<lower=1> N;  // total number of observations (integer); at least 1
  // for a given size of numebr of cells is the number of sampled depths * 30 
  // (or howver many replicates)
  real validation_error[N];     // scvi validation_error: outcome variable with N elements (real-valued)
  real umis_per_cell[N];  // umis per cell: predictor variable with N elements (real-valued)
  real ncells[N];  // number of cells: predictor variable with N elements (real-valued)

}

// the parameters to be estimated from the data
parameters { 
  real<lower = 0> intercept;              // = predicted outcome at breakpoint boundary, a constant
  real <upper = 0>  umi_slope_before;     // umi slope before the breakpoint
  real umi_slope_after;                   // umi slope after the breakpoint
  real <upper = 0>  cell_slope_before;   // cells slope before the breakpoint
  real cell_slope_after;                  // cells slope after the breakpoint
  real bp;                                // the breakpoint: the number of UMIs at which saturation begins
  real<lower = 0> before_variance;        // standard deviation of residuals error before the breakpoint
  real<lower = 0> after_variance;         // standard deviation of residuals  error after the breakpoint

                                  //  (always positive, hence <lower = 0>)
} 

// Functions of estimated parameters.
transformed parameters{
  vector[N] conditional_mean; // the estimated average validation_error for each observation
  // conditional_mean depends on whether umis per cell is before or after bp, the breakpoint

  vector[N] umis_per_cell_per_cell;

  for (i in 1:N) {
  umis_per_cell_per_cell[i] = umis_per_cell[i] / ncells[i];
  }

  for (i in 1:N) {
    if (umis_per_cell_per_cell[i] < bp) {
      conditional_mean[i] = intercept + umi_slope_before * umis_per_cell[i] + cell_slope_before * ncells[i];
    } else {
      conditional_mean[i] = intercept + umi_slope_after * umis_per_cell[i] + cell_slope_after * ncells[i];
    }
  }
}

// The model itself specifies how the data are expected to have
// been generated and what the prior expectations for the model parameters are.
model {
  // Set priors
  intercept ~ normal(10, 5);  // Average validation_error at breakpoint
  cell_slope_before ~ normal(0, 0.2);  // cell slope before breakpoint
  cell_slope_after ~ normal(0, 0.2);   // cell slope  breakpoint
  umi_slope_before ~ normal(0, 0.2);  // umi slope before breakpoint
  umi_slope_after ~ normal(0, 0.2);   // umi slope  breakpoint
  bp ~ normal(1, 1);           // Breakpoint at which saturation begins pretty wide, but around 1000-8000 umis per cell
  before_variance ~ normal(0, 0.2);        // Residual error before the breakpoint, 
  // stdev in validation error for 30 replicates across datasets ranges between 400 to 8, 300 sounds reasonable here
  after_variance ~ normal(0, 0.2);        // Residual error after the breakpoint
  
  // How the data are expected to have been generated:
  // normal distribution with mu = conditional_mean and 
  // std = error, estimated from data.
  for (i in 1:N) {
    //validation_error[i] ~ normal(conditional_mean[i], error);
    if (umis_per_cell_per_cell[i] < bp) {
      validation_error[i] ~ normal(conditional_mean[i], before_variance);
    } else {
      validation_error[i] ~ normal(conditional_mean[i], after_variance);
    }
  }
}
/*
generated quantities {
  real cell_slope_difference;      // the difference between slope_after and slope_before
  real cell_after_over_before;      // the ratio between slope_after / slope_before
  real cell_before_over_after;      // the ratio between slope_after / slope_before
  real umi_slope_difference;      // the difference between slope_after and slope_before
  real umi_after_over_before;      // the ratio between slope_after / slope_before
  real umi_before_over_after;      // the ratio between slope_after / slope_before
  real bp_umis;               // the breakpoint raised to power 2 to get UMIs and not log2(umis)

  bp_umis = pow(2,bp);
  cell_slope_difference =  cell_slope_after -  cell_slope_before;  
  cell_after_over_before = cell_slope_after /  cell_slope_before;  
  cell_before_over_after = cell_slope_before / cell_slope_after;  

  umi_slope_difference =  umi_slope_after -  umi_slope_before;  
  umi_after_over_before = umi_slope_after /  umi_slope_before;  
  umi_before_over_after = umi_slope_before / umi_slope_after;  
}

*/
































