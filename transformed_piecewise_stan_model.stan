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
  real intercept;              // = predicted outcome at breakpoint boundary, a constant
  real umi_slope_before;     // umi slope before the breakpoint
  real umi_slope_after;                   // umi slope after the breakpoint
  real cell_slope_before;   // cells slope before the breakpoint
  real cell_slope_after;                  // cells slope after the breakpoint
  real bp;                                // the breakpoint: the number of UMIs at which saturation begins
  real before_variance;        // standard deviation of residuals error before the breakpoint
  real after_variance;         // standard deviation of residuals  error after the breakpoint
  // real bp_intercept;                      // a constant to the breakpoint boundary equation so it doesnt need to touch 0
} 

// Functions of estimated parameters.
transformed parameters{

  // log2 transformed data for centering and scaling
  real log2_validation_error[N];   
  real log2_umis_per_cell[N]; 
  real log2_ncells[N];  

  real mean_log2_validation_error;   
  real mean_log2_umis_per_cell;
  real mean_log2_ncells;  

  real sd_log2_validation_error;   
  real sd_log2_umis_per_cell; 
  real sd_log2_ncells;  

  real standardized_log2_validation_error[N];   
  real standardized_log2_umis_per_cell[N]; 
  real standardized_log2_ncells[N];  

  log2_validation_error = log2(validation_error);
  log2_umis_per_cell = log2(umis_per_cell);
  log2_ncells = log2(ncells);

  mean_log2_validation_error = mean(log2_validation_error);
  mean_log2_umis_per_cell = mean(log2_umis_per_cell);
  mean_log2_ncells = mean(log2_ncells);

  sd_log2_validation_error = sd(log2_validation_error);
  sd_log2_umis_per_cell = sd(log2_umis_per_cell);
  sd_log2_ncells = sd(log2_ncells);

  standardized_log2_validation_error = (log2_validation_error -  mean_log2_validation_error)/sd_log2_validation_error ;
  standardized_log2_umis_per_cell = (log2_umis_per_cell  -  mean_log2_umis_per_cell )/sd_log2_umis_per_cell ;
  standardized_log2_ncells = (log2_ncells -  mean_log2_ncells)/sd_log2_ncells ;


  vector[N] conditional_mean; // the estimated average validation_error for each observation
  // conditional_mean depends on whether umis per cell is before or after bp, the breakpoint

  for (i in 1:N) {
    if (umis_per_cell[i] < bp) { // breakpoint only depends on UMIs

    //if (umis_per_cell[i] + bp_intercept < bp*ncells[i]) { // breakpoint equation
      conditional_mean[i] = intercept + umi_slope_before * (standardized_log2_umis_per_cell[i] - bp) + cell_slope_before * standardized_log2_ncells[i];
    } else {
      conditional_mean[i] = intercept + umi_slope_after * (standardized_log2_umis_per_cell[i] - bp) + cell_slope_after * standardized_log2_ncells[i];
    }
  }
}

// The model itself specifies how the data are expected to have
// been generated and what the prior expectations for the model parameters are.
model {
  // Set priors
  intercept ~ normal(0, 1);  // Average validation_error at breakpoint
  //bp_intercept ~ normal(0, 0.2);  // constant (in umis per cell) to add to the breakpont equation

  cell_slope_before ~ normal(0, 1);  // cell slope before breakpoint
  cell_slope_after ~ normal(0, 1);   // cell slope  breakpoint
  umi_slope_before ~ normal(0, 1);  // umi slope before breakpoint
  umi_slope_after ~ normal(0, 1);   // umi slope  breakpoint
  bp ~ normal(0, 1);           // Breakpoint at which saturation begins pretty wide, but around 1000-8000 umis per cell
  before_variance ~ normal(0, 1);        // Residual error before the breakpoint, 
  // stdev in validation error for 30 replicates across datasets ranges between 400 to 8, 300 sounds reasonable here
  after_variance ~ normal(0, 1);        // Residual error after the breakpoint
  
  // How the data are expected to have been generated:
  // normal distribution with mu = conditional_mean and 
  // std = error, estimated from data.
  for (i in 1:N) {
    //validation_error[i] ~ normal(conditional_mean[i], error);
    // if (umis_per_cell[i] + bp_intercept < bp*ncells[i]) { // breakpoint equation
    if (umis_per_cell[i] < bp) { // breakpoint only depends on UMIs
      standardized_log2_validation_error[i] ~ normal(conditional_mean[i], before_variance);
    } else {
      standardized_log2_validation_error[i] ~ normal(conditional_mean[i], after_variance);
    }
  }
}

generated quantities {
 
  real cell_slope_difference;      // the difference between slope_after and slope_before
  real cell_after_over_before;      // the ratio between slope_after / slope_before
  real cell_before_over_after;      // the ratio between slope_after / slope_before
  real umi_slope_difference;      // the difference between slope_after and slope_before
  real umi_after_over_before;      // the ratio between slope_after / slope_before
  real umi_before_over_after;      // the ratio between slope_after / slope_before
  
  real cell_slope_before_percent; 
  real cell_slope_after_percent; 
  real umi_slope_before_percent; 
  real umi_slope_after_percent; 
  real bp_umis;               // the breakpoint raised to power 2 to get UMIs and not log2(umis)
  

  bp_umis = pow(2,bp*sd_log2_umis_per_cell + mean_log2_umis_per_cell );
  cell_slope_before_percent = 1 - pow(2,cell_slope_before*sd_log2_ncells + mean_log2_ncells);
  cell_slope_after_percent = 1 - pow(2,cell_slope_after*sd_log2_ncells + mean_log2_ncells);
  umi_slope_before_percent = 1 - pow(2,umi_slope_before*sd_log2_umis_per_cell + mean_log2_umis_per_cell );
  umi_slope_after_percent = 1 - pow(2,umi_slope_after*sd_log2_umis_per_cell + mean_log2_umis_per_cell );
  
  cell_slope_difference =  cell_slope_after -  cell_slope_before;  
  cell_after_over_before = cell_slope_after /  cell_slope_before;  
  cell_before_over_after = cell_slope_before / cell_slope_after;  

  umi_slope_difference =  umi_slope_after -  umi_slope_before;  
  umi_after_over_before = umi_slope_after /  umi_slope_before;  
  umi_before_over_after = umi_slope_before / umi_slope_after;  
}
