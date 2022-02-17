/* 
#####	A.	Model 1 (FH)
##############################################
*/

data{
  int<lower=1> N; // number of domains
  row_vector[N] y; // set of N observations
  row_vector<lower=0>[N] sigma_y; // known standard deviations
  int<lower=1> dX; // number of covariates included in matrix x
  matrix[N,dX] x; //Design matrix for P predictors in model for log-variance
} /* end data block */

transformed data{
   vector[dX] zros_X;
   zros_X    = rep_vector(0,dX);
} /* end transformed parameters block */



parameters{
  row_vector[N] u; /* random effects */
  real mu; /* mean (global intercept) of u_i */
  real<lower=0> tau_mu; /* precision in prior for mu */
  real<lower=0> lambda_u; 
  vector[dX] beta;
  cholesky_factor_corr[dX] L_X; /* cholesky of correlation matrix for cols of Beta */
  vector<lower=0>[dX] sigma_X; /* vector of sd parameters for P x P, Lambda */
} /* end parameters block */


transformed parameters{
  row_vector[N] fitted_values; /* fitted values */
  row_vector[N] synth; /* synth part */
   { /* block to declare local variables, xb */
    row_vector[N] xb;
    for( i in 1:N )
    {
        xb[i]                 =  x[i] * beta;
        fitted_values[i]      =  xb[i] + u[i];
        synth[i]      =  xb[i] + mu;
    } /* end loop i over domains */
   } /* end local block to declare xb */
 } /* end transformed parameters block */

model{
  // prior N random effects, u 
  { /* local block for parameters */
    lambda_u    ~ gamma( 1.0, 1.0 );
    tau_mu      ~ gamma( 1.0, 1.0 );
    mu   ~ normal(0,inv(sqrt(tau_mu))); 
    u    ~ normal(mu,inv(sqrt(lambda_u))); /* N x 1 vectorized */
  } /* end local block for parameters */
    
  // prior for dX x 1 regression coefficients, beta, in the fixed effects part of the model
  { /* local variable block  */
    L_X          ~ lkj_corr_cholesky(4);
    sigma_X      ~ student_t(3,0,1);
    beta        ~ multi_normal_cholesky( zros_X, diag_pre_multiply(sigma_X,L_X) );
  } /* end local variable block */


  
  // observed response likelihood
  y         ~ normal( fitted_values, sigma_y);
  
} /* end model{} block */

generated quantities{
 
  row_vector[N] y_rep; 
  
  for( i in 1:N ) /* by row of y */
    {
      y_rep[i]  = normal_rng( fitted_values[i], sigma_y[i] );
    } /* end loop i over areas */
} 
