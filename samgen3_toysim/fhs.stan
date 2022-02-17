/* 
#####	FHS
##############################################
*/

data{
  int<lower=1> N; // number of domains
  row_vector[N] y; // set of N observations
  row_vector<lower=0>[N] sigma_y; // known standard deviations
  int<lower=1> dX; // number of covariates included in matrix x
  matrix[N,dX] x; //Design matrix for dX predictors in model for log-variance
  int<lower=1> P; // number of predictors in model for log-variance
  matrix[N,P] x_var; //Design matrix for P predictors in model for log-variance
  row_vector<lower=0>[N] nResp; // number of respondents
} /* end data block */

transformed data{
   row_vector[N] log_sig2y; // log of sigma_y^2, which we will model under normal() with splines 
   row_vector<lower=0>[N] sigma2_y; // known standard deviations
   vector[P] zros_P;
   vector[dX] zros_X;
   zros_X    = rep_vector(0,dX);

   zros_P    = rep_vector(0,P);


    for( i in 1:N )
    {
     sigma2_y[i]    = (square(sigma_y[i]));
     log_sig2y[i]    = log(sigma2_y[i]);
    } /* end loop i over domains */


} /* end transformed parameters block */


parameters{
  row_vector[N] u; // random effects 
  real mu; // mean (global intercept) of u_i 
  real<lower=0> tau_mu; // precision in prior for mu 

  real<lower=0> lambda_u; 
  real log_f_shape; 
  vector[dX] beta;
  cholesky_factor_corr[dX] L_X; // cholesky of correlation matrix for cols of Beta 
  vector<lower=0>[dX] sigma_X; // vector of sd parameters for P x P, Lambda 
  vector[P] beta_var;
  cholesky_factor_corr[P] L_bvar; // cholesky of correlation matrix for cols of Beta 
  vector<lower=0>[P] sigma_bvar; // vector of sd parameters for P x P, Lambda 
 
  row_vector<lower=0>[N] fitted_vrnc; // fitted values for variances

} /* end parameters block */


transformed parameters{
  row_vector[N] fitted_values; /* fitted values */
  row_vector[N] synth; /* synth part */
  row_vector<lower=0>[N] expvrnc; /* fitted values for variances*/
  row_vector<lower=0>[N] inv_fitted;
  row_vector<lower=0>[N] fitted_sd; /* fitted values for sd*/
  row_vector[N] log_fitted_vars;
  real<lower=0> f_shape; 

  f_shape=exp(log_f_shape);

   { /* block to declare local variables, xb */
    row_vector[N] xb;
    vector[N] xb_logvars;

    xb_logvars              = x_var * beta_var;

    for( i in 1:N )
    {
        xb[i]                 =  x[i] * beta;
        fitted_values[i]      =  xb[i] + u[i];
        synth[i]      =  xb[i] + mu;
        log_fitted_vars[i]    = xb_logvars[i]; //+mu_g;       
        expvrnc[i]         = exp(log_fitted_vars[i]); 
        inv_fitted[i]     =1/fitted_vrnc[i];
        fitted_sd[i]         = sqrt(fitted_vrnc[i]);
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

  // prior for P x 1 regression coefficients, beta_var, in model for log-variance
  { /* local variable block  */
    L_bvar          ~ lkj_corr_cholesky(4);
    sigma_bvar      ~ student_t(3,0,1);
    beta_var        ~ multi_normal_cholesky( zros_P, diag_pre_multiply(sigma_bvar,L_bvar) );
  } /* end local variable block */

   log_f_shape ~ student_t( 3, 0.0, 1.0 );
   fitted_vrnc              ~ inv_gamma(2,  to_vector(expvrnc));

  // observed response likelihood
  y                     ~ normal( fitted_values, fitted_sd);
  sigma2_y              ~gamma(0.5*f_shape *nResp, 0.5*f_shape * nResp .* inv_fitted);

  
} /* end model{} block */


generated quantities{
 
  row_vector[N] y_rep; 
  row_vector[N] sigma2_y_rep; 
  
  for( i in 1:N ) /* by row of y */
    {
      y_rep[i]  = normal_rng( fitted_values[i], fitted_sd[i] );
      sigma2_y_rep[i] = gamma_rng(0.5*f_shape *nResp[i], 0.5*f_shape * nResp[i] * inv_fitted[i]);
    } /* end loop i over areas */
} 
