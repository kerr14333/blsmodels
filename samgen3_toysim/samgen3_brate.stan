/* 
#####	CFHSG used to be called Cfhsindvb_reff.stan, cleaned
##############################################
*/

data{
  int<lower=1> N;                     //  number of domains
  int<lower=1> K;                     //  number of clusters
  int<lower=1> Nst;                   //  number of States 
  row_vector[N] y;                    //  set of N observations
  row_vector<lower=0>[N] sigma_y;     //  observed standard deviations
  int<lower=1> dX;                    //  number of covariates included in matrix x
  matrix[N,dX] x;                     //  Design matrix for dX predictors 
  int<lower=1> P;                     //  number of predictors in model for log-variance
  vector[P] x_var[N];                 //  Design matrix for P predictors in model for log-variance
  row_vector<lower=0>[N] nResp;       //  number of respondents in each domain (standardized)
  int<lower=1> Bd;                    //  number of predictors for clustering structure
  row_vector[Bd] B[N];                //  clustering structure predictors
  int<lower=1,upper=53> fipsindex[N]; //  ordered index of state fips
} /* end data block */


transformed data{
   row_vector<lower=0>[N] sigma2_y;   // observed variances
   vector[P] zros_P;
   vector<lower=0>[K] ones_K;
   vector[dX] zros_X;
   row_vector[Bd] zros_Bd;

   zros_Bd   = rep_row_vector(0,Bd);  // vector of zeros for sampling 1xBd  
   zros_X    = rep_vector(0,dX);
   ones_K    = rep_vector(1,K);       // dirichlet prior on alpha has equal shapes
   zros_P    = rep_vector(0,P);


   for( i in 1:N )
   {
     sigma2_y[i]     = (square(sigma_y[i]));
   } /* end loop i over domains */


} /* end transformed parameters block */




parameters{
  row_vector[N] u;                    /* random effects */
  row_vector[K] mu;                   /* mean (global intercept) of u_i */
  simplex[K] theta;                   /* mixture probabilities for model for point estimate */
  row_vector<lower=0>[K] b_rate;      /* bias correction for variances */
  
  real<lower=0> alpha;                /* DP concentration parameter on mixture model for point estimate */


  real<lower=0> tau_mu;               /* precision in prior for mu */
  real<lower=0> lambda_u;             /* precision on u */

  row_vector[Nst] u_st;               /* State random effects */
  real<lower=0> lambda_u_st;          /* precision on u_st */

  real<lower=0> f_shape;              /* parameter for shape/scale of variances distribution */
  vector[dX] beta;                    /* linear model coefficients */
  cholesky_factor_corr[dX] L_X;       /* cholesky of correlation matrix for cols of Beta */
  vector<lower=0>[dX] sigma_X;        /* vector of sd parameters for P x P, Lambda */
 
  vector[P] beta_var[K];              /* model coefficients for model for variance */
  cholesky_factor_corr[P] L_bvar;     /* cholesky of correlation matrix for cols of Beta */
  vector<lower=0>[P] sigma_bvar;      /* vector of sd parameters for P x P, Lambda */


  row_vector[Bd] mu_B[K];             /* means of B */
  cholesky_factor_corr[Bd] L_muB ;    /* cholesky of correlation matrix for mu_B   */
  row_vector<lower=0>[Bd] sigma_muB ; /* vector of sd parameters for Bd x Bd matrix */
  row_vector<lower=0>[Bd] sigma_B ;   /* vector of sd parameters for Bd x Bd matrix */
 
  row_vector<lower=0>[N] fitted_vrnc; /* fitted values for variances*/


} /* end parameters block */


transformed parameters{
  row_vector[N] fitted_values;        /* fitted values */
  row_vector<lower=0>[N] inv_fitted;  /* inverse of fitted variance */
  row_vector<lower=0>[N] fitted_sd;   /* fitted values for sd*/
  real<lower=0> tau_u2 ;              /* variances of random effects */
  real<lower=0> tau_u2_st ;           /* variances of State random effects */
  row_vector[N] xb;                   /* fitted Xb values */
 
  tau_u2=1/lambda_u;          
  tau_u2_st=1/lambda_u_st;          

 { /* block to declare local variables, xb */
  for( i in 1:N )
  {
      xb[i]                 =  x[i] * beta;
      fitted_values[i]      =  xb[i] + u_st[fipsindex[i]] + u[i];
      inv_fitted[i]         = 1/(fitted_vrnc[i]); 
      fitted_sd[i]          = sqrt(fitted_vrnc[i]);
  } /* end loop i over domains */
 } /* end local block to declare xb */


 } /* end transformed parameters block */

model{

   alpha 	~ gamma( 1.0, 1.0 );              /* DP concentration parameter for point estimate */
    // priors for mean and covariance cluster locations
   theta 	~ dirichlet( alpha/K * ones_K );  /* instantiate a truncated DP prior for point */

 { /* local block for priors for K sets of  cluster centers and variances*/
    tau_mu            ~ gamma( 1.0, 1.0 );
    lambda_u          ~ gamma( 1.0, 1.0); 
    lambda_u_st       ~ gamma( 1.0, 1.0); 
    L_bvar            ~ lkj_corr_cholesky(4);
    sigma_bvar        ~ student_t(3,0,1);
    L_muB             ~ lkj_corr_cholesky(4);
    sigma_muB         ~ student_t(3,0,1);
    sigma_B           ~ student_t(3,0,1);
    
    u_st              ~ normal(0,sqrt(tau_u2_st));
 
    for( k in 1:K )
    {
      mu[k]             ~ normal(0,inv(sqrt(tau_mu)));
      b_rate[k]         ~ gamma( 10.0, 10.0 );
      mu_B[k]           ~ multi_normal_cholesky( zros_Bd, diag_pre_multiply(sigma_muB,L_muB) );
      beta_var[k]       ~ multi_normal_cholesky( zros_P, diag_pre_multiply(sigma_bvar,L_bvar) ); /* P x 1 */
    }
 } /* end local block for cluster parameters */


{
    /* mixture prior for u[i], i = 1 ,.., N */
    for( i in 1:N ) /* by row of N x T, u and beta_var */
    {
      vector[K] ps;
      for( k in 1:K )
      {               // Random effects, u

        ps[k]       	= log(theta[k]) + normal_lpdf(u[i]| mu[k], sqrt(tau_u2)) + 
                      // True variances, fitted_vrnc
                      inv_gamma_lpdf(fitted_vrnc[i] | 2, exp(dot_product(x_var[i],beta_var[k]))) +
                      // Observed variances, sigma2_y
                      gamma_lpdf(sigma2_y[i] | 
                        0.5*f_shape *nResp[i], 0.5*f_shape * nResp[i] * inv_fitted[i]/b_rate[k]) +
                      normal_lpdf(B[i]|mu_B[k],sigma_B);

      } /* end loop k over clusters / mixture components */
      target += log_sum_exp(ps);
    } /* end loop i over case observations */
  } /* end local block in mixture prior for u */


  // prior for dX x 1 regression coefficients, beta, in the fixed effects part of the model
  { /* local variable block  */
    L_X          	~ lkj_corr_cholesky(4);
    sigma_X      	~ student_t(3,0,1);
    beta        	~ multi_normal_cholesky( zros_X, diag_pre_multiply(sigma_X,L_X) );
  } /* end local variable block */


   f_shape      	~ student_t( 3, 0, 1.0 );

    // observed response likelihood
   y                    ~ normal( fitted_values, fitted_sd);
  
} /* end model{} block */

generated quantities{
 
  row_vector[N] y_rep; 

  for( i in 1:N ) /* by row of y */
    {
      y_rep[i]	= normal_rng( fitted_values[i], fitted_sd[i] );
    }
} 
