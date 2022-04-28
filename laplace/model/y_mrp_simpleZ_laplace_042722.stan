functions{
  
  matrix K_functor(array[] int weight_index, array[] int unit, vector sigma2_lambda) {
    int n                  = num_elements(unit);
    matrix[n, n] K         = diag_matrix( sigma2_lambda[weight_index[unit]] );
    return K;
  } /* end latent response covariance function, K_functor() */
  
  real L_functor(
  vector theta, // latent Gaussian variable
  vector mu_vec, // parameter-dependent variables for the likelihood (eta)
  vector dat, // real data (mean offset)
  array[] int dat_int) { // interger data (observations)
     int n_obs                      = num_elements(dat_int) / 2; // integer valued
     array[n_obs] int y_obs         = dat_int[1:n_obs];
     array[n_obs] int ind_obs       = dat_int[(n_obs+1):2*n_obs];
     array[n_obs] real alpha_obs    = to_array_1d(mu_vec[ind_obs] + theta[ind_obs]);
     return poisson_log_lpmf(y_obs | alpha_obs);
  }  /* end marginalized likelihood function, L_functor() */
  
} /* end functions{} block */


// The input data is a vector 'y' of length 'N'.
data{
    int<lower=1> n; // number of respondents.  A respondent is a report_with record, split further if necessary to
                    // fit within detailed industry - series code (SS + 6 digit NAICS) and MSA
                    //detailed industry: ESS drop down is the (SS)+(2 Digit Naics), the sub industry is 
                    // approx 4 digit-naics, so in case of '605411' its super sector 60 
                    // and 2 digit naics 54, and its 4 digit naics is 541
    int<lower=1> T; // number of mohths to jointly model y
	  int<lower=1> K; // number of linear predictors in marginal model for y
	  int<lower=1> K_r; // number of levels of main and interaction random effects
	  int<lower=0> n_obs; 		// Number of non-missing values
    array[n_obs] int<lower=1> ind_obs;     	// Vector of non-missing indices in n*T x 1 vector, y
    int<lower=1> n_w; // number of unique weights that will be used to define post-strata
    int<lower=1> n_1; // number of main random effects in K_r = n_1 + n_2 + n_3
    int<lower=1> n_2; // number of second order interaction random effects in K_r
    int<lower=1> n_3; // number of third order interaction random effects in K_r
    array[n_2,2] int<lower=1> second_order; // Each row denotes two main effects (\in 1,..n_1)
                                        // attached to each second order random effect
    array[n_3,3] int<lower=1> third_order; // Each row denotes three main effects (\in 1,..n_1)
                                        // attached to each third order random effect
    array[n,T] int<lower=0> y; // Response variable 
    matrix[n*T, K] X; // Predictors include domain-indexed 12 month lagged values,
                    // NO Intercept
                    // Last month's by-record value
                    // Last month's by-domain value
                    // variables used to determine weight stratification, including size class
                    // assume stack 1:n units, each of t = 1,..,T months
    // Input Z in 3 vector CSR (sparse) format.
    // construct 3 vector sparse representation using rstan::extract_sparse_parts(Z[[t]]) in R (not Stan)
    // (w,v,u): w equal to non-zero values in matrix; v = column index of each non-zero value; u = 
    // index of where each row starts in w.
    array[T] int<lower=1> num_wv_vals;// contains the length of  vector of values !ne 0, w[[t]], under CSR 
                                  // (compressed storage representation.)  
    vector[sum(num_wv_vals)] w_csr; // stacked vector of values in each month, t \in 1,..,T sparse form
    array[sum(num_wv_vals)] int<lower=1> v_csr; // column location value for each w value in each month 
                                        // same length as w
    array[(n+1)*T] int<lower=1> u_csr; // row locations of w values with padding on either side for each month
    // matrix[n*T, K_r] Z; // Binary matrix of main and interaction random effect links of units
    //                     // Membership or link to random effect levels may change in each month
    //                     // e.g., main effects: detailed industry, MSA, State, Stratum
    //                     // e.g., interaction effects: industry-MSA, industry-State, Industry-Stratum
    //                     // stack n x K_r matrix T times, month-by-month, rbind
    array[n*T] int<lower=1,upper=T> month; // month assignment for each i*t case.
    array[n*T] int<lower=0,upper=n> unit; // unit assignment for each i*t case
    array[n] int<lower=1> weight_index; // post-stratum index for each respondent - takes n_w unique values
                                  // labeling 1,...,n_w must correspond to index of w
    //real<lower = 0> eta_nw; /* hyperparameter for inducing Sigma_nw ~ LKJ(eta_nw) */
    int<lower=1> n_df;
} /* end data block */


transformed data {
    array[n*T] int y_vec; // vectorized y, including missing values
    array[n_obs] int yobs_vec; // vectorized observed values of y
    vector[0] dat; /* empty dummy offset data vector for L_functor */
    //vector[n*T] lambdayctr_0 = rep_vector(1, n*T); // initial guess for lambda_center_y
    vector[n*T] lambdayctr_0;
    array[2*n_obs] int<lower=0> dat_int;
    
    // vectorize y for modeling
    for(t in 1:T)
    {
      y_vec[(t-1)*n+1:t*n] = y[1:n,t];// column-major order; stack the T (n x 1) columns
    } /* end loop t over T months */
    // extract vector of observed y's for likelihood 
    yobs_vec  = y_vec[ind_obs];
    
    //stacked y_obs and ind_obs for L_functor()
    //to allow definining likelihood only over observed y
   dat_int[1:n_obs]                       = yobs_vec;
   dat_int[(n_obs + 1):(2*n_obs)]         = ind_obs;
   
   lambdayctr_0                           = rep_vector(1,n*T);
   lambdayctr_0[ind_obs]                  = log(to_vector(yobs_vec) + 1);

    
}

// The parameters accepted by the model. 
parameters {
  // vector[n*T] lambda_y; /* latent log-mean of poisson likelihood on y */
  vector<lower=0>[K] sigma_betax; /* standard deviations of K x K covariance matrix n_w, K x 1 rows of beta_x */
  vector<lower=0>[T] sigma_betaz; /* global shrinkage scale of random effects coefficients, K_r x 1, beta_z */
  vector<lower=0>[n_1] lambda_1; /* local scale parameters for first order interaction random effects */ 
  // vector<lower=0>[2] delta; /* scale adjustment for second and third order random effects interaction */
  //                           /* scale. e.g., lambda_2[k] = delta[2] * lambda_1[second_order[k,1],
  //                                                                    lambda_1[second_order[k,2]]; */
  matrix[K,T] betaraw_x; /* coefficients indexed by stratum (as well as by predictor) */
  matrix[K_r,T] betaraw_z; /* coefficients indexed by stratum (as well as by predictor) */
  //cholesky_factor_corr[K] L_beta; /* cholesky of correlation matrix for K x K Sigma_betay */
  //vector<lower=0>[n_w] sigma_nw; // /* standard deviation in marginal distribution for y */
  vector<lower=0>[n_w] sigma2_lambda; // standard deviation of log-mean of y 
} /* end parameters block */

transformed parameters {
  vector[n*T] muvec_z;
  vector[n*T] muvec_x;
  vector[n*T] muvec_y;
  matrix[n,T] mumat_y;
  
  matrix[K,T] beta_x;
  matrix[K_r,T] beta_z;
  // vector<lower=0>[n_2] lambda_2; /* local scale parameters for second order interaction random effects */
  // vector<lower=0>[n_3] lambda_3; /* local scale parameters for third order interaction random effects */
  // vector<lower=0>[K_r] lambda; // lambda = rbind( lambda_1, lambda_2, lambda_3 )
  vector<lower=0>[n_w] sigma_lambda = sqrt(sigma2_lambda); // standard deviation of log-mean of y
  
  /* Allow dependence over n_w stratum dimension AND K predictors in fixed effects coefficients for y */
  // place an RW(1) prior on beta_y, beta_w1, beta_w2
  // e.g., beta_y[t] ~ beta_y[t-1] + N_{n_w x K}(Sigma_nw, Sigma_beta)
  for( k in 1:K )
  {
    beta_x[k,] = cumulative_sum(betaraw_x[k,]);
  } /* end loop k over K fixed effects predictors */
  
  for( k in 1:K_r )
  {
    beta_z[k,] = cumulative_sum(betaraw_z[k,]);
  } /* end loop k over K fixed effects predictors */
  
  // for scale parameters for interaction effects from those for main effects to which they link
  // lambda_2     = delta[2-1] * (lambda_1[second_order[,1]] .* lambda_1[second_order[,2]]);
  // lambda_3     = delta[3-1] * (lambda_1[third_order[,1]]  .* lambda_1[third_order[,2]]
  //                                                         .* lambda_1[third_order[,3]]);
  // 
  // lambda = append_row(append_row(lambda_1,lambda_2),lambda_3); /* K_r x T */
  
  { /* local block for inducing prior on beta_y */
    for( t in 1:T )
    {
      // beta_y[,t]    =  diag_pre_multiply(sigma_beta,L_beta) *
      //                       beta_y[,t];
      beta_z[,t]   = beta_z[,t] .* (lambda_1 * sigma_betaz[t]);
      beta_x[,t]   = diag_matrix(sigma_betax) * beta_x[,t];
    } /* end loop t over T months to scale rw1 betaint_y */
  } /* end local block */
  
  
  for( i in 1:n*T )
  {
    muvec_x[i]   = dot_product(X[i,],to_vector(beta_x[,month[i]]));
  } /* end loop i over n*T respondents */
  
  { /* local block to compute Z_t * beta_z[,t] */
    int start_t = 1;
    // Fill n*T x 1, muvec_z in n x 1 blocks for each month
    for(t in 1:T)
    {
      // extract sparse representation for n x K_r, Z_t from stacked input vectors
      vector[num_wv_vals[t]] Z_tw               = w_csr[start_t:(start_t+num_wv_vals[t]-1)];
      array[num_wv_vals[t]] int Z_tv            = v_csr[start_t:(start_t+num_wv_vals[t]-1)];
      array[n+1] int Z_tu                       = u_csr[(t-1)*(n+1)+1:t*(n+1)];
      start_t                                   += num_wv_vals[t];
      // sparse matrix multiply n x K_r, Z_t with K_r x 1, beta_z[,t]
      muvec_z[((t-1)*n + 1):t*n]                = csr_matrix_times_vector(n, K_r, Z_tw, Z_tv, Z_tu, beta_z[,t]); // n x 1 
    } /* end loop t over T months */
  } /* end local block */
  
  muvec_y              = muvec_x + muvec_z; // n*T x 1
  mumat_y              = to_matrix(muvec_y,n,T); // column major order, taking n x 1 values
  
} /* end transformed parameters block */
  
 model {
  // L_beta          ~ lkj_corr_cholesky(6);
  sigma_betax            ~ student_t(n_df,0,1);
  sigma_betaz            ~ student_t(1,0,1); /* cauchy */
  lambda_1               ~ std_normal();
  // delta                  ~ std_normal();
  to_vector(betaraw_x)   ~ std_normal();
  to_vector(betaraw_z)   ~ std_normal();
  
  sigma2_lambda          ~ gamma(1,1);
  
  
  /* Model likelihood for lambda_y, y, w */
 // lambda_y                ~ normal(muvec_y,sigma_lambda[weight_index[unit]]);
 // yobs_vec                ~ poisson_log(lambday_obs);
 
 target += laplace_marginal_lpmf(dat_int | L_functor, muvec_y, dat,
                                     lambdayctr_0, K_functor, weight_index, 
                                     unit, sigma2_lambda);

} /* end model block */  

generated quantities {
  
  
  matrix[n,T] y_rep; // replicate data for posterior predictive checks
  matrix[n,T] mu_y; // mean of y without over-dispersion parameter - extract this as smoothed y.
  
  { /* local block for vectorized quantities because Chris demanded it */
    vector[n*T] yv_rep;
    
    vector[n*T] lambda_center_y = laplace_marginal_rng(L_functor, muvec_y, dat, dat_int,
                                               lambdayctr_0, K_functor,
                                               forward_as_tuple(weight_index, unit),
                                               forward_as_tuple(weight_index, unit), sigma2_lambda);
                                               /* centered latent response on log scale */
                                            
     vector[n*T] lambda_y       = muvec_y + lambda_center_y; /* latent response on log scale */
     mu_y                       = to_matrix( exp(lambda_y), n, T ); /* mean on the data scale */
    // for( ell in 1:n*T )
    // {
    //    thetavec_y[ell] = expm1( muvec_y[ell] + 0.5*square(sigma_lambda[weight_index[unit[ell]]]) ) + 1;
    // }
    // mu_y                  = to_matrix(thetavec_y,n,T);
    
    for( i in 1:n*T )
    {
      yv_rep[i]        = poisson_log_rng( lambda_y[i] );
    } /* end loop i over n units */

    y_rep   = to_matrix(yv_rep,n,T,1); /* on data scale */
    
  } /* end local block */
  
 
} /* end generated quantities block */
