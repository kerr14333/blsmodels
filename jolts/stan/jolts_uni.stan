/* 
one month at a time
*/
  
  data{
    int<lower=1> N;                     // number of domains
    int<lower=0> y[N];                  // set of N observations (direct sample-based domain estimates)
    row_vector<lower=0>[N] cv2_y;       // observed cv2
    row_vector<lower=0>[N] nResp;       // number of respondents
    int<lower=1> P;                     // number of covariates included in matrix x
    matrix[N,P] x;                      // Design matrix for P predictors 
    row_vector<lower=0>[N] Emp;         // offset (e.g., employment)
    
    int<lower=0> N_obs; 		// Number of non-missing values
    int<lower=0> N_miss; 		// Number of missing values
    int<lower=1> ind_obs[N_obs];     	// Vector of non-missing value indices
    int<lower=1> ind_miss[N_miss];     	// Vector of missing value indices

    int<lower=1> R;                     // number of regions

    matrix[N,R] region;			// region assignment
    int<lower=0> y_r[R];                // set of R observations for regions
    row_vector<lower=0>[R] cv2_y_r;     // observed cv2, regions
    real<lower=0> cv2_y_nat;            // observed cv2, national

  } /* end data block */
  
  
  transformed data{
    vector[P] zros_b;
    matrix[N,P] logx;  
    int<lower=0> y_obs[N_obs];             // non-missing ys
    real<lower=0> nResp_nat;               // number of respondents, national
    int<lower = 0> y_nat;                  // y for national level
    row_vector<lower=0>[R] nResp_r;        // number of respondents, regions
    row_vector<lower=0>[N_obs] cv2_y_obs;  // observed cv2
    row_vector<lower=0>[N_obs] nResp_obs;  // number of respondents, observed

    zros_b       = rep_vector(0,P);
    logx         = log(x);

    y_obs	 = y[ind_obs];
    nResp_obs	 = nResp[ind_obs];
    cv2_y_obs    = cv2_y[ind_obs];

    nResp_r      = nResp * region ;
    y_nat        = sum( y_r );
    nResp_nat    = sum(nResp_r) ;

  } /* end transformed parameters block */
  
  
  parameters{
    real<lower=0> sqrt_shape; 
    real<lower=0> sqrt_shape_r; 
    real<lower=0> sqrt_shape_nat; 

    //real<lower=0> sigma_sphi;
    //real<lower=0> sigma_sphi_r;
    //real<lower=0> sigma_sphi_nat;

    real<lower=0> sigma_lam; 
    row_vector[N] lambda; 
    row_vector[N] log_epsilon_raw;

    row_vector<lower=0>[N] sqrt_phi;
    row_vector<lower=0>[R] sqrt_phi_r;
    real<lower=0> sqrt_phi_nat;
    real<lower=0> phi_beta;
    real<lower=0> phi_beta_r;
    // real<lower=0> phi_beta_nat;


    vector[P] beta;
//    cholesky_factor_corr[P] C_b; /* cholesky of correlation matrix for cols of Beta */
    vector<lower=0>[P] sigma_b; /* vector of sd parameters for P x P, Lambda */
    
    /* region */ 
    row_vector[R] log_epsilonr_raw;

    /* National */
    real log_epsilon_nat_raw;

    
  } /* end parameters block */
  
  
  transformed parameters{
    vector[N] xb;
    row_vector<lower=0>[N] fitted_y;
    row_vector<lower=0>[N] mean_y;
    row_vector<lower=0>[N_obs] mean_y_obs;
    row_vector<lower=0>[N] fitted_cv2 ;
    row_vector[N] log_epsilon; 
    row_vector<lower=0>[N] epsilon; 
    
    
// phies
    row_vector<lower=0>[N] phi;
    row_vector<lower = 0>[R] phi_r;  //phies for region
    real<lower = 0> phi_nat;  //National

    row_vector<lower=0>[N_obs] fitted_cv2_obs;
    
    /* region */

    row_vector<lower=0>[R] fitted_y_r;
    row_vector<lower=0>[R] fitted_cv2_r;
    row_vector[R] log_epsilonr; 
    row_vector<lower=0>[R] epsilonr; 
    row_vector<lower=0>[R] mean_y_r;

    /* National */
    real<lower=0> fitted_y_nat;
    real<lower = 0> fitted_cv2_nat;
    real log_epsilon_nat;
    real<lower=0> epsilon_nat;
    real<lower=0> mean_y_nat;


    real<lower=0> shape; 
    real<lower=0> shape_r; 
    real<lower=0> shape_nat; 

    real vbias;
    vbias = 1;

    shape       = square(sqrt_shape);
    shape_r     = square(sqrt_shape_r);
    shape_nat   = square(sqrt_shape_nat);

    /* states */
    phi         = square((sqrt_phi));
    phi_r       = square( sqrt_phi_r );
    phi_nat     = square( sqrt_phi_nat );

    log_epsilon = log_epsilon_raw .* (sqrt_phi) - 0.5*phi;
    epsilon     = exp(log_epsilon); 
    xb          = logx * beta ;
    fitted_y    = Emp .* exp(lambda);
    mean_y      = Emp .* exp(lambda ) .* epsilon;
    mean_y_obs  = mean_y[ind_obs];
    fitted_cv2  = inv(fitted_y) + (exp(phi) - 1) ;
    fitted_cv2_obs  = fitted_cv2[ind_obs];


  
   // add up by regions


    log_epsilonr  = log_epsilonr_raw .* sqrt_phi_r - 0.5*phi_r;
    epsilonr      = exp(log_epsilonr);
    fitted_y_r    = fitted_y * region ;
    mean_y_r      = fitted_y_r .* epsilonr;
    fitted_cv2_r  = inv(fitted_y_r) + (exp(phi_r) - 1) ;
    
   // National
    fitted_y_nat    = sum( fitted_y_r );
    log_epsilon_nat = log_epsilon_nat_raw * sqrt_phi_nat - 0.5*phi_nat;
    epsilon_nat     = exp(log_epsilon_nat);


    mean_y_nat      = fitted_y_nat * epsilon_nat;
    fitted_cv2_nat  = inv(fitted_y_nat) + ( exp(phi_nat) - 1) ;

  } /* end transformed parameters block */
  
  model{
    { /* local block for parameters */
      sigma_lam       ~ student_t( 3, 0.0, 1.0 ); 
      lambda          ~ normal( xb, sigma_lam);
      log_epsilon_raw ~ std_normal();
    } /* end local block for parameters */
      
      
    { /* local variable block  */
      //C_b           ~ lkj_corr_cholesky(4);
      //beta          ~ multi_normal_cholesky( zros_b, diag_pre_multiply(sigma_b,C_b) );
      sigma_b         ~ student_t( 3, 0.0, 1.0 );
      beta            ~ normal( 0, sigma_b );
    } /* end local variable block */
      


     sqrt_shape        ~ std_normal();
     sqrt_shape_r      ~ std_normal();
     sqrt_shape_nat    ~ std_normal();

     //sqrt_phi        ~ std_normal();            
     //sqrt_phi_r      ~ std_normal();             
     //sqrt_phi_nat    ~ std_normal();

     //phi_beta 	       ~ normal(0,0.1);            
     //phi_beta_r        ~ normal(0,0.1);            
     // phi_beta_nat      ~ normal(0,0.1);            
     phi_beta 	       ~ std_normal();            
     phi_beta_r	       ~ std_normal();            

     //sigma_sphi        ~ student_t( 3, 0.0, 1.0 );
     //sigma_sphi_r      ~ student_t( 3, 0.0, 1.0 );
     //sigma_sphi_nat    ~ student_t( 3, 0.0, 1.0 );

     //sqrt_phi 	       ~ normal(phi_beta*inv(sqrt(nResp)),sigma_sphi);            
     //sqrt_phi_r        ~ normal(phi_beta_r*inv(sqrt(nResp_r)),sigma_sphi_r);            
     //sqrt_phi_nat      ~ normal(phi_beta_nat*inv(sqrt(nResp_nat)),sigma_sphi_nat);            

     sqrt_phi 	       ~ normal(phi_beta*inv(sqrt(nResp)),0.1);            
     sqrt_phi_r        ~ normal(phi_beta_r*inv(sqrt(nResp_r)),0.1);            
     //sqrt_phi_nat      ~ normal(phi_beta_nat*inv(sqrt(nResp_nat)),0.1);            
     sqrt_phi_nat      ~ normal(0,0.1);            
 
     y_obs             ~ poisson(mean_y_obs);
     cv2_y_obs         ~ gamma (0.5*shape*nResp_obs, 0.5*shape*nResp_obs .* inv( fitted_cv2_obs ));
 
     log_epsilonr_raw  ~ std_normal();
     y_r               ~ poisson( mean_y_r);
     cv2_y_r           ~ gamma (0.5*shape_r*nResp_r, 0.5*shape_r*nResp_r .* inv( fitted_cv2_r ));

     log_epsilon_nat_raw ~ std_normal();
     y_nat               ~ poisson( mean_y_nat);
     cv2_y_nat           ~ gamma (0.5*shape_nat*nResp_nat, 0.5*shape_nat*nResp_nat * inv( fitted_cv2_nat ));


  } /* end model{} block */


generated quantities{
    row_vector<lower=0>[N] fitted_vrnc;
    row_vector<lower = 0>[R] fitted_vrnc_r;
    real<lower = 0> fitted_vrnc_nat;

    fitted_vrnc = (fitted_y .* fitted_y) .* fitted_cv2 ;
    fitted_vrnc_r = (fitted_y_r .* fitted_y_r) .* fitted_cv2_r; 
    fitted_vrnc_nat = square(fitted_y_nat) * fitted_cv2_nat; 
}
