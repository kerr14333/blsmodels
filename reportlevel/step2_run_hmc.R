library(rstan)

#output file for the fit
out_file   =   "fit/hmc_y_countimp_ispbase_010522.rds"

#compile the model
mod          <- stan_model("stan/y_countimp_ispbase_010522.stan", auto_write = TRUE)

## read in the data from the previous step ##
stan_data = readRDS( file = "data/stan_data.rds" )

## 
stan_data$y[ is.na(stan_data$y) ] = 1  
stan_data$y[ stan_data$y < 1] = 1
stan_data$y = round( stan_data$y )

stan_data$w = stan_data$w[,1]

stan_data$num_knots = 10
y_strat = as.vector( by( stan_data$y[,2], stan_data$weight_index, mean))
stan_data$knots <- unname( quantile( y_strat, probs = seq( from=0, to=1, length.out = stan_data$num_knots), na.rm = TRUE))

stan_data$spline_degree = 3


fit <- sampling(object = mod , 
             data = stan_data,  
             pars = c("beta_z", "lambda_y", "mu_y", "lambda_1", "lambda", "sigma_betaz"),
             control=list(max_treedepth=15, adapt_delta = 0.99),
             init_r = 0.1,
             core=4,
             iter = 2000, 
             chains = 3, 
             verbose=TRUE)

saveRDS(fit, file=out_file)
