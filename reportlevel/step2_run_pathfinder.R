

library(cmdstanr)
set_cmdstan_path("/work/pathfinder/cmdstan")

#### Output file ####
out_file   =   "fit/pathfinder_y_countimp_ispbase_010522.rds"
stan_data = readRDS( file = "data/stan_data.rds")


### Model
mod          <- cmdstanr::cmdstan_model("stan/y_countimp_ispbase_010522.stan")


##### Prepare the inputs ######

stan_data$y[ is.na(stan_data$y) ] = 1  

## CES y values are integer so we force the y's to be so
stan_data$y[ stan_data$y < 1 ] = 1
stan_data$y = round( stan_data$y )

stan_data$w = stan_data$w[,1]

stan_data$num_knots = 10;
y_strat = as.vector( by( stan_data$y[,2], stan_data$weight_index, mean))
stan_data$knots <- unname( quantile( y_strat,probs = seq(from=0, to=1, length.out = stan_data$num_knots),na.rm=TRUE))

stan_data$spline_degree = 3

#### Run Stan ####

fit = mod$pathfinder(algorithm = "multi", data = stan_data,
  		                               refresh = 1, num_threads = 12, num_paths = 12, psis_draws = 2000)

#saveRDS(fit, file=out_file)


