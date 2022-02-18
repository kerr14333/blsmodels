# libraries
library(cmdstanr)
library(jsonlite)

#read in data
stan_data <- read_json( "data/eight_schools.json", simplifyVector=T)

set_cmdstan_path("/work/pathfinder_testing/cmdstan")

#### Run eight schools non-centered version
mod1     <- cmdstanr::cmdstan_model( "stan/eight_schools_noncentered.stan" )

fit1 = mod1$pathfinder(algorithm = "multi", data = stan_data,
                       refresh = 1, num_threads = 20, num_paths = 20, 
                       psis_draws = 2000, init=2)

fit1_samp = mod1$sample( data = stan_data,chains=3, adapt_delta = .99 )

### Run eight schools centered version, the one that goes womp womp womp
mod2     <- cmdstanr::cmdstan_model( "stan/eight_schools_centered.stan" )

fit2 = mod2$pathfinder(algorithm = "multi", data = stan_data,
                       refresh = 1, num_threads = 20, num_paths = 20, 
                       psis_draws = 2000, num_draws=8000, init=2)

fit2_samp = mod2$sample( data = stan_data,chains=3, adapt_delta = .99)
cat("\n\n Non-Centered Pathfinder Summary \n\n" )

cat("\n\n Pathfinder Summary \n\n" )
fit1$summary()

cat("\n\n HMC Summary \n\n" )
fit1_samp$summary()

cat("\n\n Centered  \n\n" )
cat("\n\n Pathfinder Summary \n\n" )
fit2$summary()

cat("\n\n HMC Summary \n\n" )
fit2_samp$summary()
