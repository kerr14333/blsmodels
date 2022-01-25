# libraries
library(cmdstanr)
library(jsonlite)

#read in data
stan_data <- read_json( "data/eight_schools.json", simplifyVector=T)

#### Run eight schools non-centered version
mod1     <- cmdstanr::cmdstan_model( "stan/eight_schools_noncentered.stan" )

fit1 = mod1$pathfinder(algorithm = "multi", data = stan_data,
                       refresh = 1, num_threads = 12, num_paths = 12, 
                       psis_draws = 2000)

### Run eight schools centered version, the one that goes womp womp womp
mod2     <- cmdstanr::cmdstan_model( "stan/eight_schools_centered.stan" )

fit2 = mod2$pathfinder(algorithm = "multi", data = stan_data,
                       refresh = 1, num_threads = 12, num_paths = 12, 
                       psis_draws = 2000)

