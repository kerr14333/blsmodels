###############################################################################
#
# For cross-sectional (univariate) version
#
###############################################################################

#libraries
library(cmdstanr)
library(jsonlite)

### Compile model 
mod         <- cmdstanr::cmdstan_model("stan/jolts_uni.stan")

#Read data
in_dt <- read_json("data/jolts_uni.json", simplifyVector = T)

## do some data preparation
in_dt$y[is.na(in_dt$y)]=1
in_dt$y[in_dt$y==0]=1
in_dt$true_v_y[is.na(in_dt$true_v_y)]=0

fit = mod$pathfinder(algorithm = "multi", data = in_dt,
    refresh = 1, num_threads = 12, num_paths = 12, psis_draws = 2000,init=1)

