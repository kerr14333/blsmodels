
#libraries
library(cmdstanr)

#set working directory if necessary
#setwd() 

#set path of cmdstan
set_cmdstan_path("/work/laplace_testing/cmdstan")

#get data for model
stan_data <- readRDS("data/data_simpleZ.rds")

#y has some NA's in it so we set them 
#to 1 to run the model
stan_data$y[is.na(stan_data$y)]=1
stan_data$y[stan_data$y<1]=1
stan_data$y=round(stan_data$y)
stan_data$w=stan_data$w[,1]


#get model file
stan_file <- "model/y_mrp_simpleZ_laplace_042722.stan"

#compile model
mod_mrp        <- cmdstan_model(stan_file)

# tuning parameters for MCMC
num_chains <- 4
num_warm <- 500
num_post <- 500

#sample from model
fit <- mod_mrp$sample( data = stan_data,
                       chains = num_chains, 
                       parallel_chains = num_chains,
                       iter_warmup = num_warm, 
                       iter_sampling = num_post,
                       seed = 123, refresh = 0)
