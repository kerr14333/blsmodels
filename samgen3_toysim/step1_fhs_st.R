####################
## load libraries
####################
library(cmdstanr) 
library(jsonlite)

#set cmdstan path
set_cmdstan_path("/work/pathfinder_testing/cmdstan")

#read in data
dt = read_json(path="samgen3_toysim.json",simplifyVector=T)

### compile stan script
mod1     <- cmdstanr::cmdstan_model("fhs_st.stan")

### standardize data where necessary

  M = nrow(dt)
  K = 6          ## assumed number of clusters
  Nst = length(unique(dt$fipsindex))

  mY=mean(dt$Y)
  msd=sqrt(mean(dt$sigma2_y))
  y=(dt$Y-mY)/msd
  sigma_y=sqrt(dt$sigma2_y)/msd
  mX=mean(dt$X)
  #x=as.matrix(cbind(rep(1,M),(dt$X-mX)/msd))
  x=as.matrix((dt$X-mX)/msd)
  x_var=cbind(rep(1,M),(log(dt$x_var)-mean(log(dt$x_var)))/sd(log(dt$x_var)))

  oB=as.matrix(dt$B1,dt$B2)
  moB     <- colMeans(oB)
  MoB     <- matrix(rep(moB,times = nrow(oB)),nrow(oB),ncol(oB),byrow=TRUE)
  soB     <- apply(oB,2,sd)
  SoB     <- matrix(rep(soB,times = nrow(oB)),nrow(oB),ncol(oB),byrow=TRUE)
  B   <- (oB - MoB) / SoB  ## N x Bd


  ### prepare the list for the stan call
  stan_data <- list( N = M,
  	                K = K,
		           Nst = Nst,
		           y = y,
		           sigma_y = sigma_y,
		           dX = ncol(x),
		           x = x,
		           P = ncol(x_var),
		           x_var = x_var,
		           nResp = dt$nResp,
		           Bd = ncol(B),
		           B = B,
                     fipsindex = dt$fipsindex
		         )
 

  fit1 = mod1$pathfinder(algorithm = "multi", data = stan_data,
			   refresh = 1, num_threads = 12, num_paths = 12, psis_draws = 2000)

  fit1$summary()
