 
  ## load libraries
  library(cmdstanr) 
  library(jsonlite)
  
  #location/name of stan file 
  stan_file <-"samgen3_brate.stan"
  
  ## Read in simulated data that Julie created
  dt = read_json(path="samgen3_toysim.json",simplifyVector=T)

  ## Set cmdstan path (version with pathfinder)
  #cmdstanr::set_cmdstan_path("/work/pathfinder_testing/cmdstan")
  cmdstanr::set_cmdstan_path("/work/pathfinder_testing/cmdstan")


  #### Start transforming data ####

  M   <- nrow(dt)                     ## Number of estimation cells
  K   <- 6                            ## assumed number of clusters
  Nst <- length(unique(dt$fipsindex)) ## number of distinct geographic area

  # each estimation cell has two quantities an estimate Y and an estimate 
  # of survey variance sigma2_y

  ## transform Y by subtracting the mean and dividing by 
  ## the square root of average surery variances
  mY      <- mean( dt$Y )
  msd     <- sqrt( mean( dt$sigma2_y ) )
  y       <- ( dt$Y - mY )/msd
  
  # we also scale each of the survey variance estimates simliarly
  sigma_y <- sqrt( dt$sigma2_y )/msd

  # Create and scale our X matrix
  mX      <- mean( dt$X )
  x       <- as.matrix( cbind( rep(1,M), (dt$X-mX)/msd) )
  x_var   <- cbind( rep(1,M), (log(dt$x_var)-mean(log(dt$x_var)))/sd(log(dt$x_var)))

    
  oB      <- as.matrix( dt$B1, dt$B2)
  moB     <- colMeans(oB)
  MoB     <- matrix(rep(moB,times <- nrow(oB)),nrow(oB),ncol(oB),byrow<-TRUE)
  soB     <- apply(oB,2,sd)
  SoB     <- matrix(rep(soB,times <- nrow(oB)),nrow(oB),ncol(oB),byrow<-TRUE)
  B       <- (oB - MoB) / SoB  ## N x Bd

  # create data for Stan from the transformed
  # data we made above
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
  
  ## compile model
  mod1     <- cmdstanr::cmdstan_model(stan_file ) 

  #Run pathfinder
  fit1 = mod1$pathfinder(algorithm = "multi", data = stan_data,
			   refresh = 1, num_threads = 12, num_paths = 12, psis_draws = 2000)
 

  fit1$summary()
