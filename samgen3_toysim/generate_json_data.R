library(jsonlite)

#output path for data
#we write it to a json at the end
outfile <- "samgen3_toysim.json"


set.seed(1793)

#############################################################################
Nst = 50                  ### number of "states"

### assign "fips" to domains (suppose we have 5 groups of states, depending on the number of domains in a state); 
parts = c(25,15,10,5,1)
groups = length(parts)
group = Nst / groups

### number of domains
M = (Nst / groups) * sum(parts)


### standardize the vector of the "number of respondents"
Resp=rep(50,M)
tmp= max(Resp)-min(Resp)
if (tmp==0 ) { nResp = rep.int(1,times=length(Resp)) } else {nResp=(Resp-min(Resp)+1)/tmp}


X = runif(M,-2,2)         ## covariates for point estimates
x_var = exp(rnorm(M,0,1)) ## covariates for variances

true_beta = 1       ### true slopes for point ests regression
true_beta_v = 0.1   ### true slopes for variance model

tau_st = 0.8     ### std of state random effects u_st
tau = 0.6        ### std of domain random effects u_eff

theta = c(0.95,0.05)                ### true cluster probabilities
muk = c(0,2)                        ### true cluster specific intercept

thetaB = c(0.6,0.2,0.2)             ### true cluster probabnilities for bias in variances 
biask = c(2,0.3,1)                  ### multiplicative bias                       


mub = c(0,2)     ### intercept for auxiliary variable B1
bu = c(0,2,4)    ### intercept for auxiliary variable B2

lambda_v = 3    ## to control variability of variance estimates
lambda_g = 4    ## to control variability of true variances around model part


### create vectors of indexes
###################################

fipsindex = NULL  ### for states
for (g in 1:groups){
  tmp = NULL
  for (st in (1+10*(g-1)):(10*g)) 
  {
    tmp = c(tmp,rep(st,parts[g]))
  }
  fipsindex=c(fipsindex,tmp)
}

fipsindex
fipsindex = sample(fipsindex)
length(unique(fipsindex))

#####
clusterindex = NULL  ## for true point ests clusters
for (i in 1:length(theta)) 
  clusterindex = c(clusterindex,rep(i,theta[i]*M))

#####
clusterindexB = NULL  ## for true variance clusters  
for (i in 1:length(thetaB)) 
  clusterindexB = c(clusterindexB,rep(i,thetaB[i]*M))
                                


gen_dt<-function(M=M, Nst=Nst, fipsindex=fipsindex, clusterindex = clusterindex, clusterindexB = clusterindexB,
                     X=X, true_beta = true_beta, 
                     x_var = x_var, true_beta_v = true_beta_v,
                     tau_st=tau_st, tau=tau,
                     theta=theta, thetaB=thetaB,
                     muk=muk,biask=biask,
                     mub=mub,bu=bu,
                     lambda_g=lambda_g,lambda_v=lambda_v){


  fixed=true_beta*X           ## fixed effect part
  u_eff = rnorm(M,0,tau)      ## generate domain random effects
  u_st = rnorm(Nst,0,tau_st)  ## generate state random effects
  
  truth = muk[clusterindex] + fixed + u_st[fipsindex] + u_eff              ### true signal
  
  v_fixed=exp(true_beta_v*log(x_var)-(true_beta_v^2)*0.5)
  print(c("Mean v_fixed", mean(v_fixed)))
  vrnc=1/rgamma(M,(lambda_g+1),lambda_g*v_fixed)  ## true variance of point estimator

  eps=rnorm(M,0,sqrt(vrnc))         ### generate random errors
  eps1=rnorm(M,0,sqrt(vrnc))        ### generate random errors
  
  Y=truth+eps                       ### direct estimate of truth
  Y1=truth+eps1                     ### direct estimate of truth (to be used for checking the coverage)
  
  ## generate "estimates" of variances (depend on number of respondents; outliers also affect the estimates)
  sigma2_y=biask[clusterindexB]*rgamma(M,lambda_v*nResp,lambda_v*nResp/vrnc)
  df=2*lambda_v
  
  print(c("MeanVrnc",median(vrnc),mean(vrnc),median(sigma2_y),mean(sigma2_y)))
  
  
  B1 = rnorm(M,mub[clusterindex],0.5*tau) ## generate auxiliary variable (for point ests)
  B2 = rnorm(M,bu[clusterindexB],0.5*tau) ## generate auxiliary variable (for variance bias)


  dt=as.data.frame(cbind(truth,Y,X,x_var,vrnc,v_fixed,sigma2_y,nResp,df,B1,B2,Y1,fipsindex))
  names(dt)=c("truth","Y","X","x_var","vrnc","v_fixed","sigma2_y","nResp","df","B1","B2","Y1","fipsindex")
  return(dt)
}

#plot(Y,Y1)
#plot(truth,Y1)
#plot(log(vrnc),log(sigma2_y))
#
#plot(log(v_fixed),log(vrnc))

dt = gen_dt(M=M, Nst=Nst, fipsindex=fipsindex, clusterindex = clusterindex, clusterindexB = clusterindexB,
                     X=X, true_beta = true_beta, 
                     x_var = x_var, true_beta_v = true_beta_v,
                     tau_st=tau_st, tau=tau,
                     theta=theta, thetaB=thetaB,
                     muk=muk,biask=biask,
                     mub=mub,bu=bu,
                     lambda_g=lambda_g,lambda_v=lambda_v)

write_json( dt, file = outfile )


