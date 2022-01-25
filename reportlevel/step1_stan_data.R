library( stringr )
library( dplyr )
library( rstan ) #I use rstan::extract_spare_parts because we had trouble loading
                 #The Z matrix for random effects 

#####################################################################
### function from Chris to define design matrix Z and related matrices
#####################################################################
# f is an R formula
# df is the associated data frame we want to use

# Terrance used a 'horseshoe' prior which required having an
# indexing of which first order effects Z columns are associated 
# with 2nd and 3rd order random effects. The function below
# returns a list with the Z matrix and the indexing described above


z_helper <- function(f, df){

  #remove intercept if it exists
  f <- update(f, ~ . -1 )
  
  #get the terms form the formula
  my.terms <- terms(f)
  
  #the labels contain the labels for hte first, second, third order terms
  labels <- attr(my.terms, "term.labels")
  
  #get the first order variables from formula
  vars <- labels[ !grepl(pattern=":", x = labels) ]
  
  #####Create model matrix Z
  #model.matrix does an n-1 level formulation due to the intercept
  #we hack this a bit so that it gives us all n levels
  cont <- lapply(vars, function(x) contrasts( df[[x]],contrasts=F))
  names(cont) <- vars
  Z <- model.matrix( f , data = df, contrasts.arg = cont)
  
  #get column names
  Z.names <- colnames(Z)
  
  #what order interaction is each
  num_effects <- str_count(Z.names, pattern=":") + 1
  
  #list of main effects
  main <- Z.names[ which( num_effects == 1) ]
  
  #get matrix of second order interactions
  order <- Z.names[ which( num_effects == 2) ]
  sp <- str_split(string = order, ":")
  second_order <- do.call(rbind, lapply( sp, function(y) { sapply( y, function(x) which( x == main) ) }))
  
  #get matrix of third order interactions
  order <- Z.names[ which( num_effects == 3) ]
  sp <- str_split(string = order, ":")
  third_order <- do.call(rbind, lapply( sp, function(y) { sapply( y, function(x) which( x == main) ) }))
  
  ### get n_1, n_2, etc
  #create list containing counts
  count_eff <- as.list( table(num_effects) )
  
  #rename them as we need them
  names(count_eff) <- paste0("n_",names(count_eff))
  
  #remove attributes
  attr(Z,"assign")    <- NULL
  attr(Z,"contrasts") <- NULL
  
  #create list that I will return to the user
  my.list <- c(
    list( second_order = second_order,
          third_order  = third_order,
          Z = Z),
    count_eff )
  #return the final product to the user
  return( my.list )
}
#####################################################################
#####################################################################

## some basic stats
#####################################################################
stat<-function(estname,est,truth){
  return(data.frame(estnm = estname,
                    bias  = mean(est-truth),    
                    mad = mean(abs(est-truth))))
}


### Read in data files ###
pop  = readRDS(file = "data/pop.rds")
smpl = readRDS(file = "data/smpl.rds")


########################################################################
## Part II. Just some checks
#######################################################

### True domain pop values
############################################################################
TrueX_m  = as.vector( by( pop$x, pop$m, sum))     ## true domain totals at month 0
TrueY1_m = as.vector( by( pop$y1, pop$m, sum))   ## true domain totals at month 1
TrueY2_m = as.vector( by( pop$y2, pop$m, sum))   ## true domain totals at month 2
N_m      = as.vector( by( pop, pop$m, nrow))          ## number of population units per domain

### overall
TrueY1 = sum( pop$y1 )
TrueY2 = sum( pop$y2 )
print( TrueY2 - TrueY1 )
############################################################################


wx  = smpl$w_h*smpl$x
wy1 = smpl$w_h*smpl$y1
wy2 = smpl$w_h*smpl$y2


n_r = as.vector(by(smpl,smpl$r,nrow))        ## number of sample units per region
print(n_r)

n_m = as.vector(by(smpl,smpl$m,nrow))        ## number of sample units per domain
print(n_m)

n_h = as.vector(by(smpl,smpl$h,nrow))        ## number of sample units per stratum
print(n_h)

### See how some estimators for month 2 would look like
############################################################################

## (unweighted) ratio estimator for month 2; domains pop total for month 1 are used as known
Unw_Y2_m = TrueX_m*as.vector( by( smpl$y2, smpl$m, sum))/as.vector( by( smpl$x, smpl$m, sum))

## Horvitz-Thompson (expansion) estimator for month 2
HT_Y2_m = as.vector(by(wy2,smpl$m,sum))

## (weighted) ratio estimator for month 2; domains pop total for month 1 are used as known
WLR_Y2_m = TrueX_m*as.vector( by( wy2, smpl$m, sum))/as.vector( by( wx, smpl$m, sum))


stat_HT  = stat( "HT",  HT_Y2_m,  TrueY2_m)
stat_WLR = stat( "WLR", WLR_Y2_m, TrueY2_m)
stat_Unw = stat( "Unw", Unw_Y2_m, TrueY2_m)

print(rbind(stat_HT,stat_WLR,stat_Unw))


###############################################################################################
### Part III. Let us now create inputs for our current STAN models, along with some other variables
###############################################################################################


## pop totals at domain / stratum level (if known, could be used for poststratification)  
pop_mh = as.data.frame( aggregate( cbind( pop$x, pop$y1), list(pop$r,pop$m,pop$h), sum))
names(pop_mh) = c( "r", "m", "h", "X_mh", "Y1_mh")

## pop mean of x at domain x stratum, to be used as a filler in unsampled cells in our inputs for a stan model   
pop_mhx = as.data.frame( aggregate( pop$x, list(pop$r,pop$m,pop$h), mean) )
names(pop_mhx) = c("r", "m", "h", "mx_mh")
 
pop_mh=merge(pop_mh,pop_mhx,by=c("r","m","h"),sort=FALSE)
  
## unweighted sample totals at domain / stratum level
sam_mh = as.data.frame( aggregate( cbind(smpl$y2, smpl$y1, smpl$x), list(smpl$r,smpl$m,smpl$h), sum))
names( sam_mh ) = c("r","m","h","y2_mh","y1_mh","x_mh")

## compute sample unweighted ratios of y2/y1 at domain x stratum level
sam_mh$r_mh=sam_mh$y2_mh/sam_mh$x_mh
head(sam_mh)

nrow(pop_mh)
nrow(sam_mh)

### add lines for cells that are not in the sample (use them to define our Z matrix etc.)
miss_mh=merge(pop_mh,sam_mh,by=c("r","m","h"), all.x = TRUE,sort=FALSE)
miss_mh=miss_mh[is.na(miss_mh$r_mh),c(1:6)]
nrow(miss_mh)
sam_h=aggregate(smpl$w_h,list(smpl$h),mean)
names(sam_h)=c("h","w_h")
smpl=merge(smpl,miss_mh,by=c("r","m","h"),all=TRUE, sort=FALSE)
smpl=select(smpl,-w_h)
smpl=merge(smpl,sam_h,by="h",sort=FALSE)
nrow(smpl) 

### replace missing x's from pop mean
smpl$logx[is.na(smpl$x)]=log(smpl$mx_mh[is.na(smpl$x)])
smpl$x[is.na(smpl$x)]=smpl$mx_mh[is.na(smpl$x)]

### new: 11/17/2021
pop_h=aggregate(pop$x,list(pop$h),mean)
names(pop_h)=c("h","mX_h")
smpl=merge(smpl,pop_h,by="h",sort=FALSE)
nrow(smpl) 


### use function from Chris to define design matrix Z and related matrices
zdf<- smpl %>%
  as_tibble() %>%
  select(r,m,h) %>%
  mutate_all(as.factor)

zobject <- z_helper(~r*h*m, zdf)


#### let us see how one additional estimator would look like (a "prototype" for MRP; postbenchmarked to X_mh )

## form (unweighted) ratios at fine (mxh) levels, times (presumably known) population X_mh totals (everything at the mh level)   
sam_mh=merge(sam_mh,pop_mh,by=c("r","m","h"),sort=FALSE)
sam_mh$wr_mh=as.numeric(as.vector(sam_mh$X_mh))*as.numeric(as.vector(sam_mh$r_mh))
## estimated domain ratio (as sum_m(X_mh*r_mh)/sum_m(X_mh)  )
ratio_m=as.vector(by(sam_mh$wr_mh,sam_mh$m,sum))/as.vector(by(sam_mh$X_mh,sam_mh$m,sum))
## ratio estimator for month 2 (postbenchmarked to Y1_mh); a prototype for MRP
Pst_Y2_m = TrueX_m*ratio_m

## plot estimators against true pop Y2 levels
#plot(TrueY2_m,HT_Y2_m)
#points(TrueY2_m,Unw_Y2_m,pch=15,col="red")
#points(TrueY2_m,WLR_Y2_m,pch=19)
#points(TrueY2_m,Pst_Y2_m,pch=19,col="green")
#abline(0,1)

stat_HT=stat("HT",HT_Y2_m,TrueY2_m)
  stat_WLR=stat("WLR",WLR_Y2_m,TrueY2_m)
  stat_Unw=stat("Unw",Unw_Y2_m,TrueY2_m)
  stat_Pst=stat("Pst",Pst_Y2_m,TrueY2_m)
  
  print(rbind(stat_HT,stat_WLR,stat_Unw,stat_Pst))

  
  ### "linear predictors"
  n=nrow(smpl)
  x_std=(smpl$logx-mean(smpl$logx))/sd(smpl$logx)

  ### matrix X, NO intercept
  X=as.matrix(c(x_std,x_std))
  K=ncol(X)
  
  month = as.vector(c(rep(1,n),rep(2,n)))
  unit = as.vector(c(1:n,1:n))
  y = c(smpl$y1,smpl$y2)
  ind_all=1:length(y)

  # vector of non-missing indices in n*T x 1 vector, y
  ind_obs= as.vector(ind_all[!is.na(y)])
  ind_miss=sort(setdiff(ind_all,ind_obs))
  #print(ind_obs)
  #print(ind_miss)
  
  
  #### some additional variables that are NOT used in the MRP model
  #### indices to use for specific random effects 
  vecw=unique(smpl$w_h)
  n_w=length(vecw)
  weight_index=1:nrow(smpl)
  for (i in 1:n_w)
    weight_index[smpl$w_h==vecw[i]]=i
  
  #### for domain random effects
  vecm = unique(smpl$m)
  M=length(vecm)
  domain=1:nrow(smpl)
  for (i in 1:M)
    domain[ smpl$m == vecm[i] ] = i
  ########################################
  
   
  csr_rep <- rstan::extract_sparse_parts(zobject$Z)
  
  #Get length of w for each month
  
  
  # int<lower=1> num_wv_vals[T];// contains the length of  vector of values !ne 0, w[[t]], under CSR 
  # // (compressed storage representation.)  
  # vector[sum(num_wv_vals)] w; // stacked vector of values in each month, t \in 1,..,T sparse form
  # int<lower=1> v[sum(num_wv_vals)]; // column location value for each w value in each month 
  # // same length as w
  # int<lower=1> u[(n+1)*T]; // row locations of w values with padding on either side for each month

  MSA_Stratum =   as.numeric(interaction(as.factor(domain),as.factor(weight_index)))
  #print(colnames(zobject$Z))
  
  stan_data <- list(
    T=2,
    n   = nrow( smpl ),                  # Number of respondents
    M   = length( unique( smpl$m ) ),    # Number of domains
    K = ncol(X),                         # Number of columns in X
    X = X,                               # Design matrix of covariates (no intercept)
    K_r = ncol(zobject$Z),                  # Number of columns in Z 
    #Z = as.matrix(rbind(zobject$Z,zobject$Z)), # Design 0/1 matrix for random effects, including interactions
    num_wv_vals = c(length( csr_rep$w),length( csr_rep$w)),
    w_csr = c( csr_rep$w, csr_rep$w),
    v_csr = c( csr_rep$v, csr_rep$v),
    u_csr = c( csr_rep$u, csr_rep$u),
    n_1 = zobject$n_1,                      # number of main effects in Z matrix
    n_2 = zobject$n_2,                      # number of 2-way interactions
    n_3 = zobject$n_3,                      # number of 3-way interactions
    second_order = zobject$second_order,    # utility n_2x2 matrix for a specific stan script, used to design a prior for variances 
    third_order = zobject$third_order,      # utility n_3x3 matrix for a specific stan script, used to design a prior for variances
    y = as.matrix(cbind(smpl$y1,smpl$y2)),    # nxT matrix of observed y's 

    n_obs=length(ind_obs),                    # numver of non-missing y's
    ind_obs=ind_obs,                          # vector of non-missing indices in n*T x 1 vector, y    
    
    ### other variables
    
    w = as.matrix(cbind(smpl$w_h,smpl$w_h)),  # sampling weights
    #M = length( unique( domain ) ),  	      # number of unique domains
    n_w = length( unique( smpl$w_h ) ),       # number of unique weights (actually, it's the number of strata)
    MSA = domain,                             # domains indices
    MSA_Stratum = MSA_Stratum,		      # domain x stratum interaction
    certainty_index = as.numeric(round(smpl$w_h)==1), # indices for certainty strata
    weight_index = weight_index,                      # unique weights indices
    month = month,                            # month indices
    unit  = unit ,                            # units indices
    eta_nw = 6,                               # hyperparameter for inducing Sigma_nw ~ LKJ(eta_nw)
    n_df=3,                                   # degrees of freedom for inducing sigmas
    mX_h = as.vector(smpl$mX_h)
  )
  
  
  saveRDS(stan_data, "data/stan_data.rds")
