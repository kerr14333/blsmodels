
set.seed(13)

## number of domains (of equal sizes)
M = 50
## "regions"; domains are nested in regions
R = 4 

### 20 strata with varying number of pop units (assume, each stratum includes all domains)
N_mh = c(1,3, 5,10,15,20, 30,50,70, 90, 120, 140, 160, 180, 200, 250, 350,450, 550, 1000) # Num of obs per stratum per domain
N_h = M*N_mh    ## Num of obs per stratum
H = length(N_h) ## Num of strata 
N = sum(N_h)    ## Num of pop units

### set of strata specific intercepts to generate x variable (the population variable involved in design; 
## (will be used to define prob's of selection)
a_h=log(c(2000,1000,500,300,200,100,90,80,70,60,50,40,30,20,10,rep(5,5)))

## set true parameters values
sigma0_h = rep(0.4,H)  ## for x random error
tau0 = 0.2             ## for x domain effect 
phi0 = 0.1             ## for x m*h interaction random effect

sigma1_h = rep(0.4,H)  ## for y[t=1] random error
psi1 = 0.3             ## for y[t=1] stratum effect
tau1 = 0.2             ## for y[t=1] domain effect
phi1 = 0.1             ## for y[t=1] m*h interaction random effect

sigma2_h = rep(0.4,H)  ## for y[t=2] random error
psi2 = 0.3             ## for y[t=2] stratum effect
tau2 = 0.2             ## for y[t=2] domain effect
phi2 = 0.1             ## for y[t=2] m*h interaction random effect


## generate strata effects for month 1 (log scale)
b1_h = rnorm(H,0, psi1)

## generate strata effects for month 2 (log scale)
b2_h = rnorm(H,0, psi2)

## generate domain effects
u0 = rnorm(M,0,tau0)
u1 = rnorm(M,0,tau1)
u2 = rnorm(M,0,tau2)


### Generate finite population
pop = NULL
for (m in 1:M){  ## domain loop
  if (m<=20) {r=1} else      ## domains in region 1
  if (m<=30) {r=2} else      ## domains in region 2
  if (m<=40) {r=3} else      ## domains in region 3 
    {r=4}                    ## domains in region 4
  for (h in 1:H){  ## strata inside each domain
    ## generate random errors
    epsilon0 = rnorm(N_mh[h],0,sigma0_h[h])
    epsilon1 = rnorm(N_mh[h],0,sigma1_h[h])
    epsilon2 = rnorm(N_mh[h],0,sigma2_h[h])

    u0_mh = rnorm(1,0,phi0)
    u1_mh = rnorm(1,0,phi1)
    u2_mh = rnorm(1,0,phi2)
    
    ## collect variables into a frame 
    pop=as.data.frame(rbind(pop,
                            cbind(rep(r,N_mh[h]),
                                  rep(m,N_mh[h]),
                                  rep(h,N_mh[h]),
                                  rep(N_h[h],N_mh[h]),
                                  rep(sigma0_h[h],N_mh[h]),
                                  rep(sigma1_h[h],N_mh[h]),
                                  rep(sigma2_h[h],N_mh[h]),
                                  rep(a_h[h],N_mh[h]),
                                  rep(u0[m],N_mh[h]),
                                  rep(u0_mh,N_mh[h]),
                                  epsilon0,
                                  rep(b1_h[h],N_mh[h]),
                                  rep(u1[m],N_mh[h]),
                                  rep(u1_mh,N_mh[h]),
                                  epsilon1,
                                  rep(b2_h[h],N_mh[h]),
                                  rep(u2[m],N_mh[h]),
                                  rep(u2_mh,N_mh[h]),
                                  epsilon2)
    ))
  }
}
print(nrow(pop)) ## N
names(pop)=c("r","m","h","N_h","sigma0_h","sigma1_h","sigma2_h","a_h","u0","u0_mh","epsilon0","b1_h","u1","u1_mh","epsilon1",
             "b2_h","u2","u2_mh","epsilon2")


### generate population variables
pop$logx = pop$a_h + pop$u0 + pop$u0_mh + pop$epsilon0
## y[t=1] depends on x
pop$logy1 = pop$b1_h + pop$u1 + pop$u1_mh + pop$logx + pop$epsilon1
## y[t=2] depends on y[t=1]
pop$logy2 = pop$b2_h + pop$u2 + pop$u2_mh + pop$logy1 + pop$epsilon2

## set of labels
pop$LABEL = 1:N

## levels
pop$x  = exp(pop$logx)
pop$y1 = exp(pop$logy1)
pop$y2 = exp(pop$logy2)

hist(pop$y2-pop$y1)


## selection probabilities
sd_h = as.vector( by( pop$x, pop$h, sd)) ## find strata std based on "design variable" x
C = 1/sd_h[1]  ## we want the first stratum to be sampled with certainty
pi_h = C*sd_h  ## selection probabilities
w_h = 1/pi_h   ## inverse probs
print( round(w_h,4) )
print( round(pi_h,4) )
## expected number of sampled units
n_h = round( pi_h*N_h )
print( n_h )
sum(n_h)

### function: select Stratified SRS WR
select_sampleSTSRSWOR<-function(H=H,n_h=n_h,pop=pop){
  smpl<-NULL
  for (h in 1:H){
    smpl_h<-cbind(rep(h,n_h[h]),rep(w_h[h],n_h[h]),c(1:n_h[h]),sample(pop$LABEL[pop$h==h],size=n_h[h], replace=FALSE))
    smpl<-as.data.frame(rbind(smpl,smpl_h))
  }
  
  names(smpl)<-c("h","w_h","n","LABEL")
  return(smpl)
}

## Select ST SRS WR
smpl<-as.data.frame(select_sampleSTSRSWOR(H=H,n_h=n_h,pop=pop))

print(nrow(smpl))

## get some of the variable values from the finite pop for units in the sample
smpl$r = pop[pop$LABEL[smpl$LABEL],]$r            # region
smpl$m = pop[pop$LABEL[smpl$LABEL],]$m            # domain
smpl$logx  = pop[pop$LABEL[smpl$LABEL],]$logx     # log(x)
smpl$logy1 = pop[pop$LABEL[smpl$LABEL],]$logy1    # log(y1)
smpl$logy2 = pop[pop$LABEL[smpl$LABEL],]$logy2    # log(y2)

smpl$x  = pop[pop$LABEL[smpl$LABEL],]$x           # x
smpl$y1 = pop[pop$LABEL[smpl$LABEL],]$y1          # y1
smpl$y2 = pop[pop$LABEL[smpl$LABEL],]$y2          # y2

#hist(smpl$y2-smpl$y1)



### save the data
saveRDS(pop, file="data/pop.rds")
saveRDS(smpl, file="data/smpl.rds")


#####################################################################
### end of data simulation
#####################################################################
#####################################################################



########################################################################
## Part II. Just some checks
#######################################################

### True domain pop values
############################################################################
TrueX_m = as.vector(by(pop$x,pop$m,sum))     ## true domain totals at month 0
TrueY1_m = as.vector(by(pop$y1,pop$m,sum))   ## true domain totals at month 1
TrueY2_m = as.vector(by(pop$y2,pop$m,sum))   ## true domain totals at month 2
N_m = as.vector(by(pop,pop$m,nrow))          ## number of population units per domain

### overall
TrueY1 = sum(pop$y1)
TrueY2 = sum(pop$y2)
print(TrueY2-TrueY1)
############################################################################


wx = smpl$w_h*smpl$x
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
Unw_Y2_m = TrueX_m*as.vector(by(smpl$y2,smpl$m,sum))/as.vector(by(smpl$x,smpl$m,sum))
## Horvitz-Thompson (expansion) estimator for month 2
HT_Y2_m = as.vector(by(wy2,smpl$m,sum))
## (weighted) ratio estimator for month 2; domains pop total for month 1 are used as known
WLR_Y2_m = TrueX_m*as.vector(by(wy2,smpl$m,sum))/as.vector(by(wx,smpl$m,sum))

## plot domain estimates against true pop values
#plot(TrueY2_m,HT_Y2_m)
#points(TrueY2_m,Unw_Y2_m,pch=15,col="red")
#points(TrueY2_m,WLR_Y2_m,pch=19)
#abline(0,1)

## some basic stats
stat<-function(estname,est,truth){
  return(as.data.frame(cbind(estnm = estname,
                             bias  = mean(est-truth),    
                             mad   = mean(abs(est-truth))
                             )
                       )
         )
}

stat_HT=stat("HT",HT_Y2_m,TrueY2_m)
stat_WLR=stat("WLR",WLR_Y2_m,TrueY2_m)
stat_Unw=stat("Unw",Unw_Y2_m,TrueY2_m)

print(rbind(stat_HT,stat_WLR,stat_Unw))
