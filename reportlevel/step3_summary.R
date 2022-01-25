oldwd <- getwd()
require(MASS)
require(sampling)
library(rstan)

library( dplyr )



set.seed(17)


#read in necessary data
dt   = readRDS( file = "data/dataZ.rds")
pop  = readRDS( file = "data/pop.rds")
smpl = readRDS( file = "data/smpl.rds")





  ## pop totals at domain / stratum level (if known, could be used for poststratification)
  pop_mh=as.data.frame(aggregate(cbind(pop$x,pop$y1),list(pop$r,pop$m,pop$h),sum))
  names(pop_mh)=c("r","m","h","X_mh","Y1_mh")
  ## pop mean of x at domain x stratum, to be used as a filler in unsampled cells in our inputs for a stan model
  pop_mhx=as.data.frame(aggregate(pop$x,list(pop$r,pop$m,pop$h),mean))
  names(pop_mhx)=c("r","m","h","mx_mh")
  pop_mh=merge(pop_mh,pop_mhx,by=c("r","m","h"),sort=FALSE)
  head(pop_mh)

  ## unweighted sample totals at domain / stratum level
  sam_mh=as.data.frame(aggregate(cbind(smpl$y2,smpl$y1,smpl$x),list(smpl$r,smpl$m,smpl$h),sum))
  names(sam_mh)=c("r","m","h","y2_mh","y1_mh","x_mh")
  ## compute sample unweighted ratios of y2/y1 at domain x stratum level
  sam_mh$r_mh=sam_mh$y2_mh/sam_mh$x_mh
  head(sam_mh)
  sam_mh$r2_mh=sam_mh$y2_mh/sam_mh$y1_mh
  sam_mh$r1_mh=sam_mh$y1_mh/sam_mh$x_mh
  sam_mh$r2x_mh=sam_mh$y2_mh/sam_mh$x_mh
	     
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
		      #smpl$logx[is.na(smpl$x)]=log(smpl$mx_mh[is.na(smpl$x)])
		      #  smpl$x[is.na(smpl$x)]=smpl$mx_mh[is.na(smpl$x)]

 ### read fitted values from stan output
 fit=readRDS(file=paste("y_countimp_ispbase_010522.rds",sep=""))

 samps  <- extract(fit, pars = "mu_y", permuted = FALSE)[,1,]

 ## posterior means
 smpl$fitted_y1     <- colMeans(samps)[1:dt$n]
 smpl$fitted_y2     <- colMeans(samps)[(dt$n+1):(2*dt$n)]
 print(smpl[smpl$m==44,])


#  pdf("check2.pdf")
#  plot(smpl$y1,smpl$fitted_y1)
#  abline(0,1)
#  dev.off()


 pop =pop[order(pop$m),]
 smpl=smpl[order(smpl$m),]

 N_m = as.vector(by(pop,pop$m,nrow))
 n_m = as.vector(by(smpl,smpl$m,nrow))
 n_h = as.vector(by(smpl,smpl$h,nrow))

 TrueX_m = as.vector(by(pop$x,pop$m,sum))
 TrueY1_m = as.vector(by(pop$y1,pop$m,sum))
 TrueY2_m = as.vector(by(pop$y2,pop$m,sum))

 wx = smpl$w_h*smpl$x
 wy1 = smpl$w_h*smpl$y1
 wy2 = smpl$w_h*smpl$y2

 n_m = as.vector(by(smpl,smpl$m,nrow))
 n_h = as.vector(by(smpl,smpl$h,nrow))



###############################################################################################
  ## (unweighted) ratio estimator for months 1 and 2
  Unw_Y1_m = TrueX_m*as.vector(by(smpl$y1,smpl$m,sum,na.rm=TRUE))/as.vector(by(smpl$x,smpl$m,sum,na.rm=TRUE))
  Unw_Y2_m = Unw_Y1_m*as.vector(by(smpl$y2,smpl$m,sum,na.rm=TRUE))/as.vector(by(smpl$y1,smpl$m,sum,na.rm=TRUE))

  ## Horvitz-Thompson (expansion) estimator for months 1 and 2
  HT_Y1_m = as.vector(by(wy1,smpl$m,sum,na.rm=TRUE))
  HT_Y2_m = as.vector(by(wy2,smpl$m,sum,na.rm=TRUE))

  ## (weighted) ratio estimator for months 1 and 2;  
  ## month 1 estimate is used for month 2 base level; 
  ## here of course we could equivalently use X*(y2/x) since y1 cancels out in WLRY2=WLRY1*(y2/y1)=(X*(y1/x)*(y2/y1))
  ## but it's not so in general since monthly samples change due to nonresponse

  WLR_Y1_m = TrueX_m*as.vector(by(wy1,smpl$m,sum,na.rm=TRUE))/as.vector(by(wx,smpl$m,sum,na.rm=TRUE))
  WLR_Y2_m = WLR_Y1_m*as.vector(by(wy2,smpl$m,sum,na.rm=TRUE))/as.vector(by(wy1,smpl$m,sum,na.rm=TRUE))


  #### a sample based "prototype" for MRP
  ## form (unweighted) ratios at fine (mxh) levels, times (presumably known) population X_mh totals (everuthing at the mh level)
  sam_mh=merge(sam_mh,pop_mh,by=c("r","m","h"),sort=FALSE)
  sam_mh$wr1_mh=as.numeric(as.vector(sam_mh$X_mh))*as.numeric(as.vector(sam_mh$r1_mh))
  sam_mh$wr2_mh=as.numeric(as.vector(sam_mh$wr1_mh))*as.numeric(as.vector(sam_mh$r2_mh))
  sam_mh$wr2x_mh=as.numeric(as.vector(sam_mh$X_mh))*as.numeric(as.vector(sam_mh$r2x_mh))

  ## estimated domain ratio (as sum_m(X_mh*r_mh)/sum_m(X_mh)  
  ## need it since some mh cells may be missing: 
  ratio1_m=as.vector(by(sam_mh$wr1_mh,sam_mh$m,sum))/as.vector(by(sam_mh$X_mh,sam_mh$m,sum))
  ratio2_m=as.vector(by(sam_mh$wr2_mh,sam_mh$m,sum))/as.vector(by(sam_mh$wr1_mh,sam_mh$m,sum))
  ratio2x_m=as.vector(by(sam_mh$wr2x_mh,sam_mh$m,sum))/as.vector(by(sam_mh$X_mh,sam_mh$m,sum))

  ## estimates for months 1 and 2
  Pst_Y1_m = TrueX_m*ratio1_m
  Pst_Y2_m = Pst_Y1_m*ratio2_m
  Pst_Y2x_m = TrueX_m*ratio2x_m

  ### replace missing x's from pop mean
  smpl$x[is.na(smpl$x)]=smpl$mx_mh[is.na(smpl$x)]

  fitted_mh=as.data.frame(aggregate(cbind(smpl$fitted_y2,smpl$fitted_y1,smpl$x),list(smpl$r,smpl$m,smpl$h),sum))
  names(fitted_mh)=c("r","m","h","y2_mh","y1_mh","x_mh")
  fitted_mh$r1_mh=fitted_mh$y1_mh/fitted_mh$x_mh
  fitted_mh$r2_mh=fitted_mh$y2_mh/fitted_mh$y1_mh
  fitted_mh$r2x_mh=fitted_mh$y2_mh/fitted_mh$x_mh
  fitted_mh$wr1_mh=as.numeric(as.vector(pop_mh$X_mh))*as.numeric(as.vector(fitted_mh$r1_mh))
  fitted_mh$wr2_mh=as.numeric(as.vector(fitted_mh$wr1_mh))*as.numeric(as.vector(fitted_mh$r2_mh))
  fitted_mh$wr2x_mh=as.numeric(as.vector(pop_mh$X_mh))*as.numeric(as.vector(fitted_mh$r2x_mh))
  ## estimates for months 1 and 2
  fitted_Pst_Y1_m=as.vector(by(fitted_mh$wr1_mh,fitted_mh$m,sum))
  fitted_Pst_Y2_m=as.vector(by(fitted_mh$wr2_mh,fitted_mh$m,sum))
  fitted_Pst_Y2x_m=as.vector(by(fitted_mh$wr2x_mh,fitted_mh$m,sum))



  ## plot estimators against true pop Y2 levels
  pdf("plot1.pdf")
  plot(TrueY1_m,HT_Y1_m)
  points(TrueY1_m,Unw_Y1_m,pch=15,col="red")
  points(TrueY1_m,WLR_Y1_m,pch=19)
  points(TrueY1_m,Pst_Y1_m,pch=19,col="green")
  points(TrueY1_m,fitted_Pst_Y1_m,pch=19,col="blue")
  abline(0,1)
  dev.off()

  pdf("plot2.pdf")
  plot(TrueY2_m,HT_Y2_m)
  points(TrueY2_m,Unw_Y2_m,pch=15,col="red")
  points(TrueY2_m,WLR_Y2_m,pch=19)
  points(TrueY2_m,Pst_Y2_m,pch=19,col="green")
  points(TrueY2_m,fitted_Pst_Y2_m,pch=19,col="blue")
  #points(TrueY2_m,fitted_Pst_Y2x_m,pch=19,col="brown")
  abline(0,1)
  dev.off()

  ## some basic stats
  stat<-function(estname,est,truth){
	    return(as.data.frame(cbind(estnm=estname,bias=mean(est-truth),
				                                    mad = mean(abs(est-truth)))))
  }

  stat_HT1=stat("HT1",HT_Y1_m,TrueY1_m)
  stat_WLR1=stat("WLR1",WLR_Y1_m,TrueY1_m)
  stat_Unw1=stat("Unw1",Unw_Y1_m,TrueY1_m)
  stat_Pst1=stat("Pst1",Pst_Y1_m,TrueY1_m)
  stat_fitted1=stat("fitted1",fitted_Pst_Y1_m,TrueY1_m)

  print(rbind(stat_HT1,stat_WLR1,stat_Unw1,stat_Pst1,stat_fitted1))

  stat_HT2=stat("HT2",HT_Y2_m,TrueY2_m)
  stat_WLR2=stat("WLR2",WLR_Y2_m,TrueY2_m)
  stat_Unw2=stat("Unw2",Unw_Y2_m,TrueY2_m)
  stat_Pst2=stat("Pst2",Pst_Y2_m,TrueY2_m)
  #stat_Pst2x=stat("Pst2x",Pst_Y2x_m,TrueY2_m)
  stat_fitted2=stat("fitted2",fitted_Pst_Y2_m,TrueY2_m)
  #stat_fitted2x=stat("fitted2x",fitted_Pst_Y2x_m,TrueY2_m)

  print(rbind(stat_HT2,stat_WLR2,stat_Unw2,stat_Pst2,stat_fitted2))
