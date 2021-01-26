
#####################################################################################################
# Code used to compare our proposed additive transformation model to the existing method proposed by
# Cheng and Wang (2011) to current status data
#
#Z <-    # Usual design matrix for linear effects (Note first column should be 1s see pg 1740 of the paper)
#W <-    # Continuous covariate for 1st non-linear effect
#U <-    # Continuous covariate for 2nd non-linear effect
#V <-    # Observation time
#
#F(s)=1-exp(-exp(s))     # Taken as this so as to have a partial linear Cox model

#########################################################################
Simulation <- function(N,opt,alpha,v.seq,w.seq,u.seq){
  # N = 200; opt = 1
	
	##################################################
	# ================================================
	# Set time limite
	setTimeLimit(elapsed = 2400, transient = TRUE)
	# ================================================
	##################################################
	
	Nk    <- ceiling(N^(1/3))
	dat   <- data.gen.CS.alpha(N=N,alpha,opt,v.seq,w.seq,u.seq)
	rp    <- 1 - mean(dat$delta)
	d 		<- dat$delta
	V 		<- dat$Ci
	Z 		<- cbind(1,dat$Z)
	W 		<- dat$W1
	U 		<- dat$W2
	adj1  <- dat$adj1
	adj2  <- dat$adj2

	# ===================================================================
	# Our Proposed method
	mm0    <- Sys.time()
	degree <-3
	order  <-3
	fit    <- TS.fit(N,Nk,dat,degree,order,w.seq,u.seq,v.seq)
	mm1    <- Sys.time()


	# ===================================================================
	# Cheng and Wang method
	lower <- Nk - 2
	upper <- Nk + 2
  
	
	tt0   <- Sys.time()
	res 	<- BIC.Cheng.Wang.SAT(lower,upper,alpha,d,Z,V,W,U,v.seq,w.seq,u.seq)
	tt1 	<- Sys.time()
	
	# ===================================================================
	# save
	
	p     <- 2
	p1    <- length(v.seq)
	p2    <- length(w.seq)
	p3    <- length(u.seq)
	lln   <- 2*p + p1 + p2 + p3
	
	if(class(res)!="try-error"){
		t1         <- as.numeric(difftime(tt1,tt0,units='secs'))
		se.save    <- sqrt(diag(solve(res$hess))[2:(p+1)])
		b.save		 <- res$b
		Hv.save		 <- res$Hv+res$b0-adj1-adj2
		Phi.w.save <- res$phi.w + adj1
		Phi.u.save <- res$phi.u + adj2
		save1      <- c(b.save,se.save,Hv.save,Phi.w.save,Phi.u.save)
		save3      <- c(res$grid,res$err,res$b0)
	}else{
		t1         <- 0.0
		save1      <- rep(NA,lln)
		save3      <- rep(NA,3+1+1)
	
	}
	
	if(class(fit)!="try-error"){
		m1     <- as.numeric(difftime(mm1,mm0,units='secs'))
		lam2   <- fit$lam2
		save2  <- c(lam2$par[1:p],lam2$se,fit$f0_lam2,fit$f1_lam2,fit$f2_lam2)
	}else{
		m1     <- 0.0
		save2  <- rep(NA,lln)
	}
	
	save3  <- c(save3,adj1,adj2,t1,m1)
	
	return(list(save1=save1,save2=save2,save3=save3))


}

###################################################################
###################################################################
# Load necessary packages

library(Matrix);library(splines);library(numDeriv);library(MASS);library(Iso)
library(methods);library(mgcv)

cluster <- 1

if(cluster==1){
	source("/data/gpfs/home/yliu23/caseoneadd/PM2_functions.txt")
	source("/data/gpfs/home/yliu23/caseoneadd/Cheng_Yan_Add.txt")
}else{
	setwd("C:/Users/liuyanyxy/OneDrive - University of Nevada, Reno/A_UNR/Research/Survival_Analysis/Case1_add_trans_model/CaseOneAdd/")
	source("Cheng_Yan_Add.txt")
	source("PM2_functions.txt")
}



###################################################################
###################################################################
# Once data is observed 

alpha <- 0.5 # PM MODEL


v.seq <- seq(0,10,.01)
w.seq <- seq(-1,1,.01)
u.seq <- seq(-1,1,.01)
p     <- 2

# data store sequence id
ij    <- 1


jj    <- 1
tot   <- 55
while(jj<=tot){
	opt <- 1
	while(opt <= 3){
	  res_200 <- try(Simulation(N=200,opt,alpha,v.seq,w.seq,u.seq),silent=TRUE)
		if(class(res_200)!="try-error"){
		  res_400 <- try(Simulation(N=400,opt,alpha,v.seq,w.seq,u.seq),silent=TRUE)
		}else{
		  res_400 <- res_200
		}

		# ==================================================================================
		# save
		if(class(res_200)!="try-error"&class(res_400)!="try-error"){
			
			if(cluster==1){
				write(res_200$save1,file=paste0("PM_N200_save1_mod",opt,".txt"),append=TRUE)
				write(res_200$save2,file=paste0("PM_N200_save2_mod",opt,".txt"),append=TRUE)
				write(res_200$save3,file=paste0("PM_N200_save3_mod",opt,".txt"),append=TRUE)
				
				write(res_400$save1,file=paste0("PM_N400_save1_mod",opt,".txt"),append=TRUE)
				write(res_400$save2,file=paste0("PM_N400_save2_mod",opt,".txt"),append=TRUE)
				write(res_400$save3,file=paste0("PM_N400_save3_mod",opt,".txt"),append=TRUE)
			}
			#########################################################
			
			print(paste0("Finish dataset ",jj," mod ",opt))
			opt <- opt + 1
		}
	
	}
	
	jj <- jj + 1
}


###################################################################
if(FALSE){

	fit1 <- Cheng.Wang.SAT(d,V,Z,W,U,degree=3,n.knotsV=3,n.knotsW=3,n.knotsU=3,v.seq,w.seq,u.seq)
	solve(fit1$hess)[1:3]
	
	plot(v.seq,dat$ft)
	lines(v.seq,exp(fit1$Hv+fit1$b0-adj1-adj2),col="red")
	lines(v.seq,exp(fit2$f0_lam2),col="blue")


	plot(w.seq,dat$fw)
	lines(w.seq,fit1$phi.w+adj1,col="red")
	lines(w.seq,fit2$f1_lam2,col="blue")

	plot(u.seq,dat$fu)
	lines(u.seq,fit1$phi.u+adj2,col="red")
	lines(u.seq,fit2$f2_lam2,col="blue")

}
