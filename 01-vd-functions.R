#######################################################
##Load required libraries
#install.packages(c("reshape2", "plyr", "dplyr", "caTools", "refund", "mgcv", "patchwork", "gifski"))
library(reshape2)
library(plyr)
library(dplyr)
library(caTools)
library(refund)
library(mgcv)
library(patchwork)
library(gifski)

#load data from Refund
data(sofa)
#######################################################

#######################################################
###set weights (depending on domain length)
quadWeights<- function(argvals, method = "trapezoidal"){
	ret <- switch(method,
			trapezoidal = {
				D <- length(argvals)
				if(D==2){1/2*c(argvals[2] - argvals[1], argvals[D] - argvals[D-1])}else{1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])}
				},
			midpoint = c(0,diff(argvals)), 
	stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
 return(ret)  
}



#######################################################
#Creates a new data frame to use as prediction data set for covariance slice prediction	
create_new_data<-function(maxT,upperF=F,Hz=NULL,argvals=NULL,includezero=T){
	stopifnot(!is.null(Hz)||!is.null(argvals))
	if(is.null(argvals)){
		if(includezero){
			time_vals<-seq(0,maxT, by=1/Hz)			
		}else{
			time_vals<-seq(0,maxT, by=1/Hz)[-1]
		}
	}else{
		time_vals=argvals
	}
	npts<-length(time_vals)
	stime_new<-rep(time_vals,npts)
	ttime_new<-rep(time_vals,each=npts)
	maxT_new<-rep(maxT,length(stime_new))
	newd<-data.frame(stime=stime_new,ttime=ttime_new,maxT= maxT_new)
	if(upperF){
		newd<-newd[-which(newd$stime > newd$ttime),]
	}
	return(newd)
}

#######################################################
####predict covar mtx for given times M (variable times)
pcov_m<-function(covx, return_nonsym=FALSE, times=NULL, argvals=NULL, Hz=NULL, fit_method="REML", upper=FALSE, includezero=T){
	stopifnot(!is.null(Hz)||!is.null(argvals))
	newcov_list<-list()
	nonsym_list<-list()
	counter=1
	#i=times[1]
	for(i in times){
		temp<-create_new_data(i, upperF=upper, Hz=Hz, argvals=argvals, includezero=includezero)
		temp_pred<-predict(covx,temp)
		if(is.null(argvals)){
			if(includezero){
				pred_mat<-matrix(temp_pred,ncol=length(seq(0,i, by=1/Hz)))
			}else{
				pred_mat<-matrix(temp_pred,ncol=length(seq(0,i, by=1/Hz)[-1]))	
			}
		}else{
			pred_mat<-matrix(temp_pred,ncol=length(argvals))
		}
		
		if(upper){
			ndim<-dim(pred_mat)[1]
			new_mat<-matrix(0,nrow=ndim,ncol=ndim)
			new_mat[upper.tri(new_mat,diag=F)]<-pred_mat[upper.tri(new_mat,diag=F)]
			new_mat<-new_mat+t(new_mat)
			diag(new_mat)<-diag(pred_mat)
			newcov_list[[counter]]<-new_mat
		}else{
			newcov_list[[counter]]<-(pred_mat+t(pred_mat))/2
		}	
		nonsym_list[[counter]]<-pred_mat
		counter=counter+1
	}
	if(return_nonsym){
		return(list(sym_cov=newcov_list, nonsym_cov=nonsym_list))
	}else{
		return(newcov_list)
	}
}


#######################################################
###Calculate PC's for given domain lengths M (xtimes)
get_pcs_M <- function(xtimes,covfx,upperF=F,npcs=NULL, Hz = NULL,pve=0.99,includezero=T, pos_evals=T){

	nsub<-length(xtimes)

	#initialize return vectors/matrices
	efx_list<-list()
	npc_v<-integer(nsub)
	evals_list<-list()
	cumvar_list<-list()

	for(i in 1:nsub){
		maxt=xtimes[i]	
		argvals=if(includezero){seq(0,maxt,by=1/Hz)}else{seq(0,maxt,by=1/Hz)[-1]}

		#predict covar mtx for subject i (based on maxT = M)
		npc.0<-pcov_m(covfx, times=maxt,argvals=argvals,  upper=F, includezero=includezero)[[1]]

		#transform predicted mtx ()
		qw <- quadWeights(argvals)
		Wsqrt <- diag(sqrt(qw))
   		Winvsqrt <- diag(1/(sqrt(qw)))
   		V <- Wsqrt %*% npc.0 %*% Wsqrt

    	#get eigenvals/vectors from estimated cov mtx
    	evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
   		evalues2 = replace(evalues, which(evalues <= 0), 0)
   		var_pct<-evalues2/sum(evalues2)
   		cumvar<-numeric(length(var_pct))
		for(k in 1:(length(var_pct))){
			if(k==1){
				cumvar[k]=var_pct[k]	
			}else{
				cumvar[k]=cumvar[k-1]+var_pct[k]	
			}
		}
    	sub_npcs = ifelse(is.null(npcs), min(which(cumsum(evalues2)/sum(evalues2) > pve)), npcs)
    	sub_npcs = min(sub_npcs, length(argvals))
    	efunctions = matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[,seq(len = sub_npcs)], nrow = length(argvals), ncol = sub_npcs)

  
	npc_v[i]<-sub_npcs
	efx_list[[i]]<-efunctions
	evals_list[[i]]<-eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:sub_npcs]
	cumvar_list[[i]]<-cumvar	
}
	pc_list<-list(evalues=evals_list, efunctions=efx_list, npc=npc_v, cvar=cumvar_list)
	return(pc_list)
}

#######################################################
###Setup data to run trivariate smoother and estimate covariance fx
datsetup_cov<-function(xi_ord){
	nobs_x <- ddply(xi_ord,~id,summarize, nobs_x=length(y))$nobs_x
	nsub<-length(unique(xi_ord$id))
	subj_vec<-rep(1:nsub, nobs_x)
	xi_fit1<-select(xi_ord,-y,-mean)
	xdat2<-ddply(xi_fit1,~id,transform, kprod=kronecker(ydiff,ydiff))
	colnames(xdat2)[3]<-"stime"
	temp<-rep(nobs_x, nobs_x)
	ttime_x<-rep(xi_fit1$time,temp)
	xdat2$ttime<-ttime_x
	return(xdat2)
}

#######################################################
###Translate domain to [0,1] scale
translate<-function(vec){
  range=max(vec)-min(vec)
  t.vec=(vec-min(vec))/range
  return(t.vec)
}

#######################################################
###Normalize vector by integral of data
qfunc<-function(xvec,argvals){
		xvec/sqrt(trapz(argvals,xvec^2))
}
