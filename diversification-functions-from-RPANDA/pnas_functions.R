######################################################################################################################################
######################################################################################################################################
### The function fit.gen.bd fits a birth-death process with any form of time-dependence for lambda (the speciation rate) and mu (the extinction rate) to a list of phylogenies assumed to arise from the same underlying birth-death process;
### phylos is the list of phylogenies; there can be only one phylogeny (phylo) to fit, in this case use phylo<-list(phylo); 
### the possibility of fitting a list of phylogenies (assumed to arise from the same birth-death process) allows relaxing the hypothesis that all lineages trace back to a single common ancestor (see Morlon et al. PNAS 2011)
 
### mu can exceed lamba(negative net diversification rate); the phylogeny can be only partially sampled (i.e. f, the sampling fraction, can be <1)

### tot_time is the time duration of the process (if there is only 1 phylogeny to fit, the root length should be included when this information is available)

### f.lamb and f.mu are the specifications of the functional form of the variation of lambda and mu with time, with time measured from the present to the past
### The first argument is always time. The other arguments are the parameters of the variation, given as a single vector. 
### For example for an exponential time-dependence of the speciation rate, use: f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
### If the functional form is constant or exponential, specify it using e.g. cst.lamb=T, etc... as this will speed up computation

### lamb_par and mu_par give the initial parameter values used for the optimization.
### Note: the number of parameters entered in lamb_par and mu_par are the number of parameters used for the aicc calculation.

### the default option is that the sampling fraction is known (specified by f), and so fix.f=T, but this may be modifed.
### the default option is to estimate the extinction rate (fix.mu=F), but this may be modified (e.g. one can force mu to 0 to fit a pure birth model) 

### this function uses the function getLikelihood.gen.bd, Psi and Phi (codes below, see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

fit.gen.bd<-function(phylos, tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1, meth = "Nelder-Mead", cst.lamb=F, cst.mu=F, expo.lamb=F, expo.mu=F, fix.f=T, fix.mu=F)

{

### calculate the number of observations, used for the computation of aicc values ###
	nobs<-0
	max_age<-0
	for (i in 1:length(phylos))
	{
		phylo<-phylos[[i]]
		if (is.numeric(phylo))
		{nobs<-nobs+1}
		else
		{
			nobs<-nobs+2*(Ntip(phylo)-1)
			max_age<-max(max_age,max(node.age(phylo)$ages))
			}
		}
	if (tot_time>max_age){nobs<-nobs+length(phylos)}
	
### case when f is fixed and mu is estimated 	
	
	if ((fix.f==T) & (fix.mu==F))
	{
		init<-c(lamb_par,mu_par)
		p<-length(init)  
   		optimLH <- function(init) {
        	lamb_par <- init[1:length(lamb_par)]
        	mu_par <- init[(1+length(lamb_par)):length(init)]
        	f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
        	f.mu.par<-function(x){abs(f.mu(x,mu_par))}
        	LH <- getLikelihood.gen.bd(phylos,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
        	return(-LH)
    		}
    	temp<-optim(init, optimLH, method = meth)
    	res <-list(LH=-temp$value,aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$par[1:length(lamb_par)],mu_par=temp$par[(1+length(lamb_par)):length(init)])
		}

### case when f and mu are fixed
		
	else if ((fix.f==T) & (fix.mu==T))
	{
		init<-c(lamb_par)   
    		p<-length(init)
   		optimLH <- function(init) {
        	lamb_par <- init[1:length(lamb_par)]
        	f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
        	f.mu.par<-function(x){abs(f.mu(x,mu_par))}
        	LH <- getLikelihood.gen.bd(phylos,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
        	return(-LH)
    		}
    	temp<-optim(init, optimLH, method = meth)
    	res <-list(LH=-temp$value,aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$par[1:length(lamb_par)])
		}

### other cases		
	
	else{
		init<-c(lamb_par,mu_par,f)
		p<-length(init)    
    		optimLH <- function(init) {
        	lamb_par <- init[1:length(lamb_par)]
        	mu_par <- init[(1+length(lamb_par)):(length(init)-1)]
        	f<-init[length(init)]
        	f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
        	f.mu.par<-function(x){abs(f.mu(x,mu_par))}
        	LH <- getLikelihood.gen.bd(phylos,tot_time,f.lamb.par,f.mu.par,min(abs(f),1),cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
        	return(-LH)
    		}
    	temp<-optim(init, optimLH, method = meth)
    	res <-list(LH = -temp$value,aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$par[1:length(lamb_par)], mu_par=temp$par[(1+length(lamb_par)):(length(init)-1)],f=min(abs(temp$par[length(init)]),1))
		}
    
    return(res)
}

######################################################################################################################################
######################################################################################################################################
### The function getLikelihood.gen.bd computes the log-likelihood of a list of phylogenies assumed to arise from the same underlying birth-death process
### the likelihood expression is that of Morlon et al. PNAS 2011
### this function uses the function Psi and Phi (codes below, see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

getLikelihood.gen.bd<-function(phylos,tot_time,f.lamb,f.mu,f,cst.lamb=F,cst.mu=F,expo.lamb=F,expo.mu=F)

{
	
	r<-function(x){f.lamb(x)-f.mu(x)}
	indLikelihood<-c()
	nbfinal<-length(phylos)
	nbtips_tot<-0
	
	for (i in 1:nbfinal)
	{
		phylo<-phylos[[i]]
		
		if (is.numeric(phylo))
		{
			indLikelihood<-c(indLikelihood,Psi(0,tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu))
			nbtips_tot<-nbtips_tot+1}
		
		else
		{
			from_past<-cbind(phylo$edge,node.age(phylo)$ages)
			nbtips<-Ntip(phylo)
			nbtips_tot<-nbtips_tot+nbtips
			ages<-rbind(from_past[,2:3],c(nbtips+1,0))
			ages<-ages[order(ages[,1]),]
			age<-max(ages[,2])
	
				for (j in 1:(nbtips-1))
					{
						node<-(nbtips+j)
						edges<-phylo$edge[phylo$edge[,1]==node,]
						tj<-age-ages[edges[1,1],2]
						sj1<-age-ages[edges[1,2],2]
						sj2<-age-ages[edges[2,2],2]
						indLikelihood<-c(indLikelihood,f.lamb(tj)*Psi(sj1,tj,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)*Psi(sj2,tj,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu))
						}
						indLikelihood<-c(indLikelihood,Psi(age,tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu))
			}}
		
		data_lik<-prod(indLikelihood)*f^nbtips_tot
		Phi<-Phi(tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
		final_lik<-data_lik/(1-Phi)^nbfinal
		
	return(log(final_lik))
}

######################################################################################################################################
######################################################################################################################################
### This function computes Phi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

Phi<-function(t,f.lamb,f.mu,f,cst.lamb=F,cst.mu=F,expo.lamb=F,expo.mu=F)

{
	
	if ((cst.lamb==T) & (cst.mu==T))
	{
		lamb<-f.lamb(0)
		mu<-f.mu(0)
		r<-lamb-mu	
		res<-1-exp(r*t)/(1/f+lamb/r*(exp(r*t)-1))
		return(res)
		}
		
	if ((cst.lamb==T) & (expo.mu==T))
	{
		lamb0<-f.lamb(0)
		mu0<-f.mu(0)
		beta<-log(f.mu(1)/mu0)
		r.int<-function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
		r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
		r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
		res<-1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
		return(res)
		}
		
	if ((expo.lamb==T) & (cst.mu==T))
	{
		lamb0<-f.lamb(0)
		alpha<-log(f.lamb(1)/lamb0)
		mu0<-f.mu(0)
		r.int<-function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
		r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
		r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
		res<-1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
		return(res)
		}
		
	if ((expo.lamb==T) & (expo.mu==T))
	{
		lamb0<-f.lamb(0)
		alpha<-log(f.lamb(1)/lamb0)
		mu0<-f.mu(0)
		beta<-log(f.mu(1)/mu0)
		r.int<-function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0/beta*(exp(beta*y)-exp(beta*x))}
		r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
		r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
		res<-1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
		return(res)
		}
	
	else
	{
	r<-function(x){f.lamb(x)-f.mu(x)}	
	r.int<-function(x,y){integrate(Vectorize(r),x,y,stop.on.error=FALSE)$value}
	r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
	r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
	res<-1-exp(r.int(0,t))/(1/f+r.int.int(0,t))
	res}
	}
	
Psi<-function(s,t,f.lamb,f.mu,f,cst.lamb=F,cst.mu=F,expo.lamb=F,expo.mu=F)

######################################################################################################################################
######################################################################################################################################
### This function computes Phi (see Morlon et al. PNAS 2011)
######################################################################################################################################
######################################################################################################################################

{
	if ((cst.lamb==T) & (cst.mu==T))
	{
		lamb<-f.lamb(0)
		mu<-f.mu(0)
		r<-lamb-mu	
		res<-exp(r*(t-s))*(abs(1+(lamb*(exp(r*t)-exp(r*s)))/(r/f+lamb*(exp(r*s)-1))))^(-2)
		return(res)
		}
		
	if ((cst.lamb==T) & (expo.mu==T))
	{
		lamb0<-f.lamb(0)
		mu0<-f.mu(0)
		beta<-log(f.mu(1)/mu0)
		r.int<-function(x,y){lamb0*(y-x)-mu0/beta*(exp(beta*y)-exp(beta*x))}
		r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
		r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
		res<-exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
		return(res)
		}
		
	if ((expo.lamb==T) & (cst.mu==T))
	{
		lamb0<-f.lamb(0)
		alpha<-log(f.lamb(1)/lamb0)
		mu0<-f.mu(0)
		r.int<-function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0*(y-x)}
		r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
		r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
		res<-exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
		return(res)
		}
		
	if ((expo.lamb==T) & (expo.mu==T))
	{
		lamb0<-f.lamb(0)
		alpha<-log(f.lamb(1)/lamb0)
		mu0<-f.mu(0)
		beta<-log(f.mu(1)/mu0)
		r.int<-function(x,y){lamb0/alpha*(exp(alpha*y)-exp(alpha*x))-mu0/beta*(exp(beta*y)-exp(beta*x))}
		r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
		r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
		res<-exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
		return(res)
		}


	else
	{
	r<-function(x){f.lamb(x)-f.mu(x)}	
	r.int<-function(x,y){integrate(Vectorize(r),x,y,stop.on.error=FALSE)$value}
	r.int.0<-function(y){exp(r.int(0,y))*f.lamb(y)}
	r.int.int<-function(x,y){integrate(Vectorize(r.int.0),x,y,stop.on.error=FALSE)$value}
	res<-exp(r.int(s,t))*(abs(1+r.int.int(s,t)/(1/f+r.int.int(0,s))))^(-2)
	return(res)}
	}
	
		
