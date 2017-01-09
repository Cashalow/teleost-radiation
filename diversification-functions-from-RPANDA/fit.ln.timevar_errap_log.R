fit.ln.timevar_errap_log<-function (phylos, tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1, meth = "Nelder-Mead", cst.lamb=F, cst.mu=F, expo.lamb=F, expo.mu=F, fix.f=T, fix.mu=F)

### f.lamb and f.mu are the specifications of the functional form of the variation of lambda and mu with time. The first argument is always time. The other arguments are the parameters of the variation, given as a single vector. lamb_par and mu_par give the initial parameter values used for the optimization. The characterization and parameters are given from the present to  the past #### when fitting simulated trees, remove the NULL trees first using non_null_trees! ###########

###the number of parameters entered in lamb_par and mu_par are the number of parameters used for the aicc calculation, so careful here about what we enter as parameter values

{
	### calculate the number of observations ###
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
	
	#print(nobs)
	
	if ((fix.f==T) & (fix.mu==F))
	{
		init<-c(lamb_par,mu_par)
		p<-length(init)  
   		optimLH <- function(init) {
        lamb_par <- init[1:length(lamb_par)]
        mu_par <- init[(1+length(lamb_par)):length(init)]
        f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
        f.mu.par<-function(x){abs(f.mu(x,mu_par))}
        LH <- getLikelihood.ln.timevar_errap_log(phylos,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
        return(-LH)
    }
    temp<-optim(init, optimLH, method = meth)
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1) , lamb_par=temp$par[1:length(lamb_par)], mu_par=temp$par[(1+length(lamb_par)):length(init)])
		}
		
	else if ((fix.f==T) & (fix.mu==T))
	{
		init<-c(lamb_par)   
    	p<-length(init)
   		optimLH <- function(init) {
        lamb_par <- init[1:length(lamb_par)]
        f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
        f.mu.par<-function(x){abs(f.mu(x,mu_par))}
        LH <- getLikelihood.ln.timevar_errap_log(phylos,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
        return(-LH)
    }
    temp<-optim(init, optimLH, method = meth)
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$par[1:length(lamb_par)])
		}
		
	else{
	init<-c(lamb_par,mu_par,f)
	p<-length(init)   
    
    optimLH <- function(init) {
        lamb_par <- init[1:length(lamb_par)]
        mu_par <- init[(1+length(lamb_par)):(length(init)-1)]
        f<-init[length(init)]
        f.lamb.par<-function(x){abs(f.lamb(x,lamb_par))}
        f.mu.par<-function(x){abs(f.mu(x,mu_par))}
        LH <- getLikelihood.ln.timevar_errap_log(phylos,tot_time,f.lamb.par,f.mu.par,min(abs(f),1),cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
        return(-LH)
    }
    temp<-optim(init, optimLH, method = meth)
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1), lamb_par=temp$par[1:length(lamb_par)], mu_par=temp$par[(1+length(lamb_par)):(length(init)-1)],f=min(abs(temp$par[length(init)]),1))}
    
    return(res)
}