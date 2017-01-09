Psi_timevar_errap <- function(s,t,f.lamb,f.mu,f,cst.lamb=F,cst.mu=F,expo.lamb=F,expo.mu=F)
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
