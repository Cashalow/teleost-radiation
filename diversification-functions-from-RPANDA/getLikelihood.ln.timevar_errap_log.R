getLikelihood.ln.timevar_errap_log<-function(phylos,tot_time,f.lamb,f.mu,f,cst.lamb=F,cst.mu=F,expo.lamb=F,expo.mu=F)
#f.lamb and f.mu must be written from the present to the past
{
	
	r<-function(x){f.lamb(x)-f.mu(x)}
	indLikelihood<-c()
	nbfinal<-length(phylos)
	nbtips_tot<-0
	
	for (i in 1:nbfinal)
	{
		#print(i)
		phylo<-phylos[[i]]
		
		if (is.numeric(phylo))
		{
			indLikelihood<-c(indLikelihood,Psi_timevar_errap(0,tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu))
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
						indLikelihood<-c(indLikelihood,f.lamb(tj)*Psi_timevar_errap(sj1,tj,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)*Psi_timevar_errap(sj2,tj,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu))
						}
						indLikelihood<-c(indLikelihood,Psi_timevar_errap(age,tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu))
			}}
		#print(indLikelihood)
		log_data_lik<-sum(log(indLikelihood))+nbtips_tot*log(f)
		#print(nbtips_tot)
		Phi<-Phi_timevar_errap(tot_time,f.lamb,f.mu,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu)
		#final_lik<-data_lik/(1-Phi)^nbfinal*Phi^Nmin/(1-Phi)
		log_final_lik<-log_data_lik-nbfinal*log(1-Phi)
		return(log_final_lik)}
		