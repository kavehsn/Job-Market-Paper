library(MASS)
library(mvtnorm)


POS_Dep <-function(y,x,null=c(0,0),level=0.05,p=0.5,B=10000,...){

	if(is.matrix(x)==TRUE){

		if(nrow(x)!=length(y)){

			"The number of rows of matrix X and the length of vector y must be equal!"

		}else{

			n<-length(y)

			z<-matrix(data=c(y,x),nrow=n,ncol=ncol(x)+1)

			smp_size<-floor(0.1*nrow(z))

			z_Alt<-head(z,smp_size)
			z_Test<-tail(z,n-smp_size)

			y_Alt<-z_Alt[,1]
			x_Alt<-z_Alt[,2:ncol(z_Alt)]


			y_Test<-z_Test[,1]
			x_Test<-z_Test[,2:ncol(z_Test)]


			X_Test<-matrix(data=c(rep(1,times=nrow(x_Test)),x_Test),nrow=nrow(x_Test),ncol=ncol(x)+1)


			sgn_y<-1*((y_Test[1:(length(y_Test)-1)]-(X_Test[2:nrow(X_Test),]%*%null))>=0)


			betahat<-rlm(y_Alt[1:(length(y_Alt)-1)]~x_Alt[2:length(y_Alt),])


			a1<-rep(0,times=length(sgn_y))
			b1<-rep(0,times=length(sgn_y))


	
			UniPhi<-pt(X_Test[2:nrow(X_Test),]%*%(betahat$coefficients-null),df=1)
			a1[length(UniPhi)]<-log(1/((1/UniPhi[length(UniPhi)])-1))
			b1[length(UniPhi)]<-0

			mu<-c(0,0)
			Sigma<-diag(2)


			bvrateU<-matrix(data=c(X_Test[2:(nrow(X_Test)-1),]%*%(null-betahat$coefficients),X_Test[3:nrow(X_Test),]%*%(null-betahat$coefficients)),nrow=nrow(X_Test[2:(nrow(X_Test)-1),]),ncol=2)
			
			BiPhi<-rep(0,times=nrow(bvrateU))
			


			for(i in 1:nrow(bvrateU)){

				BiPhi[i]<-pmvt(lower=-Inf,upper=bvrateU[i,],delta=mu,sigma=Sigma,df=1)	
				
			}



			for(k in 1:(length(a1)-1)){

				a1[k]<-log((1-(BiPhi[k]/(1-UniPhi[k+1])))/(BiPhi[k]/(1-UniPhi[k+1])))
				b1[k]<-log((1-(((1-UniPhi[k])/UniPhi[k+1])-(BiPhi[k]/UniPhi[k+1])))/(((1-UniPhi[k])/UniPhi[k+1])-(BiPhi[k]/UniPhi[k+1])))-log((1-((BiPhi[k])/(1-UniPhi[k+1])))/((BiPhi[k])/(1-UniPhi[k+1])))

			}


			POSStat<-sum(a1*sgn_y)+sum(b1[1:(length(b1)-1)]*sgn_y[1:(length(sgn_y)-1)]*sgn_y[2:(length(sgn_y))])


			print(POSStat)
			print(paste("POS_Dependent test statistic: ",POSStat))
			invisible(POSStat)

			POSStat_Sim <- rep(0, times = B)

			for(i in 1:B){

				sgn_y_sim<-rbinom(length(sgn_y),1,p)
				POSStat_Sim[i]<-sum(a1*sgn_y_sim)+sum(b1[1:(length(b1)-1)]*sgn_y_sim[1:(length(sgn_y_sim)-1)]*sgn_y_sim[2:(length(sgn_y_sim))])

			}


			CritVal<-quantile(POSStat_Sim,1-level)

			print(paste("The critical value at ",level," level is", CritVal))

			if(POSStat>CritVal){

				print(paste("Rejected the null hypotehsis at ",level," level"))

			}else{

				print(paste("Failed to Reject the null hypothesis at ",level," level"))

			}




		}
	}else{


		if(length(x)!=length(y)){

			"The lengths of the vectors x and y must be equal!"

		}else{

			n<-length(y)

			z<-matrix(data=c(y,x),nrow=n,ncol=2)

			smp_size<-floor(0.1*nrow(z))

			z_Alt<-head(z,smp_size)
			z_Test<-tail(z,n-smp_size)

			y_Alt<-z_Alt[,1]
			x_Alt<-z_Alt[,2:ncol(z_Alt)]


			y_Test<-z_Test[,1]
			x_Test<-z_Test[,2:ncol(z_Test)]


			X_Test<-matrix(data=c(rep(1,times=length(x_Test)),x_Test),nrow=length(x_Test),ncol=2)


			sgn_y<-1*((y_Test[1:(length(y_Test)-1)]-(X_Test[2:nrow(X_Test),]%*%null))>=0)


			betahat<-rlm(y_Alt[1:(length(y_Alt)-1)]~x_Alt[2:length(y_Alt),])


			a1<-rep(0,times=length(sgn_y))
			b1<-rep(0,times=length(sgn_y))


	
			UniPhi<-pt(X_Test[2:nrow(X_Test),]%*%(betahat$coefficients-null),df=1)
			a1[length(UniPhi)]<-log(1/((1/UniPhi[length(UniPhi)])-1))
			b1[length(UniPhi)]<-0

			mu<-c(0,0)
			Sigma<-diag(2)


			bvrateU<-matrix(data=c(X_Test[2:(nrow(X_Test)-1),]%*%(null-betahat$coefficients),X_Test[3:nrow(X_Test),]%*%(null-betahat$coefficients)),nrow=nrow(X_Test[2:(nrow(X_Test)-1),]),ncol=2)
			
			BiPhi<-rep(0,times=nrow(bvrateU))
			


			for(i in 1:nrow(bvrateU)){

				BiPhi[i]<-pmvt(lower=-Inf,upper=bvrateU[i,],delta=mu,sigma=Sigma,df=1)	
				
			}



			for(k in 1:(length(a1)-1)){

				a1[k]<-log((1-(BiPhi[k]/(1-UniPhi[k+1])))/(BiPhi[k]/(1-UniPhi[k+1])))
				b1[k]<-log((1-(((1-UniPhi[k])/UniPhi[k+1])-(BiPhi[k]/UniPhi[k+1])))/(((1-UniPhi[k])/UniPhi[k+1])-(BiPhi[k]/UniPhi[k+1])))-log((1-((BiPhi[k])/(1-UniPhi[k+1])))/((BiPhi[k])/(1-UniPhi[k+1])))

			}


			POSStat<-sum(a1*sgn_y)+sum(b1[1:(length(b1)-1)]*sgn_y[1:(length(sgn_y)-1)]*sgn_y[2:(length(sgn_y))])


			print(POSStat)
			print(paste("POS_Dependent test statistic: ",POSStat))
			invisible(POSStat)

			POSStat_Sim <- rep(0, times = B)

			for(i in 1:B){

				sgn_y_sim<-rbinom(length(sgn_y),1,p)
				POSStat_Sim[i]<-sum(a1*sgn_y_sim)+sum(b1[1:(length(b1)-1)]*sgn_y_sim[1:(length(sgn_y_sim)-1)]*sgn_y_sim[2:(length(sgn_y_sim))])

			}


			CritVal<-quantile(POSStat_Sim,1-level)

			print(paste("The critical value at ",level," level is", CritVal))

			if(POSStat>CritVal){

				print(paste("Rejected the null hypotehsis at ",level," level"))

			}else{

				print(paste("Failed to Reject the null hypothesis at ",level," level"))

			}


		}


	}
}
