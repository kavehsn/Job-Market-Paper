## If you encounter problems installing the devtoold package, follow the relevant dependencies using 
## sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

install.packages("devtools")


library(devtools)

install_github("kavehsn/PredictiveDGP")


library(SignTestsDufour)
library(PredictiveDGP)

POS_Comparison<-function(T,coef_min=0,coef_max,incr,iter=5000,...){

	matSize <-length(seq(coef_min,coef_max,by=incr))
	simResult <- matrix(data=0, ,nrow=matSize, ncol=2)

	coef_min <- 0
	matCounter <- 1

	while(coef_min<=coef_max){

		Rej_POS_Dep<-0

		for(i in 1:iter){

			Z <- NormalDGP(n=T,beta=coef_min,...)

			y <- Z[,1]
			x <- Z[,2]

			Dec_POS_Dep <- POS_Dep(y,x,simul=TRUE,...)

			Rej_POS_Dep <- Rej_POS_Dep + Dec_POS_Dep  			
		}

		rejRate <- (Rej_POS_Dep/iter)*100

		coefRej <- cbind(coef_min ,rejRate)

		simResult[matCounter,] <- coefRej

		matCounter <- matCounter + 1

		perfCounter<-(coef_min/coef_max)*100

		coef_min <- coef_min + incr

		if(perfCounter%%10==0){

			cat(perfCounter,"% completed!\n")

		}

	}

	Results <<- simResult

}