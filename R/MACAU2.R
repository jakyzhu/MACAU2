########################################################################################################################
# Package: MACAU2.0
# Version: 1.0.1
# Date   : 2016-10-15
# Title  : Penalized Quasi-likelihood with Population Sturcture for Counts Data
# Authors: S.Q. Sun, J.Q. Zhu, and X. Zhou
# Contact: shiquans@umich.edu and jiaqiang@umich.edu 
#          University of Michigan, Department of Biostatistics
########################################################################################################################
macau2 <- function( RawCountDataSet, Phenotypes, Covariates=NULL, RelatednessMatrix=NULL, LibSize=NULL, 
                  fit.model="PMM", fit.method = "AI.REML", fit.maxiter=500, fit.tol=1e-5, numCore=1, 
				  filtering=TRUE, verbose=FALSE, ...) {
	if(numCore > 1){
		if(numCore>detectCores()){warning("MACAU2:: the number of cores you're setting is larger than detected cores!");numCore = detectCores()}
	}
	registerDoParallel(cores=numCore)
	
	# filtering counts
	if (filtering & fit.model == "PMM"){
	unfilterIdx <- apply( RawCountDataSet, 1, function(x) length(x[x>5])>=2 )
	CountData   <- RawCountDataSet[unfilterIdx,]
	}else{
		CountData <- RawCountDataSet
	}
	rm(RawCountDataSet)
    
	numVar <- dim(CountData)[1]
	numIDV <- dim(CountData)[2]
	
	# remove the intercept
	if(length(unique(Covariates[,1])) == 1){
		Covariates<-Covariates[,-1]
	}

	if(is.null(Covariates)){
		numCov <- 0
	}else{
		numCov     <- dim(Covariates)[2]
		Covariates <- as.matrix(Covariates)
	}
  
	cat(paste("## number of total individuals: ", numIDV,"\n"))
	cat(paste("## number of total genes/sites: ", numVar,"\n"))
	cat(paste("## number of adjusted covariates: ", numCov,"\n"))
  
  
	CountData  <- as.matrix(CountData)
	Phenotypes <- as.matrix(Phenotypes)
  
  
	if(is.null(RelatednessMatrix)){
		stop("MACAU2::please input relatedness matrix!")
	}else{
		RelatednessMatrix <- as.matrix(RelatednessMatrix)
		scalerM           <- diag(numIDV)-(rep(1,numIDV)%*%t(rep(1,numIDV)))/numIDV
		eig               <- eigen(RelatednessMatrix)
		eigval            <- eig$value
		eigvector         <- eig$vectors
		if(any(eigval<1e-10)){ 
			warning("MACAU2::the relatedness matrix is singular, it has been modified!")
			eigval[eigval<1e-10] <- runif(sum(eigval<1e-10), 0, 0.5)
			RelatednessMatrix    <- eigvector%*%diag(eigval)%*%t(eigvector)
		}
		rm(scalerM)
		rm(eig)
		rm(eigval)
		rm(eigvector)
	}

	RelatednessMatrix <- list( RelatednessMatrix, diag(numIDV))
  
	#***********************************#
	#       Poisson Mixed Model         #
	#***********************************#
	if(fit.model == "PMM"){
		cat("# fitting Poisson mixed model ... \n")
		if(is.null(LibSize)){
			LibSize <- apply(CountData, 2, sum)
			LibSize <- as.matrix(LibSize)
		}else{
			LibSize <- as.matrix(t(LibSize))
		}
	 	 
	 
		# do parallel using foreach function
		resPMM <-foreach(iVar=1:numVar,.combine=rbind)%dopar%{
			numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
			if(numCov==0){
				model0 <- glm(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), family = poisson(link="log"))
				 idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), na.action = na.omit)),
      				      rownames(model.frame(formula = CountData[iVar,]~Phenotypes + offset(log(LibSize)), na.action = na.pass)))
			}else{
				model0 <- glm(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), family = poisson(link="log"))
				 idx   <- match(rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.omit)),
                          rownames(model.frame(formula = CountData[iVar,]~Covariates + Phenotypes + offset(log(LibSize)), na.action = na.pass)))
			}
		
			if(verbose) {cat(paste("NO. Gene = ",iVar,"\n"))}
    
			tmpRelatednessMatrix <- RelatednessMatrix
			if(class(tmpRelatednessMatrix) == "matrix") {
				tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
			}else {
				for(ik in seq_len(length(tmpRelatednessMatrix)) ) {tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]}
			}
      
			names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")

			t1 <- system.time(model1 <- try( MACAU2.fit(model0, tmpRelatednessMatrix) ))
		
			if(class(model1) != "try-error"){
				if(verbose){cat(paste("MACAU2::PMM::tau = ", model1$theta,"\n"))}
				numAnalysis <- length(idx)
				beta        <- model1$coefficients[length(model1$coefficients)]
				se_beta     <- sqrt(diag(model1$cov)[length(model1$coefficients)] )
				pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
				sigma2      <- model1$theta[2]+model1$theta[3]
				h2          <- model1$theta[2]/(sigma2)
				tau1        <- model1$theta[2]
				tau2        <- model1$theta[3]
				converged   <- model1$converged
			}else{converged <- FALSE}
			
			res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta, 
							  pvalue = pvalue, h2 = h2, sigma2 = sigma2, 
							  converged = converged) 
		}# end for iVar, parallel
	
		rownames(resPMM) <- rownames(CountData)
		return(resPMM)
	}# end PMM 
	#***********************************#
	#       Binomial Mixed Model        #
	#***********************************#
	if(fit.model == "BMM"){ 
		cat("# fitting binomial mixed model ... \n")
		if(is.null(LibSize)){
			stop("MACAU2::BMM::ERROR: please input the LibSize file!!")
		}else{
			LibSize <- as.matrix(LibSize)
		}
  
		ratio               <- CountData/LibSize
		ratio[is.na(ratio)] <- 0
		flag                <- ratio>1.0
		sumflag             <- apply(flag,1, sum)
		idx                 <- which(sumflag>0)

		if (length(idx)>0){
			CountData <- CountData[-idx,]
			LibSize   <- LibSize[-idx,]
		}else{
			CountData <- CountData
			LibSize   <- LibSize
		}
		
		numVar <- dim(CountData)[1]
		numIDV <- dim(CountData)[2]	
     
		# do parallel
		resBMM <- foreach(iVar=1:numVar,.combine=rbind)%dopar%{
			numAnalysis <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
			if(verbose){cat(paste("NO. Gene/Site = ",iVar,"\n"))}
			if(sum(dim(LibSize)==dim(CountData)) != 2){
				stop("MACAU2::BMM::ERROR: the dimensions of read counts and total read counts do not match!")
			}
    
			LibSize <- as.matrix(LibSize)
    
			if(numCov == 0){
				model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,])
				idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.omit)),
                          rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Phenotypes, na.action = na.pass)))
			}else{
				model0 <- glm(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, family = binomial(link = "logit"), weights = LibSize[iVar,] )
				idx    <- match(rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.omit)),
						  rownames(model.frame(formula = CountData[iVar,]/LibSize[iVar,]~Covariates + Phenotypes, na.action = na.pass)))
			}
	 	
			model0$numTotal <- LibSize[iVar,idx]
			model0$numSucc  <- CountData[iVar,idx]
     
			redflag <- FALSE
			for( ierr in c(2:dim(model.matrix(model0))[2])){
				if(length(unique(model.matrix(model0)[,ierr])) == 1){
					warning(paste("MACAU2::BMM::the ",ierr-1,"-th column of covariates are the same for gene/site ",rownames(CountData)[iVar],"!",sep = "") )
					redflag <- TRUE
				}
			}
			if(redflag){continue;}
      
			tmpRelatednessMatrix <- RelatednessMatrix
			if(class(tmpRelatednessMatrix) == "matrix") {
				tmpRelatednessMatrix <- tmpRelatednessMatrix[idx, idx]
			}else {
				for(ik in seq_len(length(tmpRelatednessMatrix)) ) {
				tmpRelatednessMatrix[[ik]] <- tmpRelatednessMatrix[[ik]][idx, idx]
				}
			}
			names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")

			t1 <- system.time(model1 <- try( MACAU2.fit(model0, tmpRelatednessMatrix) ))
		
			if(class(model1) != "try-error"){
				if(verbose){cat(paste("MACAU2::BMM::tau = ", model1$theta,"\n"))}
				numAnalysis <- length(idx)
				beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
				se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )
				pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
				sigma2      <- model1$theta[2]+model1$theta[3]
				h2          <- model1$theta[2]/(sigma2)
				tau1        <- model1$theta[2]
				tau2        <- model1$theta[3]
				converged   <- model1$converged
			}else{converged <- FALSE}
		
			res <- data.frame(numIDV = numAnalysis, beta = beta, se_beta = se_beta, 
							pvalue = pvalue, h2 = h2, sigma2 = sigma2, 
							converged = converged)
		}# end for iVar, parallel
	
		rownames(resBMM) <- rownames(CountData)
		return(resBMM)
	}# end BMM
}# end function MACAU2


##########################################################
#           	   MACAU2 FIT FUNCTION					 #
##########################################################
MACAU2.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {
	
		names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
		if( (method.optim == "AI")&(!sum(model0$fitted.values<1e-5))) {
			fixtau.old <- rep(0, length(RelatednessMatrix)+1)
			model1 <- MACAU2.AI(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
			fixtau.new <- 1*(model1$theta < 1.01 * tol)

			while(any(fixtau.new != fixtau.old)) {
				fixtau.old <- fixtau.new
				model1 <- MACAU2.AI(model0, RelatednessMatrix, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
				fixtau.new <- 1*(model1$theta < 1.01 * tol)
			}
			
		if(!model1$converged) {
				print("Average Information REML not converged, using INLA to refit ...\n")
			   if(model0$family$family %in% c("binomial")){
			     y <- model0$numSucc
			   }else{
			     y <- model0$y
			   }
				
				x<- model.matrix(model0)
				numIDV<-length(y)	
				
				g_inla<-1:numIDV
				e_inla<-g_inla+numIDV
			
				eig <- eigen(RelatednessMatrix[[1]])
				eigval <- eig$value
				eigvector <- eig$vectors
				if(any(eigval<1e-10)){ 
				  print("the relatedness matrix is singular, it has been modified!")
				  eigval[eigval<1e-10] <- runif(sum(eigval<1e-10), 0, 0.5)
				  RelatednessMatrix[[1]] <- eigvector%*%diag(eigval)%*%t(eigvector)
				}
			
				Q<-solve(RelatednessMatrix[[1]])
				dataINLA<-data.frame(y=y,x=x)
	
				formula=y~x +f(g_inla, model="generic0", Cmatrix=Q, hyper=list(theta=list(prior="loggamma",param=c(0.5,0.1), initial=10)) )+ f(e_inla, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.1,0.1), initial=10)) )
		
				
				if(model0$family$family %in%  c("binomial")){
					model1=try(MACAU2.INLA(model0=model0,formula=formula,family="binomial",dataset=dataINLA, Q=Q, Ntrials=model0$prior.weights))		
					}
			
				if(model0$family$family %in%  c("poisson")){
					model1=try(MACAU2.INLA(model0=model0,formula=formula,family="poisson",dataset=dataINLA,Q=Q))
					}					
				}
		}else{
		   if(model0$family$family %in% c("binomial")){
		     y <- model0$numSucc
		   }else{
		     y <- model0$y
		   }
			
			x<- model.matrix(model0)
	
			numIDV<-length(y)		
			g_inla<-1:numIDV
			e_inla<-g_inla+numIDV
			
			eig <- eigen(RelatednessMatrix[[1]])
			eigval <- eig$value
			eigvector <- eig$vectors
			if(any(eigval<1e-10)){ 
			  print("the relatedness matrix is singular, it has been modified!")
			  eigval[eigval<1e-10] <- runif(sum(eigval<1e-10), 0, 0.5)
			  RelatednessMatrix[[1]] <- eigvector%*%diag(eigval)%*%t(eigvector)
			}
			Q<-solve(RelatednessMatrix[[1]])
			dataINLA<-data.frame(y=y,x=x)
			
			formula=y~x +f(g_inla, model="generic0", Cmatrix=Q, hyper=list(theta=list(prior="loggamma",param=c(0.5,0.1), initial=10)) )+ f(e_inla, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.1,0.1), initial=10)) )

			if(model0$family$family %in%  c("binomial")){
				model1=try(MACAU2.INLA(model0=model0, formula=formula, family="binomial", dataset=dataINLA, Q=Q, Ntrials=model0$prior.weights))	
				}
		
			if(model0$family$family %in%  c("poisson")){
				model1=try(MACAU2.INLA(model0=model0, formula=formula, family="poisson", dataset=dataINLA, Q=Q))
				}					
		}
	return(model1)
}

##########################################################
#              	MACAU2 INLA FUNCTION					 #
##########################################################
MACAU2.INLA<-function(model0,formula,family,dataset,Q,Ntrials=NULL,...){
	
    if(model0$family$family %in% c("binomial")){
      y <- model0$numSucc
    }else{
      y <- model0$y
    }
	
	numIDV<-length(y)		
	g_inla<-1:numIDV
	e_inla<-g_inla+numIDV
	
	formula=formula
	
	if(family %in%  c("binomial")){	
		res_inla<-try(inla(formula=formula,family="binomial",data=dataset,Ntrials=Ntrials,control.compute=list(dic=T)) )
	}
	
	if (family %in%  c("poisson")){
		res_inla<-try(inla(formula=formula,family="poisson",data=dataset,E=exp(model0$offset),control.compute=list(dic=T)) )
	}
	
	theta<-c()
	theta[1]<-1
	sigma.g_inla= inla.tmarginal(function(x) 1/x,res_inla$marginals.hyperpar$`Precision for g_inla`)
	theta[2]=inla.emarginal(function(x) x, sigma.g_inla)
	
	
	sigma.e_inla= inla.tmarginal(function(x) 1/x,res_inla$marginals.hyperpar$`Precision for e_inla`)
	theta[3]=inla.emarginal(function(x) x, sigma.e_inla)
	
	alpha<-res_inla$summary.fixed[,1]
	eta<-res_inla$summary.linear.predictors
	mu<-res_inla$summary.fitted.values
	cov<-diag((res_inla$summary.fixed[,2])^2)
	converged=(res_inla$mode$mode.status==0)
	
	model1=list(theta = theta, coefficients = alpha, linear.predictors = eta, fitted.values = mu, cov = cov, converged = converged)
	
	return(model1)
}

##########################################################
#       MACAU2 FIT AVERAGE INFORMATION FUNCTION			 #
##########################################################
MACAU2.AI <- function(model0, RelatednessMatrix, tau = rep(0, length(RelatednessMatrix)+1), fixtau = rep(0, length(RelatednessMatrix)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {

  if(model0$family$family %in% c("binomial")){
    y <- model0$numSucc
  }else{
    y <- model0$y
  }
	numIDV <- length(y)
	offset <- model0$offset
	if(is.null(offset)) {offset <- rep(0, numIDV)}
	
	family <- model0$family
	eta <- model0$linear.predictors
	mu <- model0$fitted.values
	mu.eta <- family$mu.eta(eta)
   D <- mu.eta/sqrt(model0$family$variance(mu))
  
	if(family$family %in% c("binomial")){
	  mu.eta <- model0$numTotal*mu.eta
	  D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
	  mu <- model0$numTotal*mu
	}

	Y <- eta - offset + (y - mu)/mu.eta	
	X <- model.matrix(model0)
	alpha <- model0$coef
	
	if(family$family %in% c("poisson", "binomial")) {
		tau[1] <- 1
		fixtau[1] <- 1
	}
	numK <- length(RelatednessMatrix)
	idxtau <- which(fixtau == 0)
	numK2 <- sum(fixtau == 0)
	if(numK2 > 0) {
	  tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)  
	  
		H <- tau[1]*diag(1/D^2)
		for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}
	
		Hinv <- chol2inv(chol(H))
		HinvX <- crossprod(Hinv, X)
		XHinvX <- crossprod(X, HinvX)

		P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
	
		if(class(P) == "try-error"){
			stop("Error in P matrix calculation!")
		}

		PY <- crossprod(P, Y)
		tau0 <- tau
		for(ik in 1:numK2) {
		        if(ik == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
			else {
				PAPY <- crossprod(P, crossprod(RelatednessMatrix[[idxtau[ik]-1]], PY))
				tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(P*RelatednessMatrix[[idxtau[ik]-1]]))/numIDV)
			}
		}
	} 
	
	for (iter in seq_len(maxiter)) {	
		alpha0 <- alpha
		tau0 <- tau
		model1 <- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)
		
		tau <- as.numeric(model1$tau)
		cov <- as.matrix(model1$cov)
		alpha <- as.numeric(model1$alpha)
		eta <- as.numeric(model1$eta) + offset
		
		
		mu <- family$linkinv(eta)
		mu.eta <- family$mu.eta(eta)
		D <- mu.eta/sqrt(family$variance(mu))

		if(family$family %in% c("binomial")){
		  mu.eta <- model0$numTotal*mu.eta
		  D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
		  mu <- model0$numTotal*mu
		}
		
		Y <- eta - offset + (y - mu)/mu.eta

		if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
		if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {
			
			iter <- maxiter
			break
		}
	}

	converged <- ifelse(iter < maxiter, TRUE, FALSE)
	res <- y - mu
	P <- model1$P
	return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, P = P, residuals = res, cov = cov, converged = converged))
}# end function


#########################################
#             CODE END                  #
#########################################
