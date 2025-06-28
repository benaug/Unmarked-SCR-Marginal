sim.USCR.Dcov.multisession <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,xlim=NA,ylim=NA,
           p0=NA,theta.thin=NA,K1D=NA,obstype="poisson"){
    if(length(D.beta0)!=N.session)stop("D.beta0 must be of length N.session")
    if(length(D.beta1)!=N.session)stop("D.beta1 must be of length N.session")
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K)!=N.session)stop("K must be of length N.session")
    if(length(X)!=N.session)stop("X must be of length N.session")
    if(length(buff)!=N.session)stop("buff must be of length N.session")
    if(obstype%in%c("negbin","poisson")){
      if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    }
    if(obstype=="negbin"){
      if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
    }else{
      #make a dummy to pass to sim.RT
      theta.d <- rep(theta.d,N.session)
    }
    J <- rep(NA,N.session)
    for(g in 1:N.session){
      X[[g]] <- as.matrix(X[[g]])
      J[g] <- nrow(X[[g]])
    }
    
    #trap operation
    if(any(is.na(K1D))){
      print("K1D not provided, assuming trap operation is perfect.")
      K1D <- vector("list",N.session)
      for(g in 1:N.session){
        K1D[[g]] <- rep(K[g],J[g])
      }
    }
    
    #simulate sessions one at a time
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.USCR.Dcov(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],InSS=InSS[[g]],res=res[g],
                               lam0=lam0[g],p0=p0[g],sigma=sigma[g],
                               K=K[g],X=X[[g]],obstype=obstype,theta.d=theta.d[g],xlim=xlim[g,],ylim=ylim[g,],
                               K1D=K1D[[g]])
    }
    return(data)
  }