e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

getArea <- function(X=X,buff=buff){
  N.session=length(X)
  area=rep(NA,N.session)
  for(a in 1:N.session){
    xlim <- c(min(X[[a]][,1]),max(X[[a]][,1]))+c(-buff[[a]],buff[[a]])
    ylim <- c(min(X[[a]][,2]),max(X[[a]][,2]))+c(-buff[[a]],buff[[a]])
    area[a] <- diff(xlim)*diff(ylim)
  }
  return(area)
}

sim.USCR.multisession <-
  function(N.session=NA,lambda.N=NA,lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,buff=NA,
           p0=NA,theta.thin=NA,K1D=NA,obstype="poisson"){
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
    #get areas
    area <- getArea(X=X,buff=buff)
    
    #simulate sessions one at a time
    N <- rpois(N.session,lambda.N)
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <-  sim.USCR(N=N[g],lam0=lam0[g],p0=p0[g],sigma=sigma[g],K=K[g],X=X[[g]],buff=buff[g],
                           obstype=obstype,theta.d=theta.d[g],K1D=K1D[[g]])
      data[[g]]$area <- area[g]
      data[[g]]$lambda.N <- lambda.N[g]
    }
    return(data)
  }