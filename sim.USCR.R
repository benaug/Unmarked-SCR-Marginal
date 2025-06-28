e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.USCR <-
  function(N=NA,lam0=NA,p0=NA,sigma=NA,theta.d=NA,K=NA,X=NA,buff=NA,obstype="poisson",K1D=NA){
    library(abind)
    X <- as.matrix(X)
    xlim <- c(min(X[,1]),max(X[,1]))+c(-buff,buff)
    ylim <- c(min(X[,2]),max(X[,2]))+c(-buff,buff)
    s <- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D <- e2dist(s,X)
    lamd <- lam0 * exp(-D * D/(2 * sigma * sigma))
    J <- nrow(X)
    
    #trap operation
    if(!any(is.na(K1D))){
      if(any(K1D>K)){
        stop("Some entries in K1D are greater than K.")
      }
      if(is.null(dim(K1D))){
        if(length(K1D)!=J){
          stop("K1D vector must be of length J.")
        }
      }
    }else{
      K1D <- rep(K,J)
    }
    
    #Capture individuals
    y.true <- array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(lam0))stop("must provide lam0 for bernoulli obstype")
      pd <- 1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,1:K1D[j]] <- rbinom(K1D[j],1,pd[i,j])
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,1:K1D[j]] <- rpois(K1D[j],lamd[i,j])
        }
      }
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      for(i in 1:N){
        for(j in 1:J){
          y.true[i,j,1:K1D[j]] <- rnbinom(K1D[j],mu=lamd[i,j],size=theta.d)
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    y.true.original <- y.true
    y.true <- y.true[which(rowSums(y.true)>0),,]
    n.cap <- sum(rowSums(y.true)>0)
    y.noID <- apply(y.true,c(2,3),sum)
    
    out <- list(y.noID=y.noID, #observed data
                y.true=y.true.original,n.cap=n.cap,s=s, #true data
                X=X,K=K,K1D=K1D,buff=buff,s=s,xlim=xlim,ylim=ylim,N=N)
    return(out)
  }