e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.USCR.Dcov <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,lam0=NA,p0=NA,
           sigma=NA,theta.d=NA,K=NA,X=NA,xlim=NA,ylim=NA,
           obstype="poisson",theta.thin=NA,K1D=NA){
    library(abind)
    #get expected N
    cellArea <- res^2
    lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
    lambda.N <- sum(lambda.cell)
    #simulate realized N
    N <- rpois(1,lambda.N)
    
    #recreate some Dcov things so we can pass fewer arguments into this function
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    
    # simulate a population of activity centers
    pi.cell <- lambda.cell/sum(lambda.cell)
    s.cell <- sample(1:n.cells,N,prob=pi.cell,replace=TRUE)
    #distribute activity centers uniformly inside cells
    s <- matrix(NA,nrow=N,ncol=2)
    for(i in 1:N){
      tmp <- which(cells==s.cell[i],arr.ind=TRUE) #x and y number
      s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]+res/2])
      s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]+res/2])
    }
    D <- e2dist(s,X)
    J <- nrow(X)
    lamd <- lam0 * exp(-D * D/(2 * sigma * sigma))
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
                X=X,K=K,K1D=K1D,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N)
    return(out)
  }