e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.USCR.Dcov <- function(data,inits=NA,M=NA,obstype="poisson"){
  y.noID <- data$y.noID
  X <- as.matrix(data$X)
  J <- nrow(X)
  buff <- data$buff
  K1D <- data$K1D
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  lam0 <- inits$lam0
  sigma <- inits$sigma
  
  #Augment data
  y.true <- matrix(0,M,J)
  
  #make a plausibly true complete data set to initialize data
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  for(j in 1:J){
    prob <- lamd[,j]
    prob <- prob/sum(prob)
    y.true[,j] <- y.true[,j] + rmultinom(1,y.noID[j],prob=prob)
  }
  
  #initialize z
  z.init <- 1*(rowSums(y.true)>0)
  
  #update s.init for individuals assigned samples
  idx <- which(rowSums(y.true)>0) #switch for those actually caught
  for(i in idx){
    trps <- matrix(X[y.true[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  
  #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
  e2dist  <-  function (x, y){
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  getCell  <-  function(s,res,cells){
    cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
  }
  alldists <- e2dist(s.init,data$dSS)
  alldists[,data$InSS==0] <- Inf
  for(i in 1:M){
    this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
    if(data$InSS[this.cell]==0){
      cands <- alldists[i,]
      new.cell <- which(alldists[i,]==min(alldists[i,]))
      s.init[i,] <- data$dSS[new.cell,]
    }
  }
  
  #check starting logProb
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  ll.y <- y.true*0
  for(i in 1:M){
    if(z.init[i]==1){
      ll.y[i,] <- dpois(y.true[i,],K1D*lamd[i,],log=TRUE)
    }
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K1D (if supplied by user) or problem initializing data.")

  return(list(s=s.init,z=z.init,K1D=K1D))
}