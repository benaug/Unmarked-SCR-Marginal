e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.USCR.multisession <- function(data,inits=NA,M=NA,obstype="poisson"){
  N.session <- length(data)
  if(length(M)!=N.session)stop("Must supply an M for each session.")
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.USCR(data[[g]],inits.use,M=M[g])
  }
  
  J <- unlist(lapply(data,function(x){nrow(x$X)}))
  maxM <- max(M)
  maxJ <- max(J)
  s <- array(NA,dim=c(N.session,maxM,2))
  z <- matrix(NA,N.session,maxM)
  y.noID <- matrix(NA,N.session,maxJ)
  for(g in 1:N.session){
    y.noID[g,1:J[g]] <- data[[g]]$y.noID
  }
  K1D <- matrix(NA,N.session,max(J))
  xlim <- matrix(NA,N.session,2)
  ylim <- matrix(NA,N.session,2)
  for(g in 1:N.session){
    s[g,1:M[g],] <- init.session[[g]]$s
    z[g,1:M[g]] <- init.session[[g]]$z
    K1D[g,1:J[g]] <- init.session[[g]]$K1D
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
  }
  
  #put X in ragged array
  X.new <- array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],] <- data[[g]]$X
  }
  dummy.data <- matrix(0,N.session,maxM) #dummy data not used, doesn't really matter what the values are
  
  return(list(y.noID=y.noID,s.init=s,z.init=z,K1D=K1D,J=J,X=X.new,xlim=xlim,ylim=ylim))
  
}