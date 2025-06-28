init.USCR.Dcov.multisession <- function(data,inits=NA,M=NA){
  N.session <- length(data)
  if(length(M)!=N.session)stop("Must supply an M for each session.")
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.USCR.Dcov(data[[g]],inits.use,M=M[g])
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
  n.cells <- unlist(lapply(data,function(x){x$n.cells}))
  n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
  n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
  n.cells.max <- max(n.cells)
  n.cells.x.max <- max(n.cells.x)
  n.cells.y.max <- max(n.cells.y)
  res <- unlist(lapply(data,function(x){x$res}))
  cellArea <- res^2
  xlim <- matrix(NA,N.session,2)
  ylim <- matrix(NA,N.session,2)
  x.vals <- matrix(NA,N.session,n.cells.x.max)
  y.vals <- matrix(NA,N.session,n.cells.y.max)
  dSS <- array(NA,dim=c(N.session,n.cells.max,2))
  InSS <- array(0,dim=c(N.session,n.cells.max))
  D.cov <- array(NA,dim=c(N.session,n.cells.max))
  cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
  
  for(g in 1:N.session){
    s[g,1:M[g],] <- init.session[[g]]$s
    z[g,1:M[g]] <- init.session[[g]]$z
    K1D[g,1:J[g]] <- init.session[[g]]$K1D
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
    x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
    y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
    dSS[g,1:n.cells[g],] <- data[[g]]$dSS
    InSS[g,1:n.cells[g]] <- data[[g]]$InSS
    D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
    cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
  }
  
  #put X in ragged array
  X.new <- array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],] <- data[[g]]$X
  }
  dummy.data <- matrix(0,N.session,maxM) #dummy data not used, doesn't really matter what the values are
  
  return(list(y.noID=y.noID,
              s.init=s,z.init=z,K1D=K1D,J=J,X=X.new,
              res=res,cellArea=cellArea,x.vals=x.vals,y.vals=y.vals,xlim=xlim,ylim=ylim,
              dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data))
  
}