library(nimble)
library(coda)
source("sim.USCR.multisession.R")
source("sim.USCR.R")
source("init.USCR.multisession.R")
source("init.USCR.R")
source("NimbleModel USCR Multisession Poisson Marginal.R")
source("NimbleFunctions USCR Multisession Poisson Marginal.R")
source("sSampler Poisson Multisession Marginal.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#Here, I'll simulate 3 populations with different K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session <- 3
D  <-  rep(0.2,N.session) #expected density in units of sigma and X
lam0 <- rep(0.25,N.session)
sigma <- rep(0.5,N.session) #change prior if you change sigma, set up with informative prior around 0.5
obstype <- "poisson"
K <- c(5,6,7) #number of occasions
buff <- rep(3,N.session) #state space buffer
#make trapping arrays
X1 <- expand.grid(3:11,3:11)
X2 <- expand.grid(3:12,3:12)
X3 <- expand.grid(3:13,3:13)
X <- list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area <- getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda.N <- D*area
lambda.N #expected N in each session

#Simulate some data
data <- sim.USCR.multisession(N.session=N.session,lambda.N=lambda.N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
g <- 1 #session to look at
str(data[[g]]$y.noID) #trap by occasion counts of unidentified detections

#To fit models faster, we will sum data over occasions
for(g in 1:N.session){
  data[[g]]$y.noID <- rowSums(data[[g]]$y.noID)
}

##Fit model in Nimble##
M <- c(150,150,150) #data augmentation level
X <- sapply(data,function(x){x$X})
J <- sapply(X,nrow) #number of detectors
area <- sapply(data,function(x){x$area}) #pull areas out of 

inits <- list(lam0=rep(1,N.session),sigma=rep(0.5,N.session)) #ballpark inits to build data, one per session.
nimbuild <- init.USCR.multisession(data=data,M=M,inits=inits)

#inits for nimble
N.init <- rowSums(nimbuild$z,na.rm=TRUE) #N.init must be consistent with z.init!
Niminits <- list(N=N.init,z=nimbuild$z,D=mean(N.init/area), #initialize D from N.init for faster convergence
                 s=nimbuild$s,lam0.fixed=0.5,sigma.fixed=0.5)

#constants for Nimble
constants <- list(N.session=N.session,M=M,J=J,K1D=nimbuild$K1D,xlim=nimbuild$xlim,ylim=nimbuild$ylim,area=area)

#supply data to nimble
Nimdata <- list(X=nimbuild$X,y.noID=nimbuild$y.noID)

# set parameters to monitor
parameters <- c('lam0.fixed','sigma.fixed','N','D','lambda.N')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c('lam0.fixed','sigma.fixed','D')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy=FALSE,
                      nodes=config.nodes) 

#Add z/N sampler
z.ups <- round(M*0.25) # how many z proposals per iteration per session?
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[",g,",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  bigLam.nodes <- Rmodel$expandNodeNames(paste("bigLam[",g,",",1:J[g],"]"))#only need this in calcNodes
  lam.noID.nodes <- Rmodel$expandNodeNames(paste("lam.noID[",g,",",1:J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.noID.nodes)
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.noID.nodes=y.noID.nodes,lam.nodes=lam.nodes,
                                                   lam.noID.nodes=lam.noID.nodes,N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

#add activty center sampler
for(g in 1:N.session){
  for(i in 1:M[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,J=J[g],xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#add block update for lam0, sigma, and/or D if posteriors correlated
conf$addSampler(target = c("lam0.fixed",'sigma.fixed',"D"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)

burnin <- 500
plot(coda::mcmc(mvSamples[-c(1:burnin),]))

#reminder of what some targets are
sapply(data,function(x){x$N})
sapply(data,function(x){x$lambda.N})

cor(mvSamples[-c(1:burnin),])
