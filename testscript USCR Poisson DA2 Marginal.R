library(nimble)
library(coda)
source("sim.USCR.R")
source("init.USCR.R")
source("NimbleModel USCR Poisson DA2 Marginal.R")
source("NimbleFunctions USCR Poisson DA2 Marginal.R")
source("sSampler Poisson Marginal.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
N <- 38
lam0 <- 0.25
sigma <- 0.5 #change prior if you change sigma, set up with informative prior around 0.5
K <- 10
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11)
obstype <- "poisson"

#Simulate some data
data <- sim.USCR(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#What is the observed data?
str(data$y.noID) #trap by occasion counts of unidentified detections

#To fit models faster, we will sum data over occasoins
data$y.noID <- rowSums(data$y.noID)

#Data augmentation level
M <- 200
X <- data$X
J <- nrow(X) #number of detectors

inits <- list(lam0=1,sigma=1) #ballpark inits to build data
nimbuild <- init.USCR(data=data,M=M,inits=inits)

#inits for nimble
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z), #z and N inits must be consistent
                 lambda.N=sum(nimbuild$z), #converges faster if you set N and lambda.N at similar values
                 s=nimbuild$s,sigma=inits$sigma,lam0=inits$lam0)

#constants for Nimble
constants <- list(M=M,J=J,xlim=nimbuild$xlim,ylim=nimbuild$ylim)

#supply data to nimble
Nimdata <- list(X=as.matrix(X),K1D=data$K1D,y.noID=data$y.noID)

#set parameters to monitor
parameters <- c('lambda.N','lam0','sigma','N')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
config.nodes <- c("lambda.N","lam0","sigma")
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = FALSE) 

###*required* sampler replacement
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal.
#nodes used for update
y.noID.nodes <- Rmodel$expandNodeNames(paste("y.noID[1:",J,"]"))
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M,",1:",J,"]"))
bigLam.nodes <- Rmodel$getDependencies("bigLam") #only need this in calcNodes
lam.noID.nodes <- Rmodel$getDependencies("lam.noID")
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,lam.nodes,bigLam.nodes,lam.noID.nodes,y.noID.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 lam.nodes=lam.nodes,lam.noID.nodes=lam.noID.nodes,
                                                 y.noID.nodes=y.noID.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)

#must use this activity center sampler
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,J=J,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#use block update for correlated posteriors. Can use "tries" to control how many times per iteration
conf$addSampler(target = c("lam0","sigma","lambda.N"),
                type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

mvSamples <- as.matrix(Cmcmc$mvSamples)
burnin <- 250
plot(mcmc(mvSamples[-c(1:burnin),]))

#posterior correlation
tmp <- cor(mcmc(mvSamples[-c(1:burnin),]))
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)
