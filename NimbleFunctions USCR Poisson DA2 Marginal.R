GetbigLam <- nimbleFunction(
  run = function(lam = double(2), z = double(1)){ 
    returnType(double(1))
    M <- nimDim(lam)[1]
    J <- nimDim(lam)[2]
    bigLam <- rep(0,J)
    for(i in 1:M){
      if(z[i]==1){
        bigLam <- bigLam + lam[i,]
      }
    }
    return(bigLam)
  }
)

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    J <- control$J
    z.ups <- control$z.ups
    y.noID.nodes <- control$y.noID.nodes
    lam.nodes <- control$lam.nodes
    lam.noID.nodes <- control$lam.noID.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    #track these "manually" so computations faster than nimble will do them
    bigLam.initial <- model$bigLam
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off all individuals
        if(model$N[1]==1){ #is this the last individual?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y.noID <- model$getLogProb(y.noID.nodes)
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn off
          bigLam.proposed <- bigLam.initial - model$lam[pick,] #subtract these out before calculate
          #make sure you didn't end up with any negative numbers due to machine precision
          bigLam.proposed[bigLam.proposed<0] <- 0
          model$calculate(lam.nodes[pick])
          model$bigLam <<- bigLam.proposed
          model$calculate(lam.noID.nodes)
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y.noID <- model$calculate(y.noID.nodes)
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.noID) -
            (lp.initial.N + lp.initial.y.noID)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["bigLam",1][1:J] <<- model[["bigLam"]][1:J]
            mvSaved["lam.noID",1][1:J] <<- model[["lam.noID"]][1:J]
            bigLam.initial <- bigLam.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["bigLam"]][1:J] <<- mvSaved["bigLam",1][1:J]
            model[["lam.noID"]][1:J] <<- mvSaved["lam.noID",1][1:J]
            model$calculate(y.noID.nodes)
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y.noID <- model$getLogProb(y.noID.nodes) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn on
          model$calculate(lam.nodes[pick])
          bigLam.proposed <- bigLam.initial + model$lam[pick,] #add these after calculate
          model$bigLam <<- bigLam.proposed
          model$calculate(lam.noID.nodes)
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y.noID <- model$calculate(y.noID.nodes)
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.noID) -
            (lp.initial.N + lp.initial.y.noID)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["bigLam",1][1:J] <<- model[["bigLam"]][1:J]
            mvSaved["lam.noID",1][1:J] <<- model[["lam.noID"]][1:J]
            bigLam.initial <- bigLam.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["bigLam"]][1:J] <<- mvSaved["bigLam",1][1:J]
            model[["lam.noID"]][1:J] <<- mvSaved["lam.noID",1][1:J]
            model$calculate(y.noID.nodes)
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)