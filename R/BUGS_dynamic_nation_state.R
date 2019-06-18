cat("
  model{

  # State model
  
  #first year/initial occupancy, model
  for (i in 1:nsite){ 

    #data distribution
    z[i,1] ~ dbern(muZ[i,1])

    #model
    logit(muZ[i,1]) <- eta[i] 
  }

  #subsequent years
  for (i in 1:nsite){ 
    for (t in 2:nyear){

      #data distribution
      z[i,t] ~ dbern(muZ[i,t]) 

      #overall model
      muZ[i,t] <- persist[i,t]*muZ[i,t-1] + colonization[i,t]*(1-muZ[i,t-1])

      #persistence occupancy model
      logit(persist[i,t]) <- a.persist[t] 

      #colonization occupancy model
      logit(colonization[i,t]) <- a.colonize[t] 
    }
  }   
 
  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j]) #data is Y
    Py[j]<- z[site[j],year[j]]*p[j] 

    #detection model:
    logit(p[j]) <-  a.p[year[j]] + 
                    phenol.p * yday[j] + 
                    phenol2.p * pow(yday[j], 2) + 
                    effort.p * Effort[j] +
                    single.p * singleList[j]
    } 
  
  # Derived parameters - annual occupancy
  for (t in 1:nyear) {  
    psi.fs[t] <- sum(z[1:nsite, t])/nsite
  } 

  #Priors 

  # State model priors
    
    #year effects on persistence
    for (t in 1:nyear) {
      a.persist[t] ~ dnorm(a.mu.persist, a.tau.persist)            
    }
    a.mu.persist ~ dnorm(0, 0.01)
    a.tau.persist <- 1 / (a.sd.persist * a.sd.persist)                 
    a.sd.persist ~ dt(0, 1, 1)T(0,)

    #year effects on colonization
    for (t in 1:nyear) {
    a.colonize[t] ~ dnorm(a.mu.colonize, a.tau.colonize)            
    }
    a.mu.colonize ~ dnorm(0, 0.01)
    a.tau.colonize <- 1 / (a.sd.colonize * a.sd.colonize)                 
    a.sd.colonize ~ dt(0, 1, 1)T(0,)

    #site effects
    for (i in 1:nsite) {
      eta[i] ~ dnorm(mu.eta, tau.eta)       
    } 
    mu.eta ~ dnorm(0,0.01)
    tau.eta <- 1/(sd.eta * sd.eta) 
    sd.eta ~ dt(0, 1, 1)T(0,) 

    #Observation model priors 
    #year effects
    for (t in 1:nyear) {
    a.p[t] ~ dnorm(a.mu.p, a.tau.p)            
    }
    a.mu.p ~ dnorm(0, 0.01)
    a.tau.p <- 1 / (a.sd.p * a.sd.p)                 
    a.sd.p ~ dt(0, 1, 1)T(0,) 

    #observation model covariates
    #int.p.temp ~ dunif(0,1)
    #int.p <- logit(int.p.temp)
    phenol.p ~ dnorm(0, 0.01)
    phenol2.p ~ dnorm(0, 0.01)
    effort.p ~ dnorm(0, 0.01)
    single.p ~ dnorm(0, 0.01)
  }
    ",fill=TRUE,file="R/BUGS_dynamic.txt")