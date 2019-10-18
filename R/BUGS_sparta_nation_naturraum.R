cat("
  model{

  # State model
  for (i in 1:nsite){ 
    for (t in 1:nyear){
      z[i,t] ~ dbern(muZ[i,t]) 
      logit(muZ[i,t])<- raum.t.effect[raumS[i],t] + 
                        raum.effect[raumS[i]] +                        
                        eta[i]
    }
  }   
  
  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j]) 
    Py[j]<- z[site[j],year[j]]*p[j] 

    #detection model:
    logit(p[j]) <-  a.p[year[j]] + 
                    raum.p[raum[j]] +
                    phenol.p * yday[j] + 
                    phenol2.p * pow(yday[j], 2) +
                    effort.p * Effort[j] +
                    single.p * singleList[j]
  } 

    #Derived parameters - annual occupancy
    for (t in 1:nyear) {
      for(s in 1:nstate){
        psi.raum[s,t] <- mean(z[raumS[s],t])
      }
    } 

    #Priors 
    # State model priors
    
    #year effects
    #year 1
    for(i in 1:nraum){
      raum.t.effect[i,1] ~ dnorm(0, 0.001)
    }

    #other years
    for(i in 1:nraum){
      for(t in 2:nyear){
        raum.t.effect[i,t] ~ dnorm(raum.t.effect[i,t-1], tau.a)
      }
    }
    tau.a <- 1/(sd.a * sd.a)
    sd.a ~ dt(0, 1, 1)T(0,) 


    #state effects
    for (i in 1:nraum) {
      raum.effect[i] ~ dnorm(0, tau.raum)       
    } 
    tau.raum <- 1/(sd.raum * sd.raum) 
    sd.raum ~ dt(0, 1, 1)T(0,) 

    #site effects
    for (i in 1:nsite){
      eta[i] ~ dnorm(0, tau.eta)       
    } 
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
    
    #random slopes for phenology
    #for(s in 1:nstate){
    #  phenol.p[s] ~ dnorm(mu.phenol, tau.phenol)
    #  phenol2.p[s] ~ dnorm(mu.phenol2, tau.phenol2)
    #}

    #random phenology effects
    phenol.p ~ dnorm(0,0.01)
    phenol2.p ~ dnorm(0,0.01)
    #tau.phenol <- pow(sd.phenol,-2)
    #tau.phenol2 <- pow(sd.phenol2,-2)    
    #sd.phenol ~ dt(0, 1, 1)T(0,) 
    #sd.phenol2  ~ dt(0, 1, 1)T(0,) 

    #state effects
    for (i in 1:nraum) {
      raum.p[i] ~ dnorm(0, tau.raum.p)       
    } 
    tau.raum.p <- 1/(sd.raum.p * sd.raum.p) 
    sd.raum.p ~ dt(0, 1, 1)T(0,)

    #effort effects
    effort.p ~ dnorm(0, 0.01)
    single.p ~ dnorm(0, 0.01)

  }
    ",fill=TRUE,file="R/BUGS_sparta_nation_naturraum.txt")