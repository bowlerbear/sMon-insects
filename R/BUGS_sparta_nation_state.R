cat("
  model{

  # State model
  for (i in 1:nsite){ 
    for (t in 1:nyear){
      z[i,t] ~ dbern(muZ[i,t]) 
      logit(muZ[i,t])<- int +  
                        eta[i]+
                        state.a.effect[stateS[i]] +
                        state.t.effect[stateS[i],t] 
    }
  }   
  
  ### Observation Model
  for(j in 1:nvisit) {
    y[j] ~ dbern(Py[j]) 
    Py[j]<- z[site[j],year[j]]*p[j] 

    #detection model:
    logit(p[j]) <-  a.p[year[j]] + 
                    state.p[state[j]] +
                    phenol.p[state[j]] * yday[j] + 
                    phenol2.p[state[j]] * pow(yday[j], 2) +
                    effort.p * Effort[j] +
                    single.p * singleList[j]
  } 

    #Derived parameters - annual occupancy
    for (t in 1:nyear) {
      for(s in 1:nstate){
        psi.fs[s,t] <- mean(z[stateS[s],t])
      }
    } 

    #Priors 
    # State model priors
    int ~ dnorm(0, 0.001)

    #box effects
    for (i in 1:nbox) {
      box.effect[i] ~ dnorm(0, tau.box)       
    } 
    tau.box <- 1/(sd.box * sd.box) 
    sd.box ~ dt(0, 1, 1)T(0,) 

    #state effects
    for (i in 1:nstate) {
      state.a.effect[i] ~ dnorm(0, tau.state)       
    } 
    tau.state <- 1/(sd.state * sd.state) 
    sd.state ~ dt(0, 1, 1)T(0,) 

    #state time effects
    for (i in 1:nstate) {
      for(t in 1:nyear){
        state.t.effect[i,t] ~ dnorm(0, tau.t.state)  
      }
    } 
    tau.t.state <- 1/(sd.t.state * sd.t.state) 
    sd.t.state ~ dt(0, 1, 1)T(0,) 


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
    for(s in 1:nstate){
      phenol.p[s] ~ dnorm(mu.phenol, tau.phenol)
      phenol2.p[s] ~ dnorm(mu.phenol2, tau.phenol2)
    }

    #random phenology effects
    mu.phenol ~ dnorm(0,0.001)
    mu.phenol2 ~ dnorm(0,0.001)
    tau.phenol <- pow(sd.phenol,-2)
    tau.phenol2 <- pow(sd.phenol2,-2)    
    sd.phenol ~ dt(0, 1, 1)T(0,) 
    sd.phenol2  ~ dt(0, 1, 1)T(0,) 

    #state effects
    for (i in 1:nstate) {
      state.p[i] ~ dnorm(0, tau.state.p)       
    } 
    tau.state.p <- 1/(sd.state.p * sd.state.p) 
    sd.state.p ~ dt(0, 1, 1)T(0,)

    #effort effects
    effort.p ~ dnorm(0, 0.001)
    single.p ~ dnorm(0, 0.001)

  }
    ",fill=TRUE,file="R/BUGS_sparta_nation_state.txt")