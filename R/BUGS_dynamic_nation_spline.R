cat("
  model{
    
    # State model
    
    for (i in 1:nsite){ 
    
    #first year/initial occupancy, model
    
    #data distribution
    z[i,1] ~ dbern(muZ[i,1])
    
    logit(muZ[i,1]) <- logit.mu + eta[i] + box.a[boxS[i]]
    
    for (t in 2:nyear){
    
      #data distribution
      z[i,t] ~ dbern(muZ[i,t]) 
    
      #overall model
      muZ[i,t] <- persist[i,t]*muZ[i,t-1] + colonize[i,t]*(1-muZ[i,t-1])
    
      #persistence occupancy model
      logit(persist[i,t]) <- muphi[t] + eta.persist[i] + lphi[boxS[i]]
    
      #colonization occupancy model
      logit(colonize[i,t]) <- mugam[t] + eta.colonize[i] +lgam[boxS[i]]
      }
    }   
    
    ### Observation Model
    for(j in 1:nvisit) {
      y[j] ~ dbern(Py[j]) #data is Y
      Py[j]<- z[site[j],year[j]]*p[j] 
    
      #detection model:
      logit(p[j]) <-  mup[year[j]] +
      lp[state[j]] +
      mu.phenol * yday[j] + 
      mu.phenol2 * pow(yday[j], 2) +
      phenol.s[box[j]] * yday[j] + 
      phenol2.s[box[j]] * pow(yday[j],2)+
      effort.p * Effort[j] +
      single.p * singleList[j]
    } 
    
  # Derived parameters

    #metapopulation-level
    for(t in 1:nyear){
      psi.fs[t] <- sum(z[1:nsite,t])/nsite 
      psi2[t] <- mean(z[,t])
    }
    
    #yearly-growth
    for(t in 2:nyear){  
      growthr[t] <- psi.fs[t]/psi.fs[t-1] 
    } 
    
    #overall
    mean.growth <- sum(growthr[2:nyear])/(nyear-1)
    mean.p <- mean(p)#mean detection probability
    
    #state-level
    for(s in 1:nsbox){
      for(t in 1:nyear){
        psi.state[s,t] <- mean(z[stateS[s],t]) 
      }
    }
    
    #yearly-growth state
    for(s in 1:nbox){
      for(t in 2:nyear){
        growthr.state[s,t] <- psi.state[s,t]/psi.state[s,t-1] 
      }
    }
    
    #overall state
    for(s in 1:nbox){
      mean.growth.state[s] <- mean(growthr.state[state[s],])
    }
    
    
    #Priors 
    
    # State model priors
    
    #spatial effects
    
    #intercept
    mu ~ dunif(0,1)
    logit.mu <- logit(mu)
    
    #state.a random effects
    for (s in 1:nstate) {
      state.a[s] ~ dnorm(0, tau.state.a)       
    } 
    tau.state.a ~ dgamma(.1,.1) 
    
    #temporal intercepts
    
    #detection model 
    for(i in 1:nyear){
      mup.prob[i] ~ dunif(0,1)
      mup[i] <- logit(mup.prob[i])
    }
    
    #persist/colon
    for(i in 2:nyear){
      muphi.prob[i] ~ dunif(0,1)
      muphi[i] <- logit(muphi.prob[i])
    
      mugam.prob[i] ~ dunif(0,1)
      mugam[i] <- logit(mugam.prob[i])
    }
    
    
    #random site variation
    
    tau.phi ~ dgamma(.1,.1)
    tau.gam ~ dgamma(.1,.1)
    taup ~ dgamma(.1,.1)
    
    for(s in 1:nstate){
      lphi[s] ~ dnorm(0,tau.phi)
      lgam[s]  ~ dnorm(0,tau.gam)
      lp[s]  ~ dnorm(0,taup)
    }
    
    #Observation model priors
    #intercept - see above
    #site-level variaion - see above
    
    #effort effects
    effort.p ~ dnorm(0, 0.01)
    single.p ~ dnorm(0, 0.01)
    
    #random slopes for phenology
    for(s in 1:nstate){
      phenol.s[s] ~ dnorm(0, tau.phenol)
      phenol2.s[s] ~ dnorm(0, tau.phenol2)
    }
    mu.phenol ~ dnorm(0,0.01)
    mu.phenol2 ~ dnorm(0,0.01)
    tau.phenol ~ dgamma(.1,.1)
    tau.phenol2 ~ dgamma(.1,.1) 


    #splines
    #spline model for initial occupancy
    eta <- X %*% b
    K1 <- S1[1:nspline,1:nspline] * lambda[1]  + S1[1:nspline,nspline1:nspline2] * lambda[2]
    b[1:nspline1] ~ dmnorm(zero[1:nspline1],K1) 
    
    ## smoothing parameter priors 
    for (i in 1:2) {
      lambda[i] ~ dgamma(.05,.005)
      rho[i] <- log(lambda[i])
    }

    #spline model for persistence
    eta.persist <- X %*% b.p
    K1.p <- S1[1:nspline,1:nspline] * lambda.p[1]  + S1[1:nspline,nspline1:nspline2] * lambda.p[2]
    b.p[1:nspline] ~ dmnorm(zero[1:nspline],K1.p)  
    
    ## smoothing parameter priors
    for (i in 1:2) {
      lambda.p[i] ~ dgamma(.05,.005)
      rho.p[i] <- log(lambda.p[i])
    }

    #spline model for colonization
    eta.colonize <- X %*% b.c
    K1.c <- S1[1:nspline,1:nspline] * lambda.c[1]  + S1[1:nspline,nspline1:nspline2] * lambda.c[2]
    b.c[1:nspline] ~ dmnorm(zero[1:nspline],K1.c)  
    
    ## smoothing parameter priors 
    for (i in 1:2) {
      lambda.c[i] ~ dgamma(.05,.005)
      rho.c[i] <- log(lambda.c[i])
    }

    
    }
    ",fill=TRUE,file="R/BUGS_dynamic_nation_spline.txt")


  
    


   