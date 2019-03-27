model{
  #Define prior distributions for community-level model parameters
  
  omega ~ dunif(0,1)
  
  # Mean and precision (tau) for community level distributions for species-level 
  # random effects
  
  #occupancy model
  #population effects
  mu.b0 ~ dnorm(0,0.001)
  mu.b1 ~ dnorm(0,0.001)
  mu.b2 ~ dnorm(0,0.001)
  #random species effects of year and site
  tau.b0 ~ dgamma(0.1,0.1)
  tau.b1 ~ dgamma(0.1,0.1)
  tau.b2 ~ dgamma(0.1,0.1)
  
  #detection model
  #population effect
  mu.a0 ~ dnorm(0,0.001)
  mu.a1 ~ dnorm(0,0.001)
  #random effect
  tau.a0 ~ dgamma(0.1,0.1)
  tau.a1 ~ dgamma(0.1,0.1)
  

  #Set species loop
  for (i in 1:(n+nzeroes)) {
    
    #Create priors for species i from the community level prior distributions
    
    w[i] ~ dbern(omega)
    
    #occupancy model coefficients
    b0[i] ~ dnorm(mu.b0, tau.b0)
    b1[i] ~ dnorm(mu.b1, tau.b1)
    b2[i] ~ dnorm(mu.b2, tau.b2)

    #detection model coefficients
    a0[i] ~ dnorm(mu.a0, tau.a0)
    a1[i] ~ dnorm(mu.a1, tau.a1)
    
    #Create a loop to estimate the Z matrix 
    #true occurrence for species i at point j
    
    #OCCUPANCY MODEL:
    
    #Site loop
    for (j in 1:J) {
      #predicted occurance
      logit(psi[j,i]) <- b0[i] + b1[i]*Year[j] + b2[i]*Site[j] 
      
      #predicted occurrence x indicator
      mu.psi[j,i] <- psi[j,i]*w[i]
      
      #true occurrence
      Z[j,i] ~ dbern(mu.psi[j,i])
      
      #Create a loop to estimate detection for species i at point j during
      #sampling period k. (k is a visit)
      #DETECTION MODEL
      
      #Visit loop
      for (k in 1:K[j]) {
        logit(p[j,k,i]) <-  a0[i] + a1[i]*dates[j,k]
        
        mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
        
        #X is the data
        X[j,k,i] ~ dbern(mu.p[j,k,i])
        
        Xnew[j,k,i] ~ dbin(mu.p[j,k,i], 1)
        
        #Create simulated dataset to calculate the Bayesian p-value
        d[j,k,i]<-  abs(X[j,k,i] - mu.p[j,k,i]) 
        dnew[j,k,i]<- abs(Xnew[j,k,i] - mu.p[j,k,i]) 
        d2[j,k,i]<- pow(d[j,k,i],2)  
        dnew2[j,k,i]<- pow(dnew[j,k,i],2) 
        
      }
      dsum[j,i]<- sum(d2[j,1:K[j],i]) 
      dnewsum[j,i]<- sum(dnew2[j,1:K[j],i])
      
    }}
  
  #Calculate the discrepancy measure, which is then defined as the mean(p.fit > p.fitnew)
  p.fit<-sum(dsum[1:J,1:(n+nzeroes)]) 
  p.fitnew<-sum(dnewsum[1:J,1:n+nzeroes])
  
  
  #Sum all species observed (n) and unobserved species (n0) to find the
  #total estimated richness
  n0 <- sum(w[(n+1):(n+nzeroes)])
  N <- n + n0
  
  # Create a loop to determine point level richness estimates for the 
  # whole community and for subsets or assemblages of interest.
  for(j in 1:J){
    N_site[j]<- inprod(Z[j,1:n],w[1:n])
  }
  
  #Estimation of Jaccard similarity index between sites (K is number of unique sites, K=41; M1 
  #and M2 are vectors of indices of paired 
  #sites).
  for(j in 1:J2){
    C[j]<-(inprod(Z[M1[j],1:n],Z[M2[j],1:n]))/ (sum(Z[M1[j],1:n])+sum(Z[M2[j],1:n])-inprod(Z[M1[j],1:n],Z[M2[j],1:n]))
  }
}
