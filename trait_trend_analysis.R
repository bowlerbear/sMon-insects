#linking traits to trends

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/sMon-insects/R/sparta_wrapper_functions.R')

#######################################################################

#get annual indicies
load("model-outputs/out_RpS_Expertise_trends.RData")
out2 <- combineModels(out)
outTrends <- getBUGSfits(out2,param="psi.fs") 

#######################################################################

#compare adult and juvenile trends

#######################################################################

fitTrends <- function(df,myspecies=NULL){

#get data for species
  dfS <- subset(df,Species==myspecies)
  dfS <- dfS[order(dfS$ParamNu),]
  
#combine for a bugs model
bugs.data<-list(
  Index = dfS$mean,
  SE = dfS$sd,
  n.year = nrow(dfS))


#write a bugs model to estimate trends
cat("
    model{
    
    #give indices with error
      for(t in 1:n.year){
      tau[t] <- 1/SE[t]^2
        Index[t] ~ dnorm(actPop[t],tau[t])
      }
    
    #fit traits model
    int ~ dnorm(0,0.01)
    trend ~ dnorm(0,0.01)
    error ~ dunif(0,10)

    #fit model
    for(t in 1:n.year){
      actPop[t] ~ dnorm(expPop[t],error)
      #linear predictor
      expPop[t] <- int + trend  * t  
    }
      
    
    }
    ",fill=TRUE,file="R/Bugs_trend.txt")

#run model
source('R/BUGS_misc_functions.R')

#specify parameters to monitor
params <- c("int","trend")

inits <- function(){
  list(int = rnorm(1,0,0.1), trend= rnorm(1,0,0.1))}

#run model
out1 <- jags(bugs.data, inits=inits, params, "R/Bugs_trend.txt", n.thin=nt,
             n.chains=3, n.burnin=1000,n.iter=5000,parallel=T)

}

######################################################################


#get traits data
load("altraits.RData")
#limit to those with complete cases


#####################################################################

#write a bugs model to test interaction between traits and trends
cat("
  model{

  #give indices with error
  for(i in 1:n.species){
    for(t in 1:n.year){
      actPop[i,t] ~ dnorm(Index[i,t],SE[i,t])
    }
  }

  #fit traits model
  for(i in 1:n.species){
    for(t in 1:n.year){
      actPop[i,t] ~ dnorm(expPop[i,t],error[i,t])
      #linear predictor
      expPop[i,t] <- int + 
                    year[t] + 
                    inprod(beta[],traitsDF[species[i],]) + 
                    year[t] * inprod(betaInt[],traitsDF[species[i],])      


   for(i in 1:n.covs){
      beta[i] ~ dnorm(0,0.1)
   }


   for(i in 1:n.covs){
      betaInt[i] ~ dnorm(0,0.1)
    }


    }
    ",fill=TRUE,file="R/BUGS_sparta.txt")

######################################################################
