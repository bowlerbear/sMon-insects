#linking traits to trends

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/sMon-insects/R/sparta_wrapper_functions.R')

#######################################################################

#get annual indicies

trendFiles <- list.files("derived-data/")[grepl("outTrends",list.files("derived-data/"))]
allTrends <- ldply(trendFiles, function(x){
  load(paste("derived-data",x,sep="/"))
  
  #get file metadata
  outTrends$stage <- strsplit(x,"_")[[1]][3]
  outTrends$state <- strsplit(x,"_")[[1]][4]
  outTrends$state <- strsplit(outTrends$state ,"\\.")[[1]][1]
  
  #return data frame
  return(outTrends)
})

#######################################################################

#compare adult and juvenile trends

ggplot(subset(allTrends,state=="Sa"),
       aes(x=ParamNu,y=mean,colour=stage))+
  geom_line()+
  facet_wrap(~Species)+
  theme(legend.position="none")

ggplot(subset(allTrends,state=="SH"),
       aes(x=ParamNu,y=mean,colour=stage))+
  geom_line()+
  facet_wrap(~Species)+
  theme(legend.position="none")


#######################################################################

fitTrends <- function(df){

#get data for species
  dfS <- df
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

#return predictions
return(data.frame(trend=out1$mean$trend,trend_sd=out1$sd$trend))

}

#####################################################################

#fit to each stage and state

trendEstimates <- ddply(allTrends,.(Species,state,stage),
                        function(x){
                          fitTrends(x)
                          })

######################################################################

#get traits data
load("alltraits.RData")
#limit to those with complete cases?

#########################################################################

#plot of traits vs trends
trendEstimate <- merge(trendEstimates,alltraits,by="Species")

#temperature preference plot
ggplot(trendEstimate,aes(x=TMean,y=trend))+
  geom_point()+
  facet_wrap(state~stage)

ggplot(trendEstimate,aes(x=TMean,y=trend))+
  geom_point(aes(colour=state,shape=stage))+
  stat_smooth(method="gam")

#habitat preference
ggplot(trendEstimate,aes(x=Habitat.y,y=trend))+
  geom_boxplot()+
  facet_wrap(state~stage)

#range size
ggplot(trendEstimate,aes(x=nuEuroGrids,y=trend))+
  geom_point()+
  facet_wrap(state~stage)

#generation time
ggplot(trendEstimate,aes(x=development.time..mean.,y=trend))+
  geom_point()+
  facet_wrap(state~stage)

ggplot(trendEstimate,aes(x=development.time..mean.,y=trend))+
  geom_point(aes(colour=state,shape=stage))

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
