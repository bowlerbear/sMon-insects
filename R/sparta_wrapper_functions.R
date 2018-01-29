#given the directory in which the sparta model outputs are found, put all the models in a list
getSpartaModels<-function(dir){
  list.files(dir) -> sp_mods
  
  models <- lapply(sp_mods, function(sp) {
    load(file=paste0(paste0(dir,"/"),sp))
    return(out)
  }
  )
  names(models) <- gsub(sp_mods, pa="\\.rdata", repl="")
  
  return(models)
}

#get annual predictions of occupancy for each species
annualPredictions <- function(models){
  library(sparta)
  library(plyr)
  ldply(models,function(x){
  
  #get annual predictions
  temp <- data.frame(summary(x))
  temp$Year<-as.numeric(row.names(summary(x)))
  temp$Species <- x$SPP_NAME
  
  #get RHat values
  bugsOutput <- x$BUGSoutput$summary
  bugsOutput <- data.frame(bugsOutput[grepl("psi.fs",row.names(bugsOutput)),])
  temp$Rhat <- as.numeric(bugsOutput$Rhat)
  
  return(temp)
})
  
}

#plot these predictions (restrict to species with more than 50 observations)
plotPredictions <- function(myAnnualPredictions,rawdata){
  
require(ggplot2)
ggplot(data=subset(myAnnualPredictions,Species %in% names(table(rawdata$Species))[table(rawdata$Species)>50])) +
  geom_line(aes(x = Year, mean))+
  geom_point(aes(x = Year, mean,colour = factor(Rhat<1.1)))+
  geom_ribbon(aes(x=Year, ymin = quant_025, ymax = quant_975), alpha=0.50)+
  theme_bw() +
  scale_x_continuous(labels=c(1990,1995,2000,2005,2010))+
  facet_wrap( ~ Species) +
  theme(legend.position = "none")
}


#get estimates of each species linear population trends and 95% CI
calculateTrends<-function(models){

  #get trends for each model
  trends <- lapply(models, occurrenceChange, firstYear=min(df$Year), lastYear=max(df$Year))
  names(trends) <- sapply(models,function(x) x$SPP_NAME)
  
  #convert into a data frame
  outputs <- data.frame(
    mean.trend = sapply(trends, function(x) x$mean),
    CI.lower = sapply(trends, function(x) x$CIs[1]),
    CI.upper = sapply(trends, function(x) x$CIs[2])) 
  
  #add number of total records for each species
  outputs$nuRecords <- sapply(models,function(x) x$species_observations)
  
  #return it
  return(outputs)
}