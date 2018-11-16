#given the directory in which the sparta model outputs are found, put all the models in a list
#@param dir = directory of sparta model outputs

getSpartaModels<-function(dir){
  
  if(!dir%in%list.files()) stop('no such directory found')  
  
  sp_mods <- list.files(dir)[grepl(".rdata",list.files(dir))] 
  
  models <- lapply(sp_mods, function(sp) {
    load(file=paste0(paste0(dir,"/"),sp))
    return(out)
  }
  )
  
  names(models) <- gsub(sp_mods, pa="\\.rdata", repl="")
  
  return(models)
}

#get annual predictions of occupancy for each species
#@param models = output of getSpartaModels()

annualPredictions <- function(models){
 
  require(sparta)
  require(plyr)
  
  out<-ldply(models,function(x){
  
  #get bugs output
  temp <- x$BUGSoutput$summary
  temp <- temp[grep("psi.fs\\[",row.names(temp)),]
  temp <- data.frame(temp)
  names(temp) <- gsub("X2.5.","quant_025", names(temp)) 
  names(temp) <- gsub("X97.5.","quant_975", names(temp)) 

  #add on year and species info
  temp$Year <- x$min_year:x$max_year
  temp$Species <- x$SPP_NAME  
    
  #get records per species
  temp$nuRecords <- x$species_observations
  temp$nuSites <- x$nsites
  
  #reorganise
  temp<-temp[,c("Species","Year","mean","quant_025","quant_975","Rhat","n.eff","nuRecords","nuSites")]
  
  return(temp)
})

  out<-out[,-1]#get rid of superfluous first column
  return(out)

}

#plot these predictions (restrict to species with more than 50 observations)
#@param myAnnualPredictions = output from annualPredictions()

plotPredictions <- function(myAnnualPredictions){

#decide on year breaks to nearest decade
surveyYears<-sort(unique(myAnnualPredictions$Year))
surveyYears<-surveyYears[surveyYears%%10==0]

require(ggplot2)
  
#plot
ggplot(data=subset(myAnnualPredictions,nuRecords>50)) +
  geom_line(aes(x = Year, mean))+
  geom_point(aes(x = Year, mean,colour = factor(Rhat<1.1)))+
  geom_ribbon(aes(x=Year, ymin = quant_025, ymax = quant_975), alpha=0.50)+
  theme_bw() +
  scale_x_continuous(breaks = surveyYears, labels = surveyYears)+
  facet_wrap( ~ Species) +
  theme(legend.position = "none") +
  ylab("Predicted occupancy proportion")
}


#get estimates of each species linear population trends and 95% CI
#@param models = output of getSpartaModels()
calculateTrends<-function(models){

  #get trends for each model
  trends <- lapply(models, function(x) occurrenceChange(bayesOut= x, firstYear=x$min_year, lastYear=x$max_year))
  names(trends) <- sapply(models,function(x) x$SPP_NAME)
  
  #convert into a data frame
  outputs <- data.frame(
    mean.trend = sapply(trends, function(x) x$mean),
    CI.lower = sapply(trends, function(x) x$CIs[1]),
    CI.upper = sapply(trends, function(x) x$CIs[2])) 
  
  #add number of total records for each species
  outputs$nuRecords <- sapply(models,function(x) x$species_observations)
  
  # indicates wheter the trend we see is 'significant' or not
  outputs$sig<- ifelse(outputs$CI.lower <0 & outputs$CI.upper <0, "-", 
                       ifelse(outputs$CI.lower >0 & outputs$CI.upper >0, "+", "0") )
  
  #return it
  return(outputs)
}



#Plot the MCMC chains for each occupancy parameter
#@models is the output of getSpartaModels

plotModels<-function(models,param="psi.fs\\["){
  
  #create a directory to put the traceplots
  newdir <- file.path("model-outputs", "traceplots")
  dir.create(newdir,showWarnings = FALSE)
  if (!grepl(newdir,getwd())){
    setwd(newdir)
  }
  
  #remove any previous plotting files, if there are any
  if(length(list.files())>0){
    file.remove(list.files())
  }
  
  #for each species, do the following
  lapply(models,function(x){
    require(reshape2)
    #get the chains
    df <- melt(x$BUGSoutput$sims.array)
    #give sensible names
    names(df)<-c("Index","Chain","Parameter","Estimate")
    #subset to occupancy parameter
    df <- subset(df,grepl(param,df$Parameter))
    
    #if there are many parameters (i.e. years), just plot every other year
    if(length(unique(df$Parameter))>20){
      df$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", df$Parameter))
      df <- df[df$ParamNu%%2==0,]
    }
    
    #plot it
    require(ggplot2)
    ggplot(df)+
      geom_path(aes(x=Index,y=Estimate,colour=factor(Chain)))+
      facet_wrap(~Parameter,scales="free",ncol=4)+
      theme_bw()+
      theme(legend.position="none")
    ggsave(paste0(x$SPP_NAME,"_traceplot.png"))
    
  })
  
}

getOccurrenceMatrix<-function(df){
  require(reshape2)
  out<-acast(df,visit~Species,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
  out[out>0]<-1
  return(out)
}


getListLength<-function(df){
  require(plyr)
  out <- ddply(df,.(visit,Date,Year,month,day,MTB_Q),summarise,
        nuSpecies=length(unique(Species)),
        nuRecords=length(Species),
        Richness2=mean(Richness),
        RpS = length(Species)/Richness2,
        expertise = sum(Expertise),
        samplingSites = length(unique(interaction(lon,lat))))
  
  #add on some indices
  out$siteIndex <- as.numeric(factor(out$MTB_Q))
  out$yearIndex <- as.numeric(factor(out$Year))

  #sort dataset to match with the occurrence Matrix
  out <- arrange (out,visit)
}


getBUGSfits <- function(out,param="trend.year"){
  out <- as.data.frame(out)
  out <- out[grepl(param,out$Param),]
  out$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", out$Param))
  return(out)
}

combineModels <- function(out){
  
  names(out) <- mainSpecies
  
  out2<-ldply(out,function(x){
    output <- data.frame(x$summary)
    output$Param <- row.names(x$summary)
    return(output)
  })
  
  names(out2)[1]<-"Species"
  
  return(out2)
  
}


getFirstRange <- function(x){
  x <- as.character(x)
  strsplit(x,"-")[[1]][1]
}

removeGerman <- function(obs){
  obs <- as.character(obs)
  Encoding(obs) <- "latin1" 
  obs <- iconv(obs, "latin1", "ASCII", sub="")
}

formatObservers <- function(obs){
  require(gdata)
  obs <- as.character(obs)
  obs <- as.character(unlist(sapply(obs,function(x)strsplit(x,";"))))
  obs <- trim(as.character(unlist(sapply(obs,function(x)strsplit(x," u ")))))
  obs <- trim(as.character(unlist(sapply(obs,function(x)strsplit(x," und ")))))
  obs <- trim(as.character(unlist(sapply(obs,function(x)strsplit(x," u. ")))))
  obs <- trim(as.character(unlist(sapply(obs,function(x)strsplit(x," \\(EVS)")))))
  obs <- gsub("et al.","", obs)
  obs <- gsub("et al","", obs)
  obs <- sort(unique(obs))
  trim(obs)
}

fitBugs<-function(mySpecies,effort="nuSpecies",modelfile="R/BUGS_sparta.txt"){
  
  #organise data
  bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                    nyear = length(unique(listlengthDF$yearIndex)),
                    nvisit = nrow(listlengthDF),
                    site = listlengthDF$siteIndex,
                    year = listlengthDF$yearIndex,
                    yday = listlengthDF$day - median(listlengthDF$day),
                    nuSpecies = log(listlengthDF$nuSpecies) - log(median(listlengthDF$nuSpecies)),
                    nuRecs = log(listlengthDF$nuRecords) - log(median(listlengthDF$nuRecords)),
                    nuSS = listlengthDF$samplingSites - median(listlengthDF$samplingSites),# up to 3
                    expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                    RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                    y = as.numeric(occMatrix[,mySpecies]))
  listlengthDF$Species <- bugs.data$y
  
  #specify initial values
  library(reshape2)
  zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
  zst [is.infinite(zst)] <- 0
  inits <- function(){list(z = zst)}
  
  #fit model
  bugs.data$Effort <- bugs.data[[effort]]
  bugs.data$Effort_v2 <- bugs.data$RpS
  
  source('R/BUGS_misc_functions.R')
  
  #specify parameters to monitor
  params <- c("phenol.p","phenol2.p","effort.p","effort_2.p","expert.p","ss.p","psi.fs")
  
  #run model
  out1 <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
               n.chains=3, n.burnin=2000,n.iter=10000,parallel=T)
  
  return(out1)
}
