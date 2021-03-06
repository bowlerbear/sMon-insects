#given the directory in which the sparta model outputs are found, put all the models in a list
#@param dir = directory of sparta model outputs

# getSpartaModels<-function(dir){
#   
#   if(!dir%in%list.files()) stop('no such directory found')  
#   
#   sp_mods <- list.files(dir)[grepl(".rdata",list.files(dir))] 
#   
#   models <- lapply(sp_mods, function(sp) {
#     load(file=paste0(paste0(dir,"/"),sp))
#     return(out)
#   }
#   )
#   
#   names(models) <- gsub(sp_mods, pa="\\.rdata", repl="")
#   
#   return(models)
# }

library(jagsUI)

getSpartaModels<-function(dir){
  
  sp_mods <- list.files(dir)[grepl(".rdata",list.files(dir))]
  require(plyr)
  ldply(sp_mods,function(file){
  load(file)#object is called out
  modelSummary <- data.frame(out$BUGSoutput$summary)
  modelSummary$Param <- row.names(modelSummary)
  modelSummary$File <- file
  return(modelSummary)
  }
  )
  
}

#for other models
getModelSummaries <- function(dir){
  myFiles <- list.files(dir)
  require(plyr)
  temp <- ldply(myFiles,function(x){
    out <- readRDS(paste(dir,x,sep="/"))
    outS <- data.frame(out$summary)
    outS$Param <- row.names(out$summary)
    outS$File <- x
    return(outS)
  })
  return(temp)
}



getModels <- function(dir){
  myFiles <- list.files(dir)
  require(plyr)
  temp <- ldply(myFiles,function(x){
    out <- readRDS(paste(dir,x,sep="/"))
    outS <- data.frame(out)
    outS$Param <- row.names(out)
    outS$File <- x
    return(outS)
  })
  return(temp)
}


plotTS <- function(x){
  require(ggplot2)
  g1 <- ggplot(x)+
    geom_line(aes(x=Year,y=mean))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
    facet_wrap(~Species)+
    theme_bw()+
    theme(axis.text = element_text(size =rel(0.55))) +
    ylab("Predicted occupancy")
  print(g1)
}

plotDetections <- function(x){
  require(ggplot2)
  g1 <- ggplot(x)+
    geom_line(aes(x=Year,y=mean))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
    facet_wrap(~Code,ncol=6)+
    theme_bw()+
    theme(axis.text = element_text(size =rel(0.55))) +
    ylab("Detection probability")
  print(g1)
}

plotTS_scaled <- function(x){
  require(ggplot2)
  g1 <- ggplot(x)+
    geom_line(aes(x=Year,y=mean))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
    theme_bw()+
    facet_wrap(~Code,ncol=6,scale="free")+
    theme(axis.text = element_text(size =rel(0.55))) +
    ylab("Predicted occupancy")
  print(g1)
}

plotTSwithRugs <- function(x,y){
  require(ggplot2)
  g1 <- ggplot(x)+
    geom_line(aes(x=Year,y=mean))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
    theme(axis.text = element_text(size =rel(0.55))) +
    ylab("Predicted occupancy")+
    geom_rug(data=y,aes(x=Year,y=dummy,colour=Obs),sides="b",size=2)+
    scale_colour_viridis_c(name = "Number of observations",trans="log10")+
    theme_bw()+
    facet_wrap(~Code,ncol=4,scale="free")+
    theme(legend.position = "top")
  print(g1)
}


plotTS_regional <- function(x){
  require(ggplot2)
  if(!"Naturraum" %in%  names(x)){
    g1 <- ggplot(data=x,aes(x=Year,y=mean,group=factor(Site)))+
      geom_line(aes(x=Year,y=mean,colour=factor(Site)))+
      geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.,fill=factor(Site)),alpha=0.5)+
      facet_wrap(~Code)+ylab("Predicted occupancy")+
      theme(legend.position="none")
  }else{
    g1 <- ggplot(data=x,aes(x=Year,y=mean,group=Naturraum))+
      geom_line(aes(x=Year,y=mean,colour=Naturraum))+
      geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.,fill=Naturraum),alpha=0.5)+
      facet_wrap(~Code)+ylab("Predicted occupancy")+
      theme(legend.position="none")
  }
  
  print(g1)
}


tidyTS <- function(annualDFS){
  
  ggplot(annualDFS)+
    geom_line(aes(x=Year,y=mean))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
    xlab("Year")+
    ylab("Occupancy proportion")+
    theme_bw()
  
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



tracePlot <- function(x){
  library(coda)
  
  #get the chains
  plot(out$samples,trace=TRUE,density=FALSE,ask=TRUE)
  
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


getBUGSfitsS <- function(out,param="psi.fs"){
  out <- as.data.frame(out$summary)
  if(all(!grepl("Param",names(out)))){out$Param <- row.names(out)}
  out <- out[grepl(param,out$Param),]
  out$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", out$Param))
  return(out)
}

getBUGSfits <- function(out,param="trend.year"){
  out <- as.data.frame(out)
  if(all(!grepl("Param",names(out)))){out$Param <- row.names(out)}
  out <- out[grepl(param,out$Param),]
  out$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", out$Param))
  return(out)
}


getBUGSfitsII <- function(out,param="psi.state"){
  out <- as.data.frame(out)
  if(all(!grepl("Param",names(out)))){out$Param <- row.names(out)}
  out <- out[grepl(param,out$Param),]
  out$ParamNu <- sub(".*\\[([^][]+)].*", "\\1", out$Param)
  out$Year <- as.numeric(sapply(out$ParamNu,function(x)strsplit(x,",")[[1]][2]))
  out$Site <- as.numeric(sapply(out$ParamNu,function(x)strsplit(x,",")[[1]][1]))
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
  obs <- trim(as.character(unlist(sapply(obs,function(x)strsplit(x,"& ")))))
  obs <- trim(as.character(unlist(sapply(obs,function(x)strsplit(x," \\(EVS)")))))
  obs <- gsub("et al.","", obs)
  obs <- gsub("et al","", obs)
  obs <- sort(unique(obs))
  trim(obs)
}

getCodeFromFile <- function(modelDF,myfile="out_sparta_nation_naturraum_adult_"){
  modelDF$Species <- gsub(myfile,"",modelDF$File)
  modelDF$Species <- gsub(".rds","",modelDF$Species)
  modelDF$Genus <- sapply(modelDF$Species, function(x)strsplit(x," ")[[1]][1])
  modelDF$Genus <- sapply(modelDF$Genus, function(x) substr(x,1,3))
  modelDF$Spec <- sapply(modelDF$Species, function(x)strsplit(x," ")[[1]][2])
  modelDF$Spec <- sapply(modelDF$Spec, function(x) substr(x,1,3))
  modelDF$Code <- paste(modelDF$Genus,modelDF$Spec,sep="_")
  return(modelDF)
}

getSurnames <- function(x){
  surName <- strsplit(x," ")
  surNames <- as.character(sapply(surName[[1]],function(x)trim(x)))
  nu <- length(surNames)
  return(surNames[nu])
}

checkDuplicates <- function(observerSurnames){
  dups <- observerSurnames[duplicated(observerSurnames)]
  sapply(dups,function(x)observerSurnames[grepl(x,observerSurnames)])
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
               n.chains=3, n.burnin=25000,n.iter=100000,parallel=T)
  
  return(out1)
}



fitTrends <- function(df){
  

  #combine for a bugs model
  bugs.data<-list(
    Index = df$mean,
    SE = df$sd,
    n.year = nrow(df))
  
  
  #write a bugs model to estimate trends
  cat("
      model{
      
      #give indices with error
      for(t in 1:n.year){
      tau[t] <- pow(SE[t],-2)
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
      ",fill=TRUE,file="C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/Bugs_trend.txt")
  
  #run model
  source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/BUGS_misc_functions.R')
  
  #specify parameters to monitor
  params <- c("int","trend")
  
  inits <- function(){
    list(int = rnorm(1,0,0.1), trend= rnorm(1,0,0.1))}
  
  #run model
  out1 <- jags(bugs.data, inits=inits, params, "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/Bugs_trend.txt", n.thin=nt,
               n.chains=3, n.burnin=1000,n.iter=5000,parallel=T)
  
  #return predictions
  return(data.frame(trend=out1$mean$trend,trend_sd=out1$sd$trend,lowerCI=out1$q2.5$trend,
         upperCI=out1$q97.5$trend))
  
}


fitTrendsBeta <- function(df){
  
  
  #combine for a bugs model
  bugs.data<-list(
    Index = df$mean,
    SE = df$sd,
    n.year = nrow(df))
  
  
  #write a bugs model to estimate trends
  cat("
      model{
      
      #give indices with error
      for(t in 1:n.year){
      tau[t] <- pow(SE[t],-2)
      Index[t] <- ilogit(iIndex[t])
      iIndex[t] ~ dnorm(actPop[t],logit(tau[t]))
      }
      
      #fit traits model
      int ~ dnorm(0,0.1)
      trend ~ dnorm(0,0.1)
      phi ~ dgamma(.1,.1)

      #fit model
      for(t in 1:n.year){
      actPop[t] ~ dbeta(alpha[t], beta[t]) T(0.0001,0.9999)
      alpha[t] <- mu[t] * phi
      beta[t]  <- (1-mu[t]) * phi
      #linear predictor
      logit(mu[t]) <- int + trend * t 
      }
      
      
      }
      ",fill=TRUE,file="C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/Bugs_trend.txt")
  
  #run model
  source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/BUGS_misc_functions.R')
  
  #specify parameters to monitor
  params <- c("int","trend")
  
  inits <- function(){
    list(int = rnorm(1,0,0.1), trend= rnorm(1,0,0.1))}
  
  #run model
  out1 <- jags(bugs.data, inits=inits, params, "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/Bugs_trend.txt", n.thin=nt,
               n.chains=3, n.burnin=1000,n.iter=5000,parallel=T)
  
  #return predictions
  return(data.frame(trend=out1$mean$trend,trend_sd=out1$sd$trend,lowerCI=out1$q2.5$trend,
         upperCI=out1$q97.5$trend))
  
}


addBest <- function(lm1,myresponse="trend"){
  
  #get missing covariates
  tested <- names(lm1$model)[-1]
  notTested <- mytraits[!mytraits %in% tested]
  
  new <- ldply(notTested,function(x){
    myformula <- paste(as.character(formula(lm1))[3],x,sep="+")
    myformula <- paste0(myresponse,"~",myformula)
    lm2 <- update(lm1,formula=as.formula(myformula))
    getGLMSummary(lm2)[(length(tested)+2),]
  })
  
  old <- getGLMSummary(lm1)
  
  #return full
  rbind(old,new)
  
}

getGLMSummary <- function(glm1){
  coefDF <- data.frame(summary(glm1)$coef)
  coefDF$lowerQ <- confint(glm1)[,1]
  coefDF$upperQ <- confint(glm1)[,2]
  coefDF$param <- row.names(coefDF)
  return(coefDF)
}

plotCluster <- function(x){
  require(ggplot2)
  g1 <- ggplot(x)+
    geom_line(aes(x=Year,y=mean,colour=Species))+
    facet_wrap(~cluster)+
    theme_bw()+
    ylab("predicted occupancy")+
    theme(legend.position="none")
  print(g1)
}

deriv <- function(x, y) diff(y) / diff(x)
middle_pts <- function(x) x[-1] - diff(x) / 2

getBestCluster <- function(temp,x){
  temp2 <- subset(temp,Param==x)
  firstderiv <- deriv(temp2$variable, temp2$value)
  #plot(firstderiv ~ temp2$variable[-1])
  second_d <- deriv(middle_pts(temp2$variable), firstderiv)
  #plot(second_d ~ temp2$variable[-c(1:2)])
  df <- data.frame(second_d = second_d, 
                   midpts = middle_pts(middle_pts(temp2$variable)),
                   Param=x)
  subset(df,second_d==max(second_d))
}

getClusterStats <- function(mydiss,myclustering){
  require(fpc)
  clustDF <- cluster.stats(mydiss, 
                           clustering=myclustering)
  data.frame(minCluster=clustDF$min.cluster.size,
             avBetween=clustDF$average.between,
             avWithin=clustDF$average.within,
             silWidth=clustDF$avg.silwidth,
             dunn=clustDF$dunn,
             sepIndex=clustDF$sindex)
}

plotTSclust <- function(mylist){
  require(ggplot2)
  require(reshape2)
  mylistMelt <- melt(mylist,id="Cluster")
  qplot(x=Cluster,y=value,data=mylistMelt)+
    facet_wrap(~variable,scales="free")
}


#from defunct SDMTools package

read.asc <-
  function (file, gz=FALSE) {
    #confirm ascii grid file name is specified
    if (gz) {
      if (substr(file, nchar(file) - 2, nchar(file)) != ".gz") stop("not a valid .gz file")
    } else {
      if (substr(tolower(file), nchar(file) - 3, nchar(file)) != ".asc") stop("not a valid .asc file")
    }
    
    #read in the header
    if (gz) { zz <- gzfile(file, "r") } else { zz <- file(file, "r") }
    nc <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); nc <- as.numeric(nc[[2]][1])#number of columns
    nl <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); nl <- as.numeric(nl[[2]][1])#number of rows
    xll <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); #lower left corner
    yll <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); #lower left corner
    cs <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); cs <- as.numeric(cs[[2]][1])#cell size
    nas <- scan(zz,what=list('',''),nlines=1,quiet=TRUE); nas <- as.numeric(nas[[2]][1])#nodata value
    #close the link to the file
    close(zz)
    
    #ensure xll & yll are centers of the cells
    if ((xll[[1]][1] == "xllcenter") | (xll[[1]][1] == "XLLCENTER")) { xll=as.numeric(xll[[2]][1]) } else { xll=as.numeric(xll[[2]][1])+ cs/2 }
    if ((yll[[1]][1] == "yllcenter") | (xll[[1]][1] == "YLLCENTER")) { yll=as.numeric(yll[[2]][1]) } else { yll=as.numeric(yll[[2]][1])+ cs/2 }
    
    #read in the data skipping the first six header rows
    if (gz) { zz <- gzfile(file, "r"); output <- scan(zz,nmax=nl*nc,skip=6,quiet = TRUE); close(zz);
    } else { output <- scan(file,nmax=nl*nc,skip=6, quiet = TRUE) }
    
    #convert no data to NA
    output[output == nas] <- NA
    #convert data to matrix
    output <- matrix(c(as.matrix(output)), ncol = nl, nrow = nc)
    output <- output[, ncol(output):1]
    #define the attributes
    attr(output, "xll") <- xll
    attr(output, "yll") <- yll
    attr(output, "cellsize") <- cs
    attr(output, "type") <- 'numeric'
    class(output) <- "asc"
    #return the file
    return(output)
  }

floor_decade    = function(value){ return(value - value %% 10) }

#taken from https://stat.ethz.ch/pipermail/r-sig-phylo/2013-August/003017.html
as.phylo.formula<-function (x, data = parent.frame(), ...){
  require(phytools)
  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
  if (length(x) != 2) 
    stop(err)
  if (x[[1]] != "~") 
    stop(err)
  f <- x[[2]]
  taxo <- list()
  while (length(f) == 3) {
    if (f[[1]] != "/") 
      stop(err)
    if (!is.factor(data[[deparse(f[[3]])]])) 
      stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if (length(f) > 1) 
      f <- f[[2]]
  }
  if (!is.factor(data[[deparse(f)]])) 
    stop(paste("Variable", deparse(f), "must be a factor."))
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[, 1])
  taxo.data[, 1] <- 1:nrow(taxo.data)
  f.rec <- function(subtaxo) {
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[, u])
    if (u == 1) {
      if (length(levels) != nrow(subtaxo)) 
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[, 1]))
    }
    t <- character(length(levels))
    for (l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)])
      t[l] <- paste("(", paste(x, collapse = ","), ")", sep = "")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse = ","),");", sep = "")
  phy <- read.newick(text = string) ## so that singles will be read without error
  phy$edge.length <- rep(1,nrow(phy$edge)) #here set eddge length to 1
  phy <- collapse.singles(phy)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


addMTBQ <- function(mtbqs){
  
  mtbqs$MTB_Q= ifelse(mtbqs$Quadrant == "NW",paste(mtbqs$Value, 1, sep = ""),
                      ifelse(mtbqs$Quadrant == "NO", paste(mtbqs$Value, 2, sep = ""),
                             ifelse(mtbqs$Quadrant == "SW", paste(mtbqs$Value, 3, sep = ""),
                                    paste(mtbqs$Value, 4, sep = ""))))
  return(mtbqs)
}

#general functions
logit <- function(x){
  l <- log(x/(1-x))
  return(l)
}

expit <- function(x){
  e <- exp(x)/(1+exp(x))
  return(e)
}

addTrendClassification <- function(x){
  x$Trend <- "stable"
  x$Trend[x$X2.5.>0 & x$X97.5.>0] <- "increasing"
  x$Trend[x$X2.5.<0 & x$X97.5.<0] <- "decreasing"
  return(x)
}

