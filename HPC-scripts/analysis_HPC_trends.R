
###natural-raum level############################

#read in model output
tdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum/5914536"
modelFiles <- list.files(tdir)
#modelFiles are the list of species model outputs


library(plyr)
modelTrends <- ldply(modelFiles,function(x){
  
  out <- readRDS(paste(tdir,x,sep="/"))
  outS <- data.frame(out$summary)
  outS$Param <- row.names(out$summary)
  
  
  #get bugs data for trend model
  temp <- ldply(1:7,function(i){
  bugs.data <- list(meanOcc = outS$mean[grepl(paste0("psi.raum\\[",i),outS$Param)],
                    sdOcc = outS$sd[grepl(paste0("psi.raum\\[",i),outS$Param)],
                    nyear = length(outS$mean[grepl(paste0("psi.raum\\[",i),outS$Param)]))
  
  #the below are used the linear regression model in the model file -see below
  bugs.data$sumX <- sum(1:bugs.data$nyear)
  bugs.data$sumX2 <- sum((1:bugs.data$nyear)^2)
  bugs.data$dummy <- rpois(bugs.data$nyear,3)

  #fit model
  library(rjags)
  library(R2WinBUGS)
  library(jagsUI)
  
  #define model params
  modelfile="R/BUGS_national_occupancy_trends.txt"
  ni <- 6000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3
  
  params <- c("regres.psi")#you may want to also monitor other things, e.g., totalChange,CoV
  
  #run model
  out <- jags(bugs.data, inits=NULL, params, modelfile, n.thin=nt,
              n.chains=nc, n.burnin=nb,n.iter=ni,parallel=T)
  
  data.frame(t(out$summary[1,]),file=x,naturraum=i)
  
})
  return(temp)
  
  
})


saveRDS(modelTrends,file="modelTrends_naturraum_trends.rds")


###national analysis################################

#naturraum model
tdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum/5914536"
modelFiles <- list.files(tdir)

#state model
tdir<-"C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/Odonata_adult_nation_state"
modelFiles <- list.files(tdir)

library(plyr)
modelTrends <- ldply(modelFiles,function(x){
  
  out <- readRDS(paste(tdir,x,sep="/"))
  outS <- data.frame(out$summary)
  outS$Param <- row.names(out$summary)
  
  #get bugs data for trend model
  bugs.data <- list(meanOcc = outS$mean[grepl("psi.fs",outS$Param)],
                    sdOcc = outS$sd[grepl("psi.fs",outS$Param)],
                    nyear = length(outS$mean[grepl("psi.fs",outS$Param)]))
  
  bugs.data$sumX <- sum(1:bugs.data$nyear)
  bugs.data$sumX2 <- sum((1:bugs.data$nyear)^2)
  bugs.data$dummy <- rpois(bugs.data$nyear,3)
  
  #fit model
  library(rjags)
  library(R2WinBUGS)
  library(jagsUI)
  
  #define model params
  modelfile="R/BUGS_national_occupancy_trends.txt"
  ni <- 6000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3
  
  params <- c("regres.psi")
  
  #run model
  out <- jags(bugs.data, inits=NULL, params, modelfile, n.thin=nt,
              n.chains=nc, n.burnin=nb,n.iter=ni,parallel=T)
  
  data.frame(t(out$summary[1,]),file=x)
  
})

saveRDS(modelTrends,file="model-outputs/modelTrends_trends.rds")#with naturraum models


###sparta trend#####################################

tdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6288453"
modelFiles <- list.files(tdir)

library(plyr)
modelTrends <- ldply(modelFiles,function(x){
  
  out <- readRDS(paste(tdir,x,sep="/"))
  outS <- data.frame(out$summary)
  outS$Param <- row.names(out$summary)
  
  #get bugs data for trend model
  bugs.data <- list(meanOcc = outS$mean[grepl("psi.fs",outS$Param)],
                    sdOcc = outS$sd[grepl("psi.fs",outS$Param)],
                    nyear = length(outS$mean[grepl("psi.fs",outS$Param)]))
  
  bugs.data$sumX <- sum(1:bugs.data$nyear)
  bugs.data$sumX2 <- sum((1:bugs.data$nyear)^2)
  bugs.data$dummy <- rpois(bugs.data$nyear,3)
  
  #fit model
  library(rjags)
  library(R2WinBUGS)
  library(jagsUI)
  
  #define model params
  modelfile="R/BUGS_national_occupancy_trends.txt"
  ni <- 6000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3
  
  params <- c("regres.psi","totalChange","CoV",
              "totalChangeGr","totalChangeLogit")
  
  #run model
  out <- jags(bugs.data, inits=NULL, params, modelfile, n.thin=nt,
              n.chains=nc, n.burnin=nb,n.iter=ni,parallel=T)
  
  temp <- data.frame(out$summary)
  temp$Param <- row.names(out$summary)
  temp$file <- x
  return(temp)
  
})

saveRDS(modelTrends,file="model-outputs/modelTrends_sparta_trends.rds")

