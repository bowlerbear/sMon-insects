########################################################################################

#load in regional datasets
#myfiles <- list.files("derived-data")
myfiles <- list.files("/data/idiv_ess/Odonata")

#read in and combine all adult files
adultFiles<-myfiles[grepl("adult_datafile",myfiles)]
adultFiles<-adultFiles[grepl("rds",adultFiles)]

#combine these files
library(plyr)
adultData <- ldply(adultFiles,function(x){
  #out<-readRDS(paste("derived-data",x,sep="/"))
  out<-readRDS(paste("/data/idiv_ess/Odonata",x,sep="/"))
  out$File <- x
  return(out)
})

#extract state from file name
adultData$State <- sapply(adultData$File,function(x)strsplit(x,"\\.rds")[[1]][1])
adultData$State <- sapply(adultData$State,function(x)strsplit(x,"_")[[1]][3])

#format date
library(lubridate)
adultData$Date<-as.Date(adultData$Date)
adultData$Year <- year(adultData$Date)
adultData$yday <- yday(adultData$Date)
adultData$week <- week(adultData$Date)
adultData$Month <- month(adultData$Date)
adultData$Day <- day(adultData$Date)
adultData <- subset(adultData,!is.na(Date))
adultData <- subset(adultData,yday!=1)
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)

###################################################################################
#filter to 1980 onwards

df <- subset(adultData, Year>=1980  & Year<2017)
origDF <- df

######################################################################################

#pick species
myspecies="Cordulegaster boltonii"
stage="adult"

######################################################################################

#subset to phenology by regions

#read in and combine all phenology files
 phenolFiles<-list.files("/data/idiv_ess/Odonata")[grepl("speciesDays",list.files("/data/idiv_ess/Odonata"))]
 phenolData <- ldply(phenolFiles,function(x){
   out<-read.delim(paste("/data/idiv_ess/Odonata",x,sep="/"))
   out$File <- x
   return(out)
 })

# phenolFiles<-list.files()[grepl("speciesDays",list.files())]
# phenolData <- ldply(phenolFiles,function(x){
#   out<-read.delim(x)
#   out$File <- x
#   return(out)
# })


#extract state from file name
phenolData$State <- sapply(phenolData$File,function(x)strsplit(x,"\\.txt")[[1]][1])
phenolData$State <- sapply(phenolData$State,function(x)strsplit(x,"_")[[1]][3])
phenolData <- subset(phenolData, Species==myspecies)

df <- subset(df,interaction(yday,State) %in% interaction(phenolData$day,phenolData$State))

######################################################################################

#define a visit
df$visit <- paste(df$MTB_Q,df$Date,df$Beobachter,sep="_")

#get rid of problem cells temporarily
df <- subset(df, !MTB_Q %in% c("50561","51561","55484","63012"))

#get occurence matrix  - detection and non-detection
getOccurrenceMatrix<-function(df){
  require(reshape2)
  out<-acast(df,visit~Species,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
  out[out>0]<-1
  return(out)
}
occMatrix <- getOccurrenceMatrix(df)

#get list length
getListLength<-function(df){
  require(plyr)
  out <- ddply(df,.(visit,Date,Year,Month,Day,yday,MTB_Q,State),summarise,
               nuSpecies=length(unique(Species)),
               nuRecords=length(Species),
               Richness2=mean(Richness),
               RpS = length(Species)/length(unique(Species)),
               expertise = sum(Expertise),
               samplingSites = length(unique(interaction(lat,lon))))
  
  #add on some indices
  out$yearIndex <- as.numeric(factor(out$Year))
  out$stateIndex <- as.numeric(factor(out$State))

  #sort dataset to match with the occurrence Matrix
  out <- arrange (out,visit)
}
listlengthDF <- getListLength(df)

#rows of occuMatrix match visits
#all(listlengthDF$visit==row.names(occMatrix))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

##########################################################################################

#get nationwide boxes
load("/data/idiv_ess/Odonata/mtbqsDF.RData")
#load("mtbqsDF.RData")
all(df$MTB_Q%in% mtbqsDF$MTB_Q)

#add on box info to listlength
listlengthDF$km50 <- mtbqsDF$km50[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
summary(listlengthDF$km50)

#add in dummy for small number of NAs - need to fix this..
listlengthDF$km50[is.na(listlengthDF$km50)]<-length(unique(listlengthDF$km50))+1
listlengthDF$boxIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,listlengthDF$km50)))
listlengthDF$siteIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                   listlengthDF$km50,
                                                   listlengthDF$MTB_Q)))
listlengthDF <- arrange(listlengthDF,visit)

#########################################################################################

#Examine amount of data per 50 km box
#df <- subset(df, Species==myspecies)
# origDF <- merge(origDF,mtbqsDF,by="MTB_Q",all=T)
# boxData <- ddply(origDF, .(km50),summarise,
#                  nuYears = length(unique(Year)),
#                  nuRecs = length(Species),
#                  nuSpecies = length(unique(Species)))
# 
# #plotting
# library(raster)
# library(sp)
# germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
# myGrid <- raster('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/km50grid.tif')
# germanAdmin <- spTransform(germanAdmin,projection(myGrid))
# 
# myGrid[] <- NA
# myGrid[boxData$km50] <- boxData$nuYears
# plot(myGrid)
# plot(germanAdmin,add=T)
# 
# myGrid[] <- NA
# myGrid[boxData$km50] <- log(boxData$nuRecs)
# plot(myGrid)
# plot(germanAdmin,add=T)
# 
# myGrid[] <- NA
# myGrid[boxData$km50] <- boxData$nuSpecies
# plot(myGrid)
# plot(germanAdmin,add=T)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("MTB_Q","stateIndex","siteIndex","boxIndex")])
siteInfo <- arrange(siteInfo,stateIndex,boxIndex,siteIndex)
head(siteInfo)

nrow(siteInfo)==length(unique(listlengthDF$siteIndex))

########################################################################################

#fit nation-wide model with random slopes to each box

########################################################################################

#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nstate = length(unique(siteInfo$stateIndex)),
                  nbox = length(unique(siteInfo$boxIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  state = listlengthDF$stateIndex,
                  box = listlengthDF$boxIndex,
                  year = listlengthDF$yearIndex,
                  stateS = siteInfo$stateIndex,
                  boxS = siteInfo$boxIndex,
                  yday = listlengthDF$yday - median(listlengthDF$yday),
                  nuSpecies = log(listlengthDF$nuSpecies) - log(median(listlengthDF$nuSpecies)),
                  singleList = listlengthDF$singleList,
                  shortList = listlengthDF$shortList,
                  nuRecs = log(listlengthDF$nuRecords) - log(median(listlengthDF$nuRecords)),
                  nuSS = log(listlengthDF$samplingSites) - log(median(listlengthDF$samplingSites)),# up to 3
                  expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                  RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                  y = as.numeric(occMatrix[,myspecies]))
listlengthDF$Species <- bugs.data$y

all(row.names(occMatrix)==listlengthDF$visit)

#######################################################################################

# #set up the spline info
siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
#  
# #get info on whether species was seen at each site
 speciesSite <- ddply(listlengthDF,.(siteIndex),summarise,PA=max(Species,na.rm=T))
 siteInfo$obs <- speciesSite$PA[match(siteInfo$siteIndex,speciesSite$siteIndex)]
 siteInfo <- subset(siteInfo,!is.na(obs))
 subset(siteInfo,is.na(x))#63012,50561,51561
 siteInfo <- subset(siteInfo,!is.na(x))
# 
# #fit normal gam
# library(mgcv)
# gam1 <- gam(obs ~ 1 + s(x,y), data=siteInfo, family=binomial)
# siteInfo$fits <- gam1$fitted.values
# library(ggplot2)
# qplot(x,y,data=siteInfo,colour=fits)+
#   scale_colour_gradient(low="blue",high="red")
# 
# ##########################################################################################
# 
# #set up space-time-spline info
# siteInfo <- unique(listlengthDF[,c("MTB_Q","stateIndex","siteIndex","boxIndex","yearIndex")])
# siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
# siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
# speciesSite <- ddply(listlengthDF,.(siteIndex,yearIndex),summarise,PA=max(Species,na.rm=T))
# siteInfo$obs <- speciesSite$PA[match(interaction(siteInfo$siteIndex,siteInfo$yearIndex),
#                                      interaction(speciesSite$siteIndex,speciesSite$yearIndex))]
# siteInfo <- subset(siteInfo,!is.na(obs))
# siteInfo <- subset(siteInfo,!is.na(x))
# 
# #fit normal gam
# library(mgcv)
# gam1 <- gam(obs ~ s(yearIndex) + s(x,y,yearIndex), data=siteInfo, family=binomial)
# siteInfo$fits <- gam1$fitted.values
# library(ggplot2)
# qplot(x,y,data=siteInfo,colour=fits)+
#   scale_colour_gradient(low="blue",high="red")+
#   facet_wrap(~yearIndex)

#######################################################################################

#get BUGS code using jagam
#difference in nrow(siteInfo) and bugs.data$nsite
library(mgcv)
jags.ready <- jagam(obs ~ 1 + s(x,y), data=siteInfo, family=binomial,file="jagam.txt")

#extract bits for the model
 
#including intercept
bugs.data$X = jags.ready$jags.data$X
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero
bugs.data$nspline = length(bugs.data$zero)-1
bugs.data$nspline1 = length(bugs.data$zero)
bugs.data$nspline2 = (length(bugs.data$zero)-1)*2

########################################################################################

#specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

########################################################################################

#define model params
ni <- 6000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3

#fit model
library(rjags)
library(R2WinBUGS)
library(jagsUI)

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

###########################################################################################

#modelfile="R/BUGS_dynamic_spline.txt"
modelfile="/data/idiv_ess/Odonata/BUGS_dynamic_spline.txt"

effort = "nuSpecies"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("phenol.p","phenol2.p","effort.p","psi.fs","single.p","eta",
            "a.persist","a.colonize")

#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
            n.chains=n.cores, n.burnin=3000,n.iter=10000,parallel=T)

#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_dynamicspline_",stage,"_", myspecies,".rds"))

########################################################################################