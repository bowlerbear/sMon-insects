########################################################################################

#need to include NAs for missing MTB data

Sys.time()

library(plyr)

#load in regional datasets
myfiles <- list.files("derived-data")
#myfiles <- list.files("/data/idiv_ess/Odonata")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

#combine these files
adultData <- ldply(adultFiles,function(x){
  out<-readRDS(paste("derived-data",x,sep="/"))
  #out<-readRDS(paste("/data/idiv_ess/Odonata",x,sep="/"))
  out$File <- x
  return(out)
})


#extract state from file name
adultData$State <- sapply(adultData$File,function(x)strsplit(x,"\\.rds")[[1]][1])
adultData$State <- sapply(adultData$State,function(x)strsplit(x,"_")[[1]][3])
nrow(adultData)#1023689
##########################################################################################

#add gbif data to fill gaps
gbifdata <- readRDS("derived-data/datafile_GBIF.rds")
#gbifdata <- readRDS("/data/idiv_ess/Odonata/datafile_GBIF.rds")
#nrow(gbifdata)#38191

#combine the two
adultData <- rbind(adultData,gbifdata)

######################################################################################

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

#####################################################################################

#remove MTBQs not in the shapefile = "51561" "50561" "49563"
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)
adultData <- subset(adultData, !MTB_Q %in% c("51561","50561","49563","55484","63012"))

###################################################################################
#filter to 1980 onwards

df <- subset(adultData, Year>=1980  & Year<2017)

######################################################################################

#pick species
#range expanding species
myspecies="Crocothemis erythraea"
stage="adult"

######################################################################################

#subset to phenology by regions

#read in and combine all phenology files
# phenolFiles<-list.files("/data/idiv_ess/Odonata")[grepl("speciesDays",list.files("/data/idiv_ess/Odonata"))]
# phenolData <- ldply(phenolFiles,function(x){
#   out<-read.delim(paste("/data/idiv_ess/Odonata",x,sep="/"))
#   out$File <- x
#   return(out)
# })
# 
# phenolFiles<-list.files()[grepl("speciesDays",list.files())]
# phenolData <- ldply(phenolFiles,function(x){
#  out<-read.delim(x)
#  out$File <- x
#  return(out)
# })
# 
# 
# #extract state from file name
# phenolData$State <- sapply(phenolData$File,function(x)strsplit(x,"\\.txt")[[1]][1])
# phenolData$State <- sapply(phenolData$State,function(x)strsplit(x,"_")[[1]][3])
# phenolData <- subset(phenolData, Species==myspecies)
# 
# df <- subset(df,interaction(yday,State) %in% interaction(phenolData$day,phenolData$State))

#if no phenolData for a given state, use max and min phenoldays
dfS <- subset(df, Species==myspecies)
obsPhenolData <- ddply(dfS,.(State),summarise,
                       minDay = round(quantile(yday,0.05)),
                       maxDay = round(quantile(yday,0.95)))

#expand to list all days between these days
obsPhenolData <- ddply(obsPhenolData,.(State),function(x){
  data.frame(Species=myspecies,
             day=as.numeric(x["minDay"]):as.numeric(x["maxDay"]),
             fits=NA,
             File=NA,
             State=x["State"])})

df <- subset(df,interaction(yday,State) %in% interaction(obsPhenolData$day,obsPhenolData$State))

#####################################################################################

#reduce size of the dataset

#reduce data to 5%%
#df <- df[sample(1:nrow(df),round(0.05*nrow(df))),]

#any oversampled plots???
out <- ddply(df,.(MTB_Q,Year),summarise,
             nuDates=length(unique(Date)))
#out <- arrange(out,desc(nuDates))
summary(out$nuDates)

#subset to at most 50 dates per year
nrow(df)
df <- ddply(df, .(Year,MTB_Q),function(x){
  mydates <- ifelse(length(unique(x$Date))>20,
                    sample(unique(x$Date),20),unique(x$Date))
  subset(x, Date %in% mydates)
})
nrow(df)

######################################################################################

df$visit <- paste(df$MTB_Q,df$Date,df$Beobachter,sep="_")

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
  out <- ddply(df,.(visit,Date,MTB_Q),summarise,
               nuSpecies=length(unique(Species)),
               nuRecords=length(Species),
               Richness2=mean(Richness),
               RpS = length(Species)/length(unique(Species)),
               expertise = sum(Expertise),
               samplingSites = length(unique(interaction(lat,lon))))
  
  #sort dataset to match with the occurrence Matrix
  out <- arrange (out,visit)
}
listlengthDF <- getListLength(df)

#get MTB
listlengthDF$MTB <- sapply(as.character(listlengthDF$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,(len-1))})

#rows of occuMatrix match visits
all(listlengthDF$visit==row.names(occMatrix))

##########################################################################################

#get nationwide boxes
load("mtbqsDF.RData")
#load("/data/idiv_ess/Odonata/mtbqsDF.RData")
all(df$MTB_Q%in% mtbqsDF$MTB_Q)
#unique(df$MTB_Q[!df$MTB_Q %in% mtbqsDF$MTB_Q])

#add on box info to listlength
listlengthDF$km50 <- mtbqsDF$km50[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
summary(listlengthDF$km50)
listlengthDF$km100 <- mtbqsDF$km100[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
summary(listlengthDF$km100)

#########################################################################################

#add on some indices
listlengthDF$State <- mtbqsDF$Counties[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- year(listlengthDF$Date)
listlengthDF$yday <- yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))
listlengthDF$stateIndex <- as.numeric(factor(listlengthDF$State))
listlengthDF$boxIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,listlengthDF$km50)))
listlengthDF$siteIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                   listlengthDF$MTB_Q)))
listlengthDF$mtbIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                  listlengthDF$MTB)))
#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("stateIndex","mtbIndex","siteIndex","boxIndex","MTB_Q")])
head(siteInfo)
save(siteInfo,file="siteInfo.RData")

#######################################################################################

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)
siteInfo <- arrange(siteInfo,siteIndex)
nrow(siteInfo)==max(siteInfo$siteIndex)

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
 
# #get info on whether species was seen at each site
speciesSite <- ddply(listlengthDF,.(siteIndex),summarise,PA=max(Species,na.rm=T))
siteInfo$obs <- speciesSite$PA[match(siteInfo$siteIndex,speciesSite$siteIndex)]

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
library(mgcv)
jags.ready <- jagam(obs ~ 1 + s(x,y,k=15), data=siteInfo, family=binomial,file="jagam.txt")

#extract bits for the model
bugs.data$X = jags.ready$jags.data$X[,-1]
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero[-1]
bugs.data$nspline = length(bugs.data$zero)
bugs.data$nspline1 = length(bugs.data$zero)+1
bugs.data$nspline2 = (length(bugs.data$zero))*2

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
#n.cores = 3

###########################################################################################

#modelfile="R/BUGS_dynamic_nation_spline.txt"
modelfile="/data/idiv_ess/Odonata/BUGS_dynamic_nation_spline.txt"

effort = "nuSpecies"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("z")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,n.adapt=1500,
            n.chains=n.cores, n.burnin=2000,n.iter=5000,parallel=T)
Sys.time()

#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_dynamicspline_",stage,"_", myspecies,".rds"))

########################################################################################