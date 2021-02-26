#get all libraries we need
#suppressMessages()
suppressMessages(library(rjags))
suppressMessages(library(R2WinBUGS))
suppressMessages(library(jagsUI))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))


#load the relational table of task ids and species
#speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_adult.txt"),as.is=T)
##get task id
#task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 
#get species for this task
#myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 

#set stage
stage="adult"

#set seed
set.seed(3)

#number of MCMC samples
#niterations = 1000

Sys.time()

#load in regional datasets
myfiles <- list.files("derived-data")
#myfiles <- list.files("/data/idiv_ess/Odonata")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

#exclude MTB64 files
adultFiles <- adultFiles[!grepl("MTB64",adultFiles)]

#exclude original files if an updated file is available
updatedFiles <- adultFiles[grepl("updated",adultFiles)]
updatedFiles <- gsub("_updated","",updatedFiles)
adultFiles <- adultFiles[!adultFiles %in% updatedFiles]


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
gbifdata <- readRDS("derived-data/datafile_iNaturalist.rds")
#gbifdata <- readRDS("/data/idiv_ess/Odonata/datafile_GBIF.rds")
#nrow(gbifdata)#38191

#combine the two
adultData <- rbind(adultData,gbifdata)

##########################################################################################

#format date
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

df <- subset(adultData, Year>=1990  & Year<=2016)

### mtb #########################################################

df$MTB <- sapply(as.character(df$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,len-1)})

df <- subset(df,!MTB %in% c("6301","5548","5156","5056","4956","4455"))

#####################################################################################

#remove sites visited once
# siteSummary <- ddply(df,.(MTB_Q),summarise,nuYears=length(unique(Year)))
# df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])

######################################################################################

#get nationwide boxes
load("mtbqsDF.RData")
names(mtbqsDF)[which(names(mtbqsDF)=="Value")] <- "MTB"
#load("/data/idiv_ess/Odonata/mtbqsDF.RData")

#this should be TRUE...
all(df$MTB%in% mtbqsDF$MTB)
#unique(df$MTB[!df$MTB %in% mtbqsDF$MTB])

df$State <- mtbqsDF$Counties[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$State))

#check missing naturruam data
mtbqsDF$Natur[is.na(mtbqsDF$Natur)] <- mtbqsDF$MTB_Natur[is.na(mtbqsDF$Natur)]
df$Natur <- mtbqsDF$Natur[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$Natur))
unique(df$MTB[is.na(df$Natur)])

mtbqsDF$CoarseNatur[is.na(mtbqsDF$CoarseNatur)] <- mtbqsDF$MTB_CoarseNatur[is.na(mtbqsDF$CoarseNatur)]
df$CoarseNatur <- mtbqsDF$CoarseNatur[match(df$MTB,mtbqsDF$MTB)]
unique(df$MTB[is.na(df$CoarseNatur)])

#####################################################################################

#subset by average phenology across whole germany

# dfS <- subset(df, Species==myspecies)
# 
# #across all species
# obsPhenolData <- summarise(df,
#                            minDay = round(quantile(yday,0.05)),
#                            maxDay = round(quantile(yday,0.95)))
# df <- subset(df, yday > obsPhenolData$minDay & yday < obsPhenolData$maxDay)

######################################################################################

#define a visit
df$visit <- paste(df$MTB,df$Date,sep="_")
#df$visit <- paste(df$MTB,df$Date,df$Beobachter,sep="_")

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
  out <- ddply(df,.(visit,Date,MTB),summarise,
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

#rows of occuMatrix match visits
all(listlengthDF$visit==row.names(occMatrix))

#######################################################################################

#add on some indices
listlengthDF$State <- mtbqsDF$Counties[match(listlengthDF$MTB,mtbqsDF$MTB)]
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- year(listlengthDF$Date)
listlengthDF$yday <- yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))
listlengthDF$stateIndex <- as.numeric(factor(listlengthDF$State))
listlengthDF$boxIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,listlengthDF$km50)))

#mtb as site
listlengthDF$siteIndex <- as.numeric(factor(listlengthDF$MTB))

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB,mtbqsDF$MTB)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB,mtbqsDF$MTB)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
subset(listlengthDF,is.na(Naturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("siteIndex","MTB","nnIndex","cnIndex")])
head(siteInfo)

######################################################################################

#pick species

#practise species
#myspecies="Aeshna cyanea"
#myspecies = "Crocothemis erythraea"
myspecies = "Sympetrum danae"
#stage="adult"

########################################################################################

#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nraum = length(unique(siteInfo$nnIndex)),
                  ncraum = length(unique(siteInfo$cnIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  raum = listlengthDF$nnIndex,
                  craum = listlengthDF$cnIndex,
                  year = listlengthDF$yearIndex,
                  craumS = siteInfo$cnIndex,
                  raumS = siteInfo$nnIndex,
                  yday = listlengthDF$yday - median(listlengthDF$yday),
                  yday2 = listlengthDF$yday^2 - median(listlengthDF$yday^2),
                  nuSpecies = log(listlengthDF$nuSpecies) - log(median(listlengthDF$nuSpecies)),
                  singleList = listlengthDF$singleList,
                  shortList = listlengthDF$shortList,
                  nuRecs = log(listlengthDF$nuRecords) - log(median(listlengthDF$nuRecords)),
                  nuSS = log(listlengthDF$samplingSites) - log(median(listlengthDF$samplingSites)),# up to 3
                  expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                  RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                  y = as.numeric(occMatrix[,myspecies]))

#add species data
listlengthDF$Species <- bugs.data$y

#check this is TRUE
all(row.names(occMatrix)==listlengthDF$visit)

#######################################################################################

# #set up the spline info
siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB,mtbqsDF$MTB)]
siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB,mtbqsDF$MTB)]

#centre and scale them?
#medX <- median(siteInfo$x)
#medY <- median(siteInfo$y)
#siteInfo$x <- (siteInfo$x - medX)/10000
#siteInfo$y <- (siteInfo$y - medY)/10000
#summary(siteInfo$x)
#summary(siteInfo$y)  
#save(siteInfo,file="siteInfo.RData") 

# #get info on whether species was seen at each site
speciesSite <- ddply(listlengthDF,.(siteIndex),summarise,
                     PA=max(Species,na.rm=T))
siteInfo$obs <- speciesSite$PA[match(siteInfo$siteIndex,speciesSite$siteIndex)]

# #fit normal gam
# library(mgcv)
# gam1 <- gam(obs ~ 1 + s(x,y), data=siteInfo, family=binomial)
# siteInfo$fits <- gam1$fitted.values
# library(ggplot2)
# qplot(x,y,data=siteInfo,colour=fits)+
#   scale_colour_gradient(low="blue",high="red")
# 
#######################################################################################

#get BUGS code using jagam
library(mgcv)
jags.ready <- jagam(obs ~ 1 + s(x,y,k=8), data=siteInfo, family=binomial,file="splines/jagam.txt")

#extract bits for the model
bugs.data$X = jags.ready$jags.data$X
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero
#bugs.data$nspline = length(bugs.data$zero)
#bugs.data$nspline1 = length(bugs.data$zero)+1
#bugs.data$nspline2 = (length(bugs.data$zero))*2

saveRDS(bugs.data,file="splines/bugs.data.rds")

#expand to include NA values ##########################################################

#take list of all MTBQ as the full data frame
siteInfo_NAs <- mtbqsDF
siteInfo_NAs <- subset(siteInfo_NAs,!duplicated(MTB))

#remove without county or badenwuttenberg?
#siteInfo_NAs <- subset(siteInfo_NAs,!is.na(Counties))
#siteInfo_NAs <- subset(siteInfo_NAs,Counties!="Baden-Württemberg")

#add on species observations
siteInfo_NAs$obs <- siteInfo$obs[match(siteInfo_NAs$MTB,siteInfo$MTB)]
table(siteInfo_NAs$obs)
sum(is.na(siteInfo_NAs$obs))/nrow(siteInfo_NAs)

#randomly insect 0 or 1 into the missing values
siteInfo_NAs$obs[is.na(siteInfo_NAs$obs)] <- sample(c(0,1),size=sum(is.na(siteInfo_NAs$obs)),replace=T)

#reorganise into same order as bugs.data object
siteInfo_NAs$siteIndex <- siteInfo$siteIndex[match(siteInfo_NAs$MTB,siteInfo$MTB)]
siteInfo_NAs <- arrange(siteInfo_NAs,siteIndex)
siteInfo_NAs$siteIndex <- 1:nrow(siteInfo_NAs)#to give values to those with NAs

saveRDS(siteInfo_NAs,file="splines/siteInfo_NAs.rds")

library(mgcv)
#scale coordinates
siteInfo_NAs$x_MTB <- (siteInfo_NAs$x_MTB-medX)/10000
siteInfo_NAs$y_MTB <- (siteInfo_NAs$y_MTB-medY)/10000

gam.ready <- gam(obs ~ 1 + s(x_MTB,y_MTB,k=15), 
                    data=siteInfo_NAs, family=binomial)
jags.ready <- jagam(obs ~ 1 + s(x_MTB,y_MTB,k=15), 
                    data=siteInfo_NAs, family=binomial,file="splines/jagam.txt")
#jags.ready <- jagam(obs ~ 1 + s(x, y, k = 60, bs = 'ds', m = c(1, 0.5)), 
#                    data=siteInfo_NAs, family=binomial,file="jagam.txt")

saveRDS(jags.ready,file="splines/jags.ready_NAs.rds")

bugs.data$nsite_Full <- nrow(siteInfo_NAs)
bugs.data$X_Full = jags.ready$jags.data$X
bugs.data$S1_Full = jags.ready$jags.data$S1
bugs.data$zero_Full = jags.ready$jags.data$zero

saveRDS(bugs.data,file="splines/bugs.data_NAs.rds")

# ##########################################################################################
# 
# set up space-time-spline info
siteInfo <- unique(listlengthDF[,c("MTB","siteIndex","yearIndex")])
siteInfo$x <- mtbqsDF$x_MTB[match(siteInfo$MTB,mtbqsDF$MTB)]
siteInfo$y <- mtbqsDF$y_MTB[match(siteInfo$MTB,mtbqsDF$MTB)]


#add species data
listlengthDF$Species <- as.numeric(occMatrix[,myspecies])

speciesSite <- ddply(listlengthDF,.(siteIndex,yearIndex),summarise,
                     obs = max(Species,na.rm=T),
                     nuVisits = length(visit),
                     nuSingles = sum(singleList))

siteInfo <- merge(siteInfo,speciesSite,by=c("siteIndex","yearIndex")) 

#fit normal gam
library(mgcv)
gam1 <- gam(obs ~ s(yearIndex) + s(x,y,yearIndex,k=50), #s(nuVisits,k=4) + s(nuSingles,k=4),
            data=siteInfo, family=binomial)
siteInfo$fits <- gam1$fitted.values

# library(ggplot2)
# qplot(x,y,data=siteInfo,colour=fits)+
#   scale_colour_viridis_c(option = "magma",direction = -1)+
#   facet_wrap(~yearIndex)+
#   theme_bw()

#predict model to whole range
library(tidyverse)
siteInfo_NAs <- mtbqsDF
siteInfo_NAs <- subset(siteInfo_NAs,!duplicated(MTB))
nuSites <- nrow(siteInfo_NAs)
siteInfo_NAs <- siteInfo_NAs %>% slice(rep(1:n(), each = length(unique(df$Year))))
siteInfo_NAs$yearIndex <- rep(1:length(unique(df$Year)),nuSites)
siteInfo_NAs$nuVisits <- 5
siteInfo_NAs$nuSingles <- 2
siteInfo_NAs$preds <- predict(gam1,siteInfo_NAs,type="response")

# qplot(x,y,data=siteInfo_NAs,colour=preds)+
#   facet_wrap(~yearIndex)+
#   scale_colour_viridis_c(option = "magma", direction = -1)+
#   theme_void()

#loop over all years
myYears <- sort(unique(df$Year))
for(i in 1:length(myYears)){
  
  qplot(x,y,data=subset(siteInfo_NAs,yearIndex==i),colour=preds)+
    scale_colour_viridis_c("Occupancy",
                           option = "magma",
                           limits=c(0,0.8),
                           #limits=c(0,0.45),
                           direction = -1)+
    theme_void()+
    ggtitle(myYears[i])
  
  ggsave(paste0("gifs/year",i,".png"),width=3,height=2.6)
}

### space time function ###############################################################

#for each species, predict 4 times
library(mgcv)
library(tidyverse)

fitSpeciesGAM <- function(myspecies){

# set up space-time-spline info
siteInfo <- unique(listlengthDF[,c("MTB","siteIndex","yearIndex")])
siteInfo$x <- mtbqsDF$x_MTB[match(siteInfo$MTB,mtbqsDF$MTB)]
siteInfo$y <- mtbqsDF$y_MTB[match(siteInfo$MTB,mtbqsDF$MTB)]
  
#add species data
listlengthDF$Species <- as.numeric(occMatrix[,myspecies])

speciesSite <- ddply(listlengthDF,.(siteIndex,yearIndex),summarise,
                     obs = max(Species,na.rm=T),
                     nuVisits = length(visit),
                     nuSingles = sum(singleList))
siteInfo <- merge(siteInfo,speciesSite,by=c("siteIndex","yearIndex")) 

#fit normal gam
gam1 <- gam(obs ~ s(yearIndex) + s(x,y,yearIndex) + 
              s(nuVisits,k=4) + 
              s(nuSingles,k=4),
            data=siteInfo, family=binomial)

#predict model to whole range
siteInfo_NAs <- mtbqsDF
siteInfo_NAs <- subset(siteInfo_NAs,!duplicated(MTB))
nuSites <- nrow(siteInfo_NAs)
siteInfo_NAs <- siteInfo_NAs %>% slice(rep(1:n(), each = length(unique(df$Year))))
siteInfo_NAs$yearIndex <- rep(1:length(unique(df$Year)),nuSites)
siteInfo_NAs$nuVisits <- 5
siteInfo_NAs$nuSingles <- 2
siteInfo_NAs$Species <- myspecies
siteInfo_NAs$preds <- predict(gam1,siteInfo_NAs,type="response")
return(siteInfo_NAs)

}


#run for all species
specieslist = read.delim("model-auxfiles/speciesTaskID_adult.txt")
output <- ldply(specieslist$Species,function(x)fitSpeciesGAM(x))
saveRDS(output,file="model-outputs/simpleGAMS_1990_2017.rds")

########################################################################################

#specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
#zst [is.infinite(zst)] <- 0
#inits <- function(){list(z = zst)}

#fill in the blanks more cleverly
zst [is.infinite(zst)] <- NA
replace_na_with_last<-function(x,a=!is.na(x)){
  x[which(a)[c(1,1:sum(a))][cumsum(a)+1]]
}

#inits <- function(){list(z = zst)}
for(i in 1:nrow(zst)){
  zst[i,] <- replace_na_with_last(zst[i,])
}

inits <- function(){list(z = zst)}
saveRDS(zst,file="splines/zst.rds")

# inits <- function(){list(z = zst,
#                          #state.a = runif(bugs.data$nstate,0.0001,0.01),
#                          #lphi = runif(bugs.data$nstate,0,0.01),
#                          #lgam = runif(bugs.data$nstate,0,0.01),
#                          effort.p = runif(1,0,0.01),
#                          mu.phenol = runif(1,-0.01,0.01),
#                          mu.phenol2 = runif(1,-0.01,0.01))}

### initial values including NAs
zst_NAs <- matrix(data=0,nrow=bugs.data$nsite,ncol=bugs.data$nyear)
zst_NAs[1:dim(zst)[1],1:dim(zst)[2]] <- zst
inits <- function(){list(z = zst)}

saveRDS(zst_NAs,file="splines/zst_NAs.rds")

### run model #########################################################################

#see HPC_static_jagam or
#HPC_tempcorr_jagam

########################################################################################

#get nationwide boxes
load("mtbqsDF.RData")

### static model #######################################################################

#running HPC_static_jagam 
out <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Static_spline/6992005/outSummary_static_spline.rds")
siteInfo <- arrange(siteInfo,siteIndex)
siteInfo$preds <- out[1:9699,1]

#add on x and y of the mtbq
siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]

ggplot(siteInfo)+
  geom_point(aes(x=x,y=y,color=preds))+
  scale_color_viridis_c()
#looks nice!!!

#static model that includes NAs
out <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Static_spline/6994377/outSummary_static_spline_NAs.rds")
siteInfo_NAs$preds <- out[1:9699,1]

#add on x and y of the mtbq
siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]

ggplot(siteInfo)+
  geom_point(aes(x=x,y=y,color=preds))+
  scale_color_viridis_c()
#looks nice!!!

### dynamic model #################################################################

#temp corr model
out <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/splines/outSummary_tempcorr_spline_NAs.rds")
#get first siteInfo object in script

#when it saved the summary
source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')
temp <- getBUGSfitsII(out,param="psi")
temp$MTB <- siteInfo_NAs$MTB[match(temp$Site,siteInfo_NAs$siteIndex)]
temp$x <- mtbqsDF$x[match(temp$MTB,mtbqsDF$MTB)]
temp$y <- mtbqsDF$y[match(temp$MTB,mtbqsDF$MTB)]

ggplot(temp)+
  geom_point(aes(x=x,y=y,color=mean))+
  scale_color_viridis_c()+
  facet_wrap(~Year)


#for jagsUI basic
out <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/splines/outSummary_tempcorr_spline_NAs.rds")

########################################################################################