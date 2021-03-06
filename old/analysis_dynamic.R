########################################################################################

#load in regional datasets
myfiles <- list.files("derived-data")


#read in and combine all adult files
adultFiles<-myfiles[grepl("adult_datafile",myfiles)]
adultFiles<-adultFiles[grepl("rds",adultFiles)]

#combine these files
library(plyr)
adultData <- ldply(adultFiles,function(x){
  out<-readRDS(paste("derived-data",x,sep="/"))
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

######################################################################################

#pick species
myspecies="Aeshna cyanea"

######################################################################################

#subset to phenology by regions

#read in and combine all phenology files
phenolFiles<-list.files()[grepl("speciesDays",list.files())]
phenolData <- ldply(phenolFiles,function(x){
  out<-read.delim(x)
  out$File <- x
  return(out)
})

#extract state from file name
phenolData$State <- sapply(phenolData$File,function(x)strsplit(x,"\\.txt")[[1]][1])
phenolData$State <- sapply(phenolData$State,function(x)strsplit(x,"_")[[1]][3])
phenolData <- subset(phenolData, Species==myspecies)

df <- subset(df,interaction(yday,State) %in% interaction(phenolData$day,phenolData$State))

######################################################################################

#define a visit
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
load("mtbqsDF.RData")
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
df <- subset(df, Species==myspecies)
df <- merge(df,mtbqsDF,by="MTB_Q",all.x=T)
boxData <- ddply(df, .(km50),summarise,
                 nuYears = length(unique(Year)),
                 nuRecs = length(Species))

#plotting
library(raster)
library(sp)
germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
myGrid <- raster('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/km50grid.tif')
germanAdmin <- spTransform(germanAdmin,projection(myGrid))
myGrid[] <- NA
myGrid[boxData$km50] <- boxData$nuYears
plot(myGrid)
plot(germanAdmin,add=T)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("stateIndex","siteIndex","boxIndex","km50")])
siteInfo <- arrange(siteInfo,stateIndex,boxIndex,siteIndex)
head(siteInfo)

########################################################################################

#fit dynamic occupancy model

########################################################################################

#get model output
library(ggplot2)
modelSummary <- readRDS("model-outputs/modelSummary_dynamic_Odonata_adult_Sa.rds")

#get species names
modelSummary$Species <- sapply(modelSummary$File, function(x)strsplit(x,"_")[[1]][5])
modelSummary$Species <- gsub(".rds","",modelSummary$Species)

#look at overall occupancy
modelSummary_psi.fs <- subset(modelSummary, grepl("psi.fs",modelSummary$Param)) 
modelSummary_psi.fs$Index <- sapply(modelSummary_psi.fs$Param,function(x){sub(".*\\[([^][]+)].*", "\\1", x)})
modelSummary_psi.fs$Index <- as.numeric(modelSummary_psi.fs$Index)
ggplot(modelSummary_psi.fs)+
  geom_point(aes(x=Index,y=mean))+
  geom_ribbon(aes(x=Index,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Species)

#colonisation
modelSummary_psi.fs <- subset(modelSummary, grepl("colonize",modelSummary$Param)) 
modelSummary_psi.fs$Index <- sapply(modelSummary_psi.fs$Param,function(x){sub(".*\\[([^][]+)].*", "\\1", x)})
modelSummary_psi.fs$Index <- as.numeric(modelSummary_psi.fs$Index)
ggplot(modelSummary_psi.fs)+
  geom_point(aes(x=Index,y=mean))+
  #geom_ribbon(aes(x=Index,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Species,scales="free")

#persistence
modelSummary_psi.fs <- subset(modelSummary, grepl("persist",modelSummary$Param)) 
modelSummary_psi.fs$Index <- sapply(modelSummary_psi.fs$Param,function(x){sub(".*\\[([^][]+)].*", "\\1", x)})
modelSummary_psi.fs$Index <- as.numeric(modelSummary_psi.fs$Index)
ggplot(modelSummary_psi.fs)+
  geom_point(aes(x=Index,y=mean))+
  #geom_ribbon(aes(x=Index,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Species,scales="free")

########################################################################################

#compare with standard model
library(ggplot2)
modelSummary <- readRDS("model-outputs/Odonata_modelSummary_flightperiod/5098409/modelSummary_Odonata_adult_Sa.rds")

#get species names
modelSummary$Species <- sapply(modelSummary$File, function(x)strsplit(x,"_")[[1]][6])
modelSummary$Species <- gsub(".rds","",modelSummary$Species)

#look at overall occupancy
modelSummary_psi.fs <- subset(modelSummary, grepl("psi.fs",modelSummary$Param)) 
modelSummary_psi.fs$Index <- sapply(modelSummary_psi.fs$Param,function(x){sub(".*\\[([^][]+)].*", "\\1", x)})
modelSummary_psi.fs$Index <- as.numeric(modelSummary_psi.fs$Index)
ggplot(modelSummary_psi.fs)+
  geom_point(aes(x=Index,y=mean))+
  geom_ribbon(aes(x=Index,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Species)

#similar!!
########################################################################################
