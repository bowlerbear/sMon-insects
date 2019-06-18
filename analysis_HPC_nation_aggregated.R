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
  out<-readRDS(paste("derived-data",x,sep="/"))
  #out<-readRDS(paste("/data/idiv_ess/Odonata",x,sep="/"))
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
#df$visit <- paste(df$MTB_Q,df$Date,df$Beobachter,sep="_")
df$visit <- paste(df$MTB_Q,df$Date,sep="_")

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

#rows of occuMatrix match visits
#all(listlengthDF$visit==row.names(occMatrix))

#add on some indices
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- year(listlengthDF$Date)
listlengthDF$yday <- yday(listlengthDF$Date)
listlengthDF$State <- df$State[match(listlengthDF$MTB_Q,df$MTB_Q)]
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))
listlengthDF$stateIndex <- as.numeric(factor(listlengthDF$State))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

##########################################################################################

#get nationwide boxes
#load("mtbqsDF.RData")
load("/data/idiv_ess/Odonata/mtbqsDF.RData")
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

# Examine amount of data per 50 km box
# df <- subset(df, Species==myspecies)
# df <- merge(df,mtbqsDF,by="MTB_Q",all.x=T)
# boxData <- ddply(df, .(km50),summarise,
#                  nuYears = length(unique(Year)),
#                  nuRecs = length(Species))
# 
# #plotting
# library(raster)
# library(sp)
# germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
# myGrid <- raster('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/km50grid.tif')
# germanAdmin <- spTransform(germanAdmin,projection(myGrid))
# myGrid[] <- NA
# myGrid[boxData$km50] <- boxData$nuYears
# plot(myGrid)
# #myGrid[boxData$km50] <- boxData$nuRecs
# #plot(myGrid)
# plot(germanAdmin,add=T)

#######################################################################################

#get coordinates of mtb boxes

#get MTB Q
#MTBshapefile <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",layer="MTBQ_25833")
#crs specified as +proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs
#MTBshapefile@data$Q <- NA
#MTBshapefile@data$Q[which(MTBshapefile@data$Quadrant=="NW")]<-1
#MTBshapefile@data$Q[which(MTBshapefile@data$Quadrant=="NO")]<-2
#MTBshapefile@data$Q[which(MTBshapefile@data$Quadrant=="SW")]<-3
#MTBshapefile@data$Q[which(MTBshapefile@data$Quadrant=="SO")]<-4
#MTBshapefile@data$MTBQ <- paste0(as.character(MTBshapefile@data$Value),
#                                 as.character(MTBshapefile@data$Q))

#make spline code

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("stateIndex","siteIndex","boxIndex")])
siteInfo <- arrange(siteInfo,stateIndex,boxIndex,siteIndex)
head(siteInfo)

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

#all(row.names(occMatrix)==listlengthDF$visit)

#######################################################################################

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

modelfile="/data/idiv_ess/Odonata/BUGS_sparta_nation.txt"
#modelfile="R/BUGS_sparta_nation.txt"

effort = "nuSpecies"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("int","box.a.effect","effort.p","single.p")

#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
            n.chains=3, n.burnin=10000,n.iter=50000,parallel=T)

#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_nation_",stage,"_",myspecies,".rds"))

########################################################################################