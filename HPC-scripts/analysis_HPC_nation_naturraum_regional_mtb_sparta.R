#get all libraries we need
#suppressMessages()
suppressMessages(library(rjags))
suppressMessages(library(R2WinBUGS))
suppressMessages(library(jagsUI))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(gdata))

#load the relational table of task ids and species
speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_adult.txt"),as.is=T)
#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)]
#myspecies <- "Crocothemis erythraea"

#set stage
stage="adult"

#set seed
set.seed(3)

#number of MCMC samples
niterations = 50000

Sys.time()

#load in regional datasets
#myfiles <- list.files("derived-data")
myfiles <- list.files("/data/idiv_ess/Odonata")

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
  #out<-readRDS(paste("derived-data",x,sep="/"))
  out<-readRDS(paste("/data/idiv_ess/Odonata",x,sep="/"))
  out$File <- x
  return(out)
})

#extract state from file name
adultData$State <- sapply(adultData$File,function(x)strsplit(x,"\\.rds")[[1]][1])
adultData$State <- sapply(adultData$State,function(x)strsplit(x,"_")[[1]][3])
nrow(adultData)#1147558

##########################################################################################

#add gbif data to fill gaps
#gbifdata <- readRDS("derived-data/datafile_iNaturalist.rds")
gbifdata <- readRDS("/data/idiv_ess/Odonata/datafile_iNaturalist.rds")
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
#filter to 1990 onwards

df <- subset(adultData, Year>=1990  & Year<2018)

### sort MTBs #########################################################

#get MTB
df$MTB <- sapply(as.character(df$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,len-1)})

df <- subset(df,!MTB %in% c("6301","5548","5156","5056","4956","4455"))

######################################################################################

#get nationwide boxes
#load("mtbqsDF.RData")
load("/data/idiv_ess/Odonata/mtbqsDF.RData")
names(mtbqsDF)[which(names(mtbqsDF)=="Value")] <- "MTB"

all(df$MTB%in% mtbqsDF$MTB)
unique(df$MTB)[!unique(df$MTB) %in% mtbqsDF$MTB]

df$State <- mtbqsDF$Counties[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$State))

#check missing naturruam data
mtbqsDF$Natur[is.na(mtbqsDF$Natur)] <- mtbqsDF$MTB_Natur[is.na(mtbqsDF$Natur)]
df$Natur <- mtbqsDF$Natur[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$Natur))

mtbqsDF$CoarseNatur[is.na(mtbqsDF$CoarseNatur)] <- mtbqsDF$MTB_CoarseNatur[is.na(mtbqsDF$CoarseNatur)]
df$CoarseNatur <- mtbqsDF$CoarseNatur[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$CoarseNatur))

mtbqsDF$MTB_MidNatur <- gdata::trim(mtbqsDF$MTB_MidNatur)
df$MidNatur <- mtbqsDF$MTB_MidNatur[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$MidNatur))

##########################################################################################

#check we have data for all midnatur
length(unique(mtbqsDF$MTB_MidNatur))#85
summaryDF <- ddply(df,.(MidNatur),summarise,nuRecs=length(MidNatur),nuYears=length(unique(Year)))
nrow(summaryDF)
#we have data from all 85
#and pretty good time frames

# ####################################################################################

#subset by average phenology across whole germany
#dfS <- subset(df, Species==myspecies)

#same for all species
# obsPhenolData <- summarise(df,
#                            minDay = round(quantile(yday,0.05)),
#                            maxDay = round(quantile(yday,0.95)))
# df <- subset(df, yday > obsPhenolData$minDay & yday < obsPhenolData$maxDay)

######################################################################################

#remove sites visited once
#siteSummary <- ddply(df,.(MTB_Q),summarise,nuYears=length(unique(Year)))
#df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])

#####################################################################################

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

# #get MTB
# listlengthDF$MTB <- sapply(as.character(listlengthDF$MTB_Q),function(x){
#   len <- nchar(x)
#   substr(x,1,(len-1))})

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

#MTB as the site index
listlengthDF$siteIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                   listlengthDF$MTB)))
listlengthDF$mtbIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                  listlengthDF$MTB)))#redundant

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB,mtbqsDF$MTB)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
#subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB,mtbqsDF$MTB)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
#subset(listlengthDF,is.na(Naturraum))

listlengthDF$MidNaturraum <- mtbqsDF$MTB_MidNatur[match(listlengthDF$MTB,mtbqsDF$MTB)]
listlengthDF$MidNaturraum[listlengthDF$Naturraum=="NeiÃ¡egebiet"] <- "OBERLAUSITZ"
listlengthDF$mnIndex <- as.numeric(factor(listlengthDF$MidNaturraum))
subset(listlengthDF,is.na(MidNaturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("mtbIndex","siteIndex","MTB","nnIndex","cnIndex","mnIndex","MidNaturraum")])
head(siteInfo)
#saveRDS(siteInfo,file="siteInfo_midnaturraum.rds")

#######################################################################################

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)
siteInfo <- arrange(siteInfo,siteIndex)

########################################################################################

#fit nation-wide model with random slopes to each box

########################################################################################

#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nraum = length(unique(siteInfo$nnIndex)),
                  nmraum = length(unique(siteInfo$mnIndex)),
                  ncraum = length(unique(siteInfo$cnIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  raum = listlengthDF$nnIndex,
                  mraum = listlengthDF$mnIndex,
                  craum = listlengthDF$cnIndex,
                  year = listlengthDF$yearIndex,
                  craumS = siteInfo$cnIndex,
                  raumS = siteInfo$nnIndex,
                  mraumS = siteInfo$mnIndex,
                  yday = listlengthDF$yday - median(listlengthDF$yday),
                  yday2 = listlengthDF$yday^2 - median(listlengthDF$yday^2),
                  nuSpecies = log(listlengthDF$nuSpecies) - median(log(listlengthDF$nuSpecies)),
                  singleList = listlengthDF$singleList,
                  shortList = listlengthDF$shortList,
                  nuRecs = log(listlengthDF$nuRecords) - median(log(listlengthDF$nuRecords)),
                  nuSS = log(listlengthDF$samplingSites) - median(log(listlengthDF$samplingSites)),# up to 3
                  expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                  RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                  y = as.numeric(occMatrix[,myspecies]))

listlengthDF$Species <- bugs.data$y
all(row.names(occMatrix)==listlengthDF$visit)

#the below are used the linear regression model in the model file -see below
bugs.data$sumX <- sum(1:bugs.data$nyear)
bugs.data$sumX2 <- sum((1:bugs.data$nyear)^2)

#########################################################################

#index each MTB to naturraum
StrIdx <- array(data=0, dim = c(bugs.data$nsite,
                                bugs.data$nyear,
                                bugs.data$ncraum))
for(i in 1:bugs.data$nsite){
  StrIdx[i,,bugs.data$craumS[i]] <- 1
}
bugs.data$crIdx <- StrIdx

#number of sites per coarseroam
nsite_cr <- ddply(siteInfo,.(cnIndex),summarise,nuSites = length(unique(siteIndex)))
bugs.data$nsite_cr <- nsite_cr$nuSites         

#index each MTB to midnaturraum
StrIdx <- array(data=0, dim = c(bugs.data$nsite,
                                bugs.data$nyear,
                                bugs.data$nmraum))
for(i in 1:bugs.data$nsite){
  StrIdx[i,,bugs.data$mraumS[i]] <- 1
}
bugs.data$mrIdx <- StrIdx

#number of sites per coarseroam
nsite_mr <- ddply(siteInfo,.(mnIndex),summarise,nuSites = length(unique(siteIndex)))
bugs.data$nsite_mr <- nsite_mr$nuSites  

#########################################################################

#specify initial values
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
zst [is.infinite(zst)] <- NA

#fill in the blanks more cleverly
replace_na_with_last<-function(x,a=!is.na(x)){
  x[which(a)[c(1,1:sum(a))][cumsum(a)+1]]
}

#inits <- function(){list(z = zst)}
for(i in 1:nrow(zst)){
  zst[i,] <- replace_na_with_last(zst[i,])
}

inits <- function(){list(z = zst)}

# inits <- function(){list(z = zst,
#                          mu.p = runif(1,-2,2),
#                          single.p = runif(1,0,0.001),
#                          effort.p = runif(1,0,0.001),
#                          mu.phenol = runif(1,-0.001,0.001),
#                          mu.phenol2 = runif(1,-0.001,0.001))}

########################################################################################

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
#n.cores = 3
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 

#############################################################################

#choose model file
#modelfile="/data/idiv_ess/Odonata/BUGS_sparta_nation_naturraum.txt"
#modelfile="/data/idiv_ess/Odonata/BUGS_sparta_regional_midnaturraum.txt"
#modelfile="/data/idiv_ess/Odonata/BUGS_sparta_regional_midnaturraumtrends.txt"
modelfile="/data/idiv_ess/Odonata/BUGS_sparta_regional_naturraum.txt"

effort = "shortList"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("mean.p","mup","cr.a","muZ.cra","regres.psi")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=20,
            n.chains=n.cores, n.burnin=round(niterations/2),
            n.iter=niterations,parallel=T)
Sys.time()

#save as output file - for regional/dynamic model
saveRDS(out,file=paste0("out_sparta_regional_nation_naturraum_",stage,"_",myspecies,".rds"))
#saveRDS(out,file=paste0("out_sparta_regional_nation_midnaturraumtrends_",stage,"_",myspecies,".rds"))
#saveRDS(out,file=paste0("out_sparta_regional_nation_finenaturraum_",stage,"_",myspecies,".rds"))

########################################################################################




