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
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))
#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)]

#set stage
stage="adult"

#set seed
set.seed(3)

#number of MCMC samples
niterations = 40000

Sys.time()

#load in regional datasets
#myfiles <- list.files("derived-data")
myfiles <- list.files("/data/idiv_ess/Odonata")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

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
nrow(adultData)#1023689

##########################################################################################

#add gbif data to fill gaps
#gbifdata <- readRDS("derived-data/datafile_GBIF.rds")
gbifdata <- readRDS("/data/idiv_ess/Odonata/datafile_GBIF.rds")
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

df <- subset(adultData, Year>=1990  & Year<2017)

###remove missing mtbs #########################################################

#get MTB
df$MTB <- sapply(as.character(df$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,len-1)})

df <- subset(df, !MTB %in% c("6301","5548","5156","5056","4956",
                             "6644","8525","8144","4602","6204","4455","5840"))

######################################################################################

#get nationwide boxes
#load("mtbqsDF.RData")
load("/data/idiv_ess/Odonata/mtbqsDF.RData")

names(mtbqsDF)[which(names(mtbqsDF)=="Value")] <- "MTB"

all(df$MTB%in% mtbqsDF$MTB)

df$State <- mtbqsDF$Counties[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$State))

#check missing naturruam data
mtbqsDF$Natur[is.na(mtbqsDF$Natur)] <- mtbqsDF$MTB_Natur[is.na(mtbqsDF$Natur)]
df$Natur <- mtbqsDF$Natur[match(df$MTB,mtbqsDF$MTB)]
df$Natur <- trim(df$Natur)
sum(is.na(df$Natur))

mtbqsDF$CoarseNatur[is.na(mtbqsDF$CoarseNatur)] <- mtbqsDF$MTB_CoarseNatur[is.na(mtbqsDF$CoarseNatur)]
df$CoarseNatur <- mtbqsDF$CoarseNatur[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$CoarseNatur))

mtbqsDF$MTB_MidNatur <- trim(mtbqsDF$MTB_MidNatur)
df$MidNatur <- mtbqsDF$MTB_MidNatur[match(df$MTB,mtbqsDF$MTB)]
sum(is.na(df$MidNatur))
#1
#OBERLAUSITZ
df$MidNatur[is.na(df$MidNatur)] <- "OBERLAUSITZ"
df$MidNatur <- trim(df$MidNatur)
mtbqsDF$MTB_MidNatur[mtbqsDF$MTB_Natur=="NeiÃ¡egebiet"] <- "OBERLAUSITZ"

##########################################################################################

#check we have data for all midnatur
summaryDF <- ddply(df,.(MidNatur),summarise,nuRecs=length(MidNatur),nuYears=length(unique(Year)))
#we have data from all 85
#and pretty good time frames

####################################################################################

#subset by average phenology across whole germany
#dfS <- subset(df, Species==myspecies)

#same for all species
# obsPhenolData <- summarise(df,
#                            minDay = round(quantile(yday,0.05)),
#                            maxDay = round(quantile(yday,0.95)))
# df <- subset(df, yday > obsPhenolData$minDay & yday < obsPhenolData$maxDay)

#####################################################################################

#remove sites visited once
#siteSummary <- ddply(df,.(MTB_Q),summarise,nuYears=length(unique(Year)))
#df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])

#####################################################################################

#define a visit
df$visit <- paste(df$Natur,df$Date,df$Beobachter,sep="_")

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
  out <- ddply(df,.(visit,Date,Natur),summarise,
               nuSpecies=length(unique(Species)),
               nuRecords=length(Species),
               nuMTBQs=length(unique(MTB_Q)),
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
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- year(listlengthDF$Date)
listlengthDF$yday <- yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))

#add natur raum indices
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$Natur,mtbqsDF$Natur)]
subset(listlengthDF,is.na(CoarseNaturraum))
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))

listlengthDF$MidNaturraum <- mtbqsDF$MTB_MidNatur[match(listlengthDF$Natur,mtbqsDF$Natur)]
subset(listlengthDF,is.na(MidNaturraum))
listlengthDF$mnIndex <- as.numeric(factor(listlengthDF$MidNaturraum))

listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Natur))

#MTB as the site index
listlengthDF$siteIndex <- listlengthDF$nnIndex

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("siteIndex","MidNaturraum","Natur",
                                   "nnIndex","cnIndex","mnIndex")])
head(siteInfo)
max(siteInfo$siteIndex)==nrow(siteInfo)

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
                  nmraum = length(unique(siteInfo$mnIndex)),
                  ncraum = length(unique(siteInfo$cnIndex)),
                  nraum = length(unique(siteInfo$nnIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  raum = listlengthDF$nnIndex,
                  mraum = listlengthDF$mnIndex,
                  craum = listlengthDF$cnIndex,
                  year = listlengthDF$yearIndex,
                  raumS = siteInfo$nnIndex,
                  craumS = siteInfo$cnIndex,
                  mraumS = siteInfo$mnIndex,
                  yday = listlengthDF$yday - median(listlengthDF$yday),
                  yday2 = listlengthDF$yday^2 - median(listlengthDF$yday^2),
                  nuSpecies = log(listlengthDF$nuSpecies) - median(log(listlengthDF$nuSpecies)),
                  singleList = listlengthDF$singleList,
                  shortList = listlengthDF$shortList,
                  nuMTBQs = log(listlengthDF$nuMTBQs) - median(log(listlengthDF$nuMTBQs)),
                  nuRecs = log(listlengthDF$nuRecords) - median(log(listlengthDF$nuRecords)),
                  nuSS = log(listlengthDF$samplingSites) - median(log(listlengthDF$samplingSites)),# up to 3
                  expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                  RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                  y = as.numeric(occMatrix[,myspecies]))

listlengthDF$Species <- bugs.data$y
all(row.names(occMatrix)==listlengthDF$visit)

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
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

#############################################################################

#choose model file
modelfile="/data/idiv_ess/Odonata/BUGS_sparta_regional_nation_midnaturraum.txt"
#modelfile="R/BUGS_sparta_regional_nation_midnaturraum.txt"

effort = "shortList"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("mean.p","mup","muZ")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=10,
            n.chains=n.cores, n.burnin=round(niterations/2),
            n.iter=niterations,parallel=T)

Sys.time()

#save as output file - for regional/dynamic model
saveRDS(out,file=paste0("out_sparta_regional_nation_raumunit_",stage,"_",myspecies,".rds"))

########################################################################################