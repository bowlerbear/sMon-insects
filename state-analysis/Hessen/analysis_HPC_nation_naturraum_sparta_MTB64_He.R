#get all libraries we need
#suppressMessages()
suppressMessages(library(rjags))
suppressMessages(library(R2WinBUGS))
suppressMessages(library(jagsUI))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))


#load the relational table of task ids and species
speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_adult.txt"),as.is=T)
#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")) 
#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 

#set stage
stage="adult"

#set seed
set.seed(3)

#number of MCMC samples
niterations = 50000

Sys.time()

#load in regional datasets
#those that are updated and in 64 units
adultFiles <- "adult_datafile_He_updated_MTB64.rds"

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
nrow(adultData)#110392

##########################################################################################

#add gbif data to fill gaps
#gbifdata <- readRDS("derived-data/datafile_GBIF.rds")
#gbifdata <- readRDS("/data/idiv_ess/Odonata/datafile_GBIF.rds")
#nrow(gbifdata)#38191

#combine the two
#adultData <- rbind(adultData,gbifdata)

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

#####################################################################################

#remove MTBQs not in the shapefile = "51561" "50561" "49563"
#remove MTBQS without naurraum 81443, 44553, 58401
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)
adultData <- subset(adultData, !MTB_Q %in% c("51561","50561","49563","55484","63012",
                                             "81443","44553","58401"))

###################################################################################

#filter to 1980 onwards

#df <- subset(adultData, Year>=1980  & Year<2017)
df <- subset(adultData, Year>=2008  & Year<2020)

####################################################################################

#get nationwide boxes
#load("mtbqsDF.RData")
load("/data/idiv_ess/Odonata/mtbqsDF.RData")
all(df$MTB_Q%in% mtbqsDF$MTB_Q)
#unique(df$MTB_Q[!df$MTB_Q %in% mtbqsDF$MTB_Q])

df$State <- mtbqsDF$Counties[match(df$MTB_Q,mtbqsDF$MTB_Q)]
sum(is.na(df$State))

#check missing naturruam data
mtbqsDF$Natur[is.na(mtbqsDF$Natur)] <- mtbqsDF$MTB_Natur[is.na(mtbqsDF$Natur)]
df$Natur <- mtbqsDF$Natur[match(df$MTB_Q,mtbqsDF$MTB_Q)]
sum(is.na(df$Natur))

mtbqsDF$CoarseNatur[is.na(mtbqsDF$CoarseNatur)] <- mtbqsDF$MTB_CoarseNatur[is.na(mtbqsDF$CoarseNatur)]
df$CoarseNatur <- mtbqsDF$CoarseNatur[match(df$MTB_Q,mtbqsDF$MTB_Q)]
sum(is.na(df$CoarseNatur))

mtbqsDF$MTB_MidNatur <- gdata::trim(mtbqsDF$MTB_MidNatur)
df$MidNatur <- mtbqsDF$MTB_MidNatur[match(df$MTB_Q,mtbqsDF$MTB_Q)]
sum(is.na(df$MidNatur))

#####################################################################################

#subset by average phenology across whole germany
dfS <- subset(df, Species==myspecies)
obsPhenolData <- summarise(dfS,
                           minDay = round(quantile(yday,0.05)),
                           maxDay = round(quantile(yday,0.95)))
df <- subset(df, yday > obsPhenolData$minDay & yday < obsPhenolData$maxDay)

######################################################################################

#remove sites visited once
siteSummary <- ddply(df,.(MTB_Q),summarise,nuYears=length(unique(Year)))
df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])

#####################################################################################

#define a visit
df$visit <- paste(df$MTB64,df$Date,df$Beobachter,sep="_")

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
  out <- ddply(df,.(visit,Date,MTB,MTB_Q,MTB64),summarise,
               nuSpecies=length(unique(Species)),
               nuRecords=length(Species))
  
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

#site is the MTB64
listlengthDF$siteIndex <- as.numeric(factor(listlengthDF$MTB64))

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
#subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
#subset(listlengthDF,is.na(Naturraum))

listlengthDF$MidNaturraum <- mtbqsDF$MTB_MidNatur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$mnIndex <- as.numeric(factor(listlengthDF$MidNaturraum))
#subset(listlengthDF,is.na(MidNaturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("siteIndex","nnIndex","cnIndex","mnIndex")])

#######################################################################################

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)
siteInfo <- arrange(siteInfo,siteIndex)

raumInfo <- unique(siteInfo[,c("nnIndex","cnIndex")])

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
                  y = as.numeric(occMatrix[,myspecies]))

listlengthDF$Species <- bugs.data$y

all(row.names(occMatrix)==listlengthDF$visit)

#the below are used the linear regression model in the model file -see below
bugs.data$sumX <- sum(1:bugs.data$nyear)
bugs.data$sumX2 <- sum((1:bugs.data$nyear)^2)

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

########################################################################################

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
#n.cores = 3
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 

###########################################################################################

#choose model file
modelfile="/data/idiv_ess/Odonata/BUGS_sparta_bundesland_naturraum.txt"

effort = "shortList"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("psi.fs","regres.psi","mean.p","mup")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=5,
            n.chains=n.cores, n.burnin=round(niterations/2),
            n.iter=niterations,parallel=T)

Sys.time()

#save as output file 
saveRDS(out,file=paste0("out_sparta_He_MTB64_naturraum_",stage,"_",myspecies,".rds"))

########################################################################################




