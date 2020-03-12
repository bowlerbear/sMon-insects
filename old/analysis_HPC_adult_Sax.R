#SELECT:

#decide on stage
stage="adult"
#stage="juv"

#decide on state
#state="Sa"
#state="Bav"
#state="NRW"
state = "Sax"
#state="SH"

#effort term
effort="nuSpecies"

####################################################################################################
#load data frame

df <- readRDS(paste0("/data/idiv_ess/Odonata/",stage,"_datafile_",state,".rds"))

#load the relational table of task ids and species
speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_",stage,"_",state,".txt"),as.is=T)
#species are labelled 1:length(species)

#get task id
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 

#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 

#format dates
library(lubridate)
df$Date <- as.Date(df$Date,file="%Y-%m-%d")
df$month <- month(df$Date)
df$day <- yday(df$Date)
df$Year <- year(df$Date)

#missing any data with missing values
df <- subset(df,!is.na(day))
df <- subset(df,!is.na(Year))
df <- subset(df,!is.na(month))

#filter months and years
yearDF <- read.delim(paste0("/data/idiv_ess/Odonata/yearDF_",stage,"_",state,".txt"),as.is=T)
df <- subset(df, Year>=(yearDF$Year) & Year<2017)

#restrict data differently by stage
#if(stage=="juv"){
#  df <- subset(df, month %in% c(4:8))
#}else if(stage=="adult"){
#  df <- subset(df, month %in% c(5:9))
#}

#get phenology info
speciesDays <- read.delim(paste0("/data/idiv_ess/Odonata/speciesDays_",stage,"_",state,".txt"),as.is=T)
speciesDays <- subset(speciesDays,Species==myspecies)

#remove surveys that were outside the active pheneology period of the species
df <- subset(df, day %in% speciesDays$day)

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
  out <- ddply(df,.(visit,Date,Year,month,day,MTB_Q),summarise,
               nuSpecies=length(unique(Species)),
               nuRecords=length(Species),
               Richness2=mean(Richness),
               RpS = length(Species)/Richness2,
               expertise = median(Expertise),
               samplingSites = length(unique(interaction(lon,lat))))
  
  #add on some indices
  out$siteIndex <- as.numeric(factor(out$MTB_Q))
  out$yearIndex <- as.numeric(factor(out$Year))
  
  #sort dataset to match with the occurrence Matrix
  out <- arrange (out,visit)
}
listlengthDF <- getListLength(df)

#rows of occuMatrix match visits
#all(listlengthDF$visit==row.names(occMatrix))

#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                    nyear = length(unique(listlengthDF$yearIndex)),
                    nvisit = nrow(listlengthDF),
                    site = listlengthDF$siteIndex,
                    year = listlengthDF$yearIndex,
                    yday = listlengthDF$day - median(listlengthDF$day),
                    singleList = ifelse(listlengthDF$nuSpecies>1,0,1),
                    nuSpecies = log(listlengthDF$nuSpecies) - log(median(listlengthDF$nuSpecies)),
                    nuRecs = log(listlengthDF$nuRecords) - log(median(listlengthDF$nuRecords)),
                    nuSS = listlengthDF$samplingSites - median(listlengthDF$samplingSites),# up to 3
                    expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                    RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                    y = as.numeric(occMatrix[,myspecies]))
listlengthDF$Species <- bugs.data$y
  
#specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}
  
#define model params
modelfile="/data/idiv_ess/Odonata/BUGS_sparta.txt"
ni <- 6000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3

#fit model
library(rjags)
library(R2WinBUGS)
library(jagsUI)

#get core info
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")
bugs.data$Effort <- bugs.data[[effort]]
bugs.data$Effort_v2 <- bugs.data$RpS
  
#specify parameters to monitor
params <- c("phenol.p","phenol2.p","effort.p","psi.fs","single.p")
  
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
               n.chains=n.cores, n.burnin=40000,n.iter=200000,parallel=T)
  
#save as output file
saveRDS(out,file=paste0("out_flightperiod_",effort,"_",stage,"_",state,"_", myspecies,".rds"))
