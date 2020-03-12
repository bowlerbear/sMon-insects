#SELECT:

#decide on stage
stage="adult"
#stage="juv"

#decide on state
#state="Sa"
#state="Bav"
state="MV"
#state = "Sax"
#state="SH"

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
df <- subset(df,!(month==1&day==1))

#missing any data with missing values
df <- subset(df,!is.na(day))
df <- subset(df,!is.na(Year))
df <- subset(df,!is.na(month))

#filter months and years
yearDF <- read.delim(paste0("/data/idiv_ess/Odonata/yearDF_",stage,"_",state,".txt"),as.is=T)
df <- subset(df, Year>=(min(yearDF$Year)) & Year<2017)

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
               RpS = length(Species)/length(unique(Species)),
               expertise = sum(Expertise),
               samplingSites = length(unique(interaction(lat,lon))))
  
  #add on some indices
  out$siteIndex <- as.numeric(factor(out$MTB_Q))
  out$yearIndex <- as.numeric(factor(out$Year))
  
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

###########################################################################################
#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                    nyear = length(unique(listlengthDF$yearIndex)),
                    nvisit = nrow(listlengthDF),
                    site = listlengthDF$siteIndex,
                    year = listlengthDF$yearIndex,
                    yday = listlengthDF$day - median(listlengthDF$day),
                    nuSpecies = log(listlengthDF$nuSpecies) - log(median(listlengthDF$nuSpecies)),
                    singleList = listlengthDF$singleList,
                    shortList = listlengthDF$shortList,
                    nuRecs = log(listlengthDF$nuRecords) - log(median(listlengthDF$nuRecords)),
                    nuSS = log(listlengthDF$samplingSites) - log(median(listlengthDF$samplingSites)),# up to 3
                    expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                    RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                    y = as.numeric(occMatrix[,myspecies]))
listlengthDF$Species <- bugs.data$y
  
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
#effort (1): linear length and single list
bugs.data$Effort <- bugs.data[["nuSpecies"]]
modelfile="/data/idiv_ess/Odonata/BUGS_sparta.txt"

#specify parameters to monitor
params <- c("phenol.p","phenol2.p","effort.p","psi.fs","single.p")
  
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
               n.chains=n.cores, n.burnin=40000,n.iter=200000,parallel=T)
  
#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_nuSpecies_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################
#effort (2): short and single list
bugs.data$Effort <- bugs.data[["shortList"]]
modelfile="/data/idiv_ess/Odonata/BUGS_sparta.txt"

#specify parameters to monitor
params <- c("phenol.p","phenol2.p","effort.p","psi.fs","single.p")

#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
            n.chains=n.cores, n.burnin=40000,n.iter=200000,parallel=T)

#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_speciesList_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################
#effort (3): linear length and single list, and expertise
bugs.data$Effort <- bugs.data[["nuSpecies"]]
bugs.data$Effort2 <- bugs.data[["expertise"]]

modelfile="/data/idiv_ess/Odonata/BUGS_sparta_2effort.txt"

#specify parameters to monitor
params <- c("phenol.p","phenol2.p","effort.p","psi.fs","single.p")

#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
            n.chains=n.cores, n.burnin=40000,n.iter=200000,parallel=T)

#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_expertise_",stage,"_",state,"_", myspecies,".rds"))

############################################################################################
#effort (4): linear length and single list, and number of records
bugs.data$Effort <- bugs.data[["shortList"]]
bugs.data$Effort2 <- bugs.data[["nuRecs"]]

modelfile="/data/idiv_ess/Odonata/BUGS_sparta_2effort.txt"

#specify parameters to monitor
params <- c("phenol.p","phenol2.p","effort.p","effort2.p","psi.fs","single.p")

#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
            n.chains=n.cores, n.burnin=40000,n.iter=200000,parallel=T)

#save as output file
saveRDS(data.frame(out$summary),file=paste0("outSummary_nuRecs_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################
