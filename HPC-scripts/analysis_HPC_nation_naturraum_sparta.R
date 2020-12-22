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

speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/problemspeciesTaskID_adult.txt"),as.is=T)


#get task id
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 
#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 

#set stage
stage="adult"

#set seed
set.seed(3)

#number of MCMC samples
niterations = 5000

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

#####################################################################################

#tidy and subset MTBQs

adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)
#adultData <- subset(adultData, !MTB_Q %in% c("51561","50561","49563","55484","63012",
#                                             "81443","44553","58401"))

adultData <- subset(adultData, !MTB_Q %in% c("63012","55484","51561","50561","49563","44553"))

###################################################################################
#filter to 1980 onwards

df <- subset(adultData, Year>=1980  & Year<2017)
#df <- subset(adultData, Year>=1980  & Year<2020)

#and subset to Hessen
#df <- subset(df, File=="adult_datafile_He_updated.rds")

#write list of species###############################################################

#summaryInfo <- ddply(df, .(Species), summarise, nuRecs=length(Species)) 
#run model for species with at least 50 records
#subset(summaryInfo,nuRecs<50)
#Species nuRecs
#38    Gomphus simillimus     17
#44    Lestes macrostigma      3
#57 Onychogomphus uncatus      1

#speciesList <- summaryInfo$Species[summaryInfo$nuRecs>=50]
#speciesDF <- data.frame(Species=speciesList,
#                        TaskID=1:length(speciesList))

#write.table(speciesDF,file=paste0("speciesTaskID_","adult",".txt"),
#            sep="\t",row.names=FALSE)

######################################################################################

#pick species

#practise species
#myspecies="Aeshna cyanea"
#stage="adult"

######################################################################################

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

##########################################################################################

#subset to phenology by regions

#read in and combine all phenology files
#phenolFiles<-list.files()[grepl("speciesDays",list.files())]
#phenolData <- ldply(phenolFiles,function(x){
#  out<-read.delim(x)
#  out$File <- x
#  return(out)
#})

# phenolFiles<-list.files("/data/idiv_ess/Odonata")[grepl("speciesDays",list.files("/data/idiv_ess/Odonata"))]
# phenolData <- ldply(phenolFiles,function(x){
#   out<-read.delim(paste("/data/idiv_ess/Odonata",x,sep="/"))
#   out$File <- x
#   return(out)
# })
# 
# #extract state from file name
# phenolData$State <- sapply(phenolData$File,function(x)strsplit(x,"\\.txt")[[1]][1])
# phenolData$State <- sapply(phenolData$State,function(x)strsplit(x,"_")[[1]][3])
# phenolData <- subset(phenolData, Species==myspecies)

#if no phenolData for a given state, use max and min phenoldays
# dfS <- subset(df, Species==myspecies)
# obsPhenolData <- ddply(dfS,.(State),summarise,
#                        minDay = round(quantile(yday,0.05)),
#                        maxDay = round(quantile(yday,0.95)))
# 
# #expand to list all days between these days
# obsPhenolData <- ddply(obsPhenolData,.(State),function(x){
#                   data.frame(Species=myspecies,
#                              day=as.numeric(x["minDay"]):as.numeric(x["maxDay"]),
#                              fits=NA,
#                              File=NA,
#                              State=x["State"])})
# 
# #remove any states already in the phenolData file
# #then rbind them
# #or just use these estimates??yes
# 
# df <- subset(df,interaction(yday,State) %in% interaction(obsPhenolData$day,obsPhenolData$State))
# ####################################################################################

#subset by average phenology across whole germany
dfS <- subset(df, Species==myspecies)
obsPhenolData <- summarise(dfS,
                           minDay = round(quantile(yday,0.05)),
                           maxDay = round(quantile(yday,0.95)))
df <- subset(df, yday > obsPhenolData$minDay & yday < obsPhenolData$maxDay)

#####################################################################################

#reduce size of the dataset

#reduce data to 5%%
#df <- df[sample(1:nrow(df),round(0.05*nrow(df))),]

# #any oversampled plots???
# out <- ddply(df,.(MTB_Q,Year),summarise,nuDates = length(unique(Date)))
# #out <- arrange(out,desc(nuDates))
# summary(out$nuDates)
# 
# #subset to at most 100 dates per year
# nrow(df)
# df <- ddply(df, .(Year,MTB_Q),function(x){
#   mydates <- ifelse(length(unique(x$Date))>50,
#                     sample(unique(x$Date),50),unique(x$Date))
#   subset(x, Date %in% mydates)
# })
# nrow(df)

######################################################################################

#remove sites visited once
siteSummary <- ddply(df,.(MTB_Q),summarise,nuYears=length(unique(Year)))
df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])

#####################################################################################

#define a visit
#df$visit <- paste(df$MTB_Q,df$Date,sep="_")
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
               nuRecords=length(Species))
               #Richness2=mean(Richness),
               #RpS = length(Species)/length(unique(Species)),
               #expertise = sum(Expertise),
               #samplingSites = length(unique(interaction(lat,lon))))
  
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

#######################################################################################

#add on some indices

#add on box info to listlength
listlengthDF$State <- mtbqsDF$Counties[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- year(listlengthDF$Date)
listlengthDF$yday <- yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))
listlengthDF$stateIndex <- as.numeric(factor(listlengthDF$State))
listlengthDF$siteIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                   listlengthDF$MTB_Q)))
listlengthDF$mtbIndex <- as.numeric(factor(paste0(listlengthDF$stateIndex,
                                                  listlengthDF$MTB)))

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
#subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
#subset(listlengthDF,is.na(Naturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#######################################################################################
#get summary site info data

siteInfo <- unique(listlengthDF[,c("stateIndex","mtbIndex","siteIndex","boxIndex",
                                   "MTB_Q","nnIndex","cnIndex")])
head(siteInfo)

#######################################################################################

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)
siteInfo <- arrange(siteInfo,siteIndex)

raumInfo <- unique(siteInfo[,c("nnIndex","cnIndex")])

########################################################################################

# #get matrix of site versus state
# 
siteInfo$dummy <- 1
siteStates <- acast(siteInfo,siteIndex~stateIndex,value.var="dummy")
siteStates[is.na(siteStates)] <- 0

#get number of sites for each state
statesSiteNu <- as.numeric(colSums(siteStates))
# 
#get number of sites per
nsite_cr <- table(siteInfo$cnIndex)

# #########################################################################################
# #get matrix of site versus raum
# 
siteInfo$dummy <- 1
siteRaums <- acast(siteInfo,siteIndex~cnIndex,value.var="dummy")
siteRaums[is.na(siteRaums)] <- 0

#get number of sites for each state
raumSiteNu <- as.numeric(colSums(siteRaums))

########################################################################################

#fit nation-wide model with random slopes to each box

########################################################################################

#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nstate = length(unique(siteInfo$stateIndex)),
                  nraum = length(unique(siteInfo$nnIndex)),
                  ncraum = length(unique(siteInfo$cnIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  state = listlengthDF$stateIndex,
                  raum = listlengthDF$nnIndex,
                  craum = listlengthDF$cnIndex,
                  year = listlengthDF$yearIndex,
                  stateS = siteInfo$stateIndex,
                  craumS = siteInfo$cnIndex,
                  craumR = raumInfo$cnIndex,
                  nsite_cr = as.numeric(nsite_cr),
                  raumS = siteInfo$nnIndex,
                  yday = listlengthDF$yday - median(listlengthDF$yday),
                  yday2 = listlengthDF$yday^2 - median(listlengthDF$yday^2),
                  nuSpecies = log(listlengthDF$nuSpecies) - median(log(listlengthDF$nuSpecies)),
                  singleList = listlengthDF$singleList,
                  shortList = listlengthDF$shortList,
                  nuRecs = log(listlengthDF$nuRecords) - median(log(listlengthDF$nuRecords)),
                  #nuSS = log(listlengthDF$samplingSites) - median(log(listlengthDF$samplingSites)),# up to 3
                  #expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                  #RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                  y = as.numeric(occMatrix[,myspecies]))
listlengthDF$Species <- bugs.data$y

all(row.names(occMatrix)==listlengthDF$visit)

#the below are used the linear regression model in the model file -see below
bugs.data$sumX <- sum(1:bugs.data$nyear)
bugs.data$sumX2 <- sum((1:bugs.data$nyear)^2)

#add on number of annual visits each year
annualVisits <- ddply(listlengthDF,.(yearIndex),summarise,nuVisits=length(unique(visit)))
bugs.data$annualVisits <- annualVisits$nuVisits

#########################################################################

#indices for the number of detections

#for each i, sum into t 
StrIdx <- array(data=0, dim = c(bugs.data$nvisit,bugs.data$nyear))
for(i in 1:bugs.data$nvisit){
  StrIdx[i,bugs.data$year[i]] <- 1
}
bugs.data$StrIdx <- StrIdx

# siteyear <- unique(data.frame(site = bugs.data$site, 
#                               year = bugs.data$year))
# 
# bugs.data$siteyear <- siteyear
# bugs.data$nsiteyear <- nrow(siteyear)
# 
# #for each i, sum into t 
# StrIdx <- array(data=0, dim = c(bugs.data$nvisit,bugs.data$nsiteyear))
# for(i in 1:100){
#   StrIdx[i,which(bugs.data$site[i] == siteyear$site &  bugs.data$year[i] == siteyear$year)] <- 1
# }
# bugs.data$StrIdx <- StrIdx

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

###########################################################################################

#choose model file
#modelfile="/data/idiv_ess/Odonata/BUGS_sparta_regional_nation_naturraum.txt"
#modelfile="/data/idiv_ess/Odonata/BUGS_sparta_nation_naturraum_phenologyChange.txt"

modelfile="/data/idiv_ess/Odonata/BUGS_sparta_nation_naturraum.txt"

effort = "shortList"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("psi.fs","regres.psi","mean.p","mup","annual.p")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=10,
            n.chains=n.cores, n.burnin=round(niterations/2),
            n.iter=niterations,parallel=T)

Sys.time()

#save as output file - for regional/dynamic model
saveRDS(out,file=paste0("out_sparta_nation_naturraum_",stage,"_",myspecies,".rds"))

#saveRDS(out,file=paste0("out_sparta_nation_naturraum_phenologyChange_",stage,"_",myspecies,".rds"))

#saveRDS(out,file=paste0("out_sparta_nation_naturraum_ppc_",stage,"_",myspecies,".rds"))

########################################################################################





