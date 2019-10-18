########################################################################################

#get libraries
#suppressMessages()
suppressMessages(library(rjags))
suppressMessages(library(R2WinBUGS))
suppressMessages(library(jagsUI))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))

#get random number seed info (contine seed if one has been made)
seed = as.integer(Sys.getenv("SEED", as.integer(Sys.time())))
#seed

#load the relational table of task ids and species
speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_adult.txt"),as.is=T)
#get task id
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 
#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 

#intermediate, unconverged model file name
intermediate_file = paste0("intermediate_dynamic_nation_state_",myspecies,".rds")
intermediate_file

#number of MCMC samples
niterations = 100

#when are we?
Sys.time()

#where are we?
getwd()

#get job chain info
Job_Chain = as.integer(Sys.getenv("Job_Chain","1"))
Job_Chain

#first check whether there is an unconverged model file
out <- if(file.exists(intermediate_file) & Job_Chain>1){
  
  print("i am continuting from the previous run")
  
  intfile <- readRDS(intermediate_file)
  
  update(intfile,n.iter=niterations,seed = seed)
  
} else if(file.exists(intermediate_file) & Job_Chain==1){
  
  #if not, start data processing to fit first model
  print("i am starting from stratch")
  
#print("my seed is", seed)
write.table(seed,file="seed.txt",sep="\t",col.names=FALSE,row.names=FALSE)
  
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

#####################################################################################

#remove MTBQs not in the shapefile = "51561" "50561" "49563"
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)
adultData <- subset(adultData, !MTB_Q %in% c("51561","50561","49563","55484","63012"))

###################################################################################
#filter to 1980 onwards

df <- subset(adultData, Year>=1980  & Year<2017)

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

#remove any states already in the phenolData file
#then rbind them
#or just use these estimates??yes

df <- subset(df,interaction(yday,State) %in% interaction(obsPhenolData$day,obsPhenolData$State))

#####################################################################################

#reduce size of the dataset

#reduce data to 5%%
df <- df[sample(1:nrow(df),round(0.05*nrow(df))),]

#any oversampled plots???
#out <- ddply(df,.(MTB_Q,Year),summarise,nuDates=length(unique(Date)))
#out <- arrange(out,desc(nuDates))
#summary(out$nuDates)

#subset to at most 50 dates per year
nrow(df)
#df <- ddply(df, .(Year,MTB_Q),function(x){
#  mydates <- ifelse(length(unique(x$Date))>20,
#                    sample(unique(x$Date),20),unique(x$Date))
#  subset(x, Date %in% mydates)
#})
nrow(df)

######################################################################################

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
#load("mtbqsDF.RData")
load("/data/idiv_ess/Odonata/mtbqsDF.RData")
all(df$MTB_Q%in% mtbqsDF$MTB_Q)
#unique(df$MTB_Q[!df$MTB_Q %in% mtbqsDF$MTB_Q])

#add on box info to listlength
listlengthDF$km50 <- mtbqsDF$km50[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
summary(listlengthDF$km50)

#######################################################################################

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

siteInfo <- unique(listlengthDF[,c("stateIndex","mtbIndex","siteIndex","boxIndex")])
head(siteInfo)

########################################################################################

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)
siteInfo <- arrange(siteInfo,siteIndex)

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

########################################################################################

#fit nation-wide model with random slopes to each box

########################################################################################

#organise data for BUGS model
bugs.data <- list(nsite = length(unique(listlengthDF$siteIndex)),
                  nyear = length(unique(listlengthDF$yearIndex)),
                  nstate = length(unique(siteInfo$stateIndex)),
                  nbox = length(unique(siteInfo$boxIndex)),
                  nmtb = length(unique(siteInfo$mtbIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  mtb = listlengthDF$mtbIndex,
                  state = listlengthDF$stateIndex,
                  box = listlengthDF$boxIndex,
                  year = listlengthDF$yearIndex,
                  stateS = siteInfo$stateIndex,
                  boxS = siteInfo$boxIndex,
                  mtbS = siteInfo$mtbIndex,
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

#another check
all(row.names(occMatrix)==listlengthDF$visit)

#######################################################################################

#specify initial values
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
zst [is.infinite(zst)] <- 0
inits <- function(){list(z = zst)}

########################################################################################

#define model params
nt <- 1

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
#n.cores = 3
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

###########################################################################################

#modelfile="/data/idiv_ess/Odonata/BUGS_sparta_nation_state.txt"
#modelfile="R/BUGS_dynamic_nation_state.txt"
modelfile="/data/idiv_ess/Odonata/BUGS_dynamic_nation_state.txt"

effort = "nuSpecies"
bugs.data$Effort <- bugs.data[[effort]]

#view data
bugs.data$nstate
table(bugs.data$stateS)
unique(listlengthDF$State)
unique(listlengthDF$stateIndex)

#specify parameters to monitor
#params <- c("int","state.a.effect","state.t.effect","effort.p","single.p","psi.fs")
#params <- c("mean.growth","mean.p","state.persist","state.colonize")
params <- c("mean.growth","mean.p","state.mean.growth",
            "psi.fs","psi.state")


#run model
jags(bugs.data, inits=inits, params, modelfile, n.thin=2,
            n.chains=n.cores, n.burnin=niterations/4,
            n.iter=niterations,parallel=T,
            seed = seed)
#up to here, we are doing first iterations


} else if(file.exists(intermediate_file) & Job_Chain==1){
  quit(status = 1)
} else if (!file.exists(intermediate_file) & Job_Chain>1){
  quit(status = 1)
}

########################################################################################


#20 mins with state random effect - 50 iterations -all data
#14 mins with state random effect - 50 iterations -limit to 20
#45 mins - 500 iteractions - limit to 20 (we want 400 times this) - 300 hours
#40 mins - 500 iterations - visit as date/mtbq - #do 100000

#save as output file
#saveRDS(data.frame(out$summary),file=paste0("outSummary_dynamic_nation_state_",stage,"_",myspecies,".rds"))

########################################################################################

indicator <- data.frame(out$summary)
print(head(indicator))

Sys.time()

if (mean(indicator$Rhat<1.1,na.rm=T)>0.9){
  saveRDS(data.frame(out$summary),file=paste0("outSummary_dynamic_nation_state_",myspecies,".rds"))
  quit(status=0)#successful finish R since model has converged
}else{
  saveRDS(out,file=intermediate_file)
  quit(status=42)#not converged  - abnormal exit
}


