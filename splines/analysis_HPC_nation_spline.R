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
niterations = 50000

Sys.time()

#load in regional datasets
myfiles <- list.files("derived-data")
#myfiles <- list.files("/data/idiv_ess/Odonata")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

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
gbifdata <- readRDS("derived-data/datafile_GBIF.rds")
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

#####################################################################################

#remove MTBQs not in the shapefile = "51561" "50561" "49563"
#remove MTBQS without naurraum 81443, 44553, 58401
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)
adultData <- subset(adultData, !MTB_Q %in% c("51561","50561","49563","55484","63012",
                                             "81443","44553","58401"))

###################################################################################
#filter to 1980 onwards

df <- subset(adultData, Year>=1985  & Year<2017)


###remove baden wuttemberg #########################################################

df <- subset(df,State!="Baden-Württemberg")
df <- subset(df,!is.na(State))

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
myspecies = "Crocothemis erythraea"
#stage="adult"

#####################################################################################

#remove sites visited once
siteSummary <- ddply(df,.(MTB_Q),summarise,nuYears=length(unique(Year)))
df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])

######################################################################################

#get nationwide boxes
load("mtbqsDF.RData")
#load("/data/idiv_ess/Odonata/mtbqsDF.RData")
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

# dfS <- subset(df, Species==myspecies)
# 
# #across all species
# obsPhenolData <- summarise(df,
#                            minDay = round(quantile(yday,0.05)),
#                            maxDay = round(quantile(yday,0.95)))
# df <- subset(df, yday > obsPhenolData$minDay & yday < obsPhenolData$maxDay)

#####################################################################################

#reduce size of the dataset

#reduce data to 5%%
#df <- df[sample(1:nrow(df),round(0.05*nrow(df))),]

#any oversampled plots???
# out <- ddply(df,.(MTB_Q,Year),summarise,nuDates = length(unique(Date)))
# #out <- arrange(out,desc(nuDates))
# summary(out$nuDates)
# 
# #subset to at most 50 dates per year
# nrow(df)
# df <- ddply(df, .(Year,MTB_Q),function(x){
#   mydates <- ifelse(length(unique(x$Date))>30,
#                     sample(unique(x$Date),30),unique(x$Date))
#   subset(x, Date %in% mydates)
# })
# nrow(df)

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

#######################################################################################

#add on some indices

#add on box info to listlength
listlengthDF$km50 <- mtbqsDF$km50[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
summary(listlengthDF$km50)

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

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB_Q,mtbqsDF$MTB_Q)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
subset(listlengthDF,is.na(Naturraum))

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

#get matrix of site versus state

siteInfo$dummy <- 1
siteRaums <- acast(siteInfo,siteIndex~cnIndex,value.var="dummy")
siteRaums[is.na(siteRaums)] <- 0

#get number of sites for each state
raumSiteNu <- as.numeric(colSums(siteRaums))

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
                  nraum = length(unique(siteInfo$nnIndex)),
                  ncraum = length(unique(siteInfo$cnIndex)),
                  nbox = length(unique(siteInfo$boxIndex)),
                  nvisit = nrow(listlengthDF),
                  site = listlengthDF$siteIndex,
                  state = listlengthDF$stateIndex,
                  raum = listlengthDF$nnIndex,
                  craum = listlengthDF$cnIndex,
                  box = listlengthDF$boxIndex,
                  year = listlengthDF$yearIndex,
                  stateS = siteInfo$stateIndex,
                  craumS = siteInfo$cnIndex,
                  raumS = siteInfo$nnIndex,
                  boxS = siteInfo$boxIndex,
                  yday = listlengthDF$yday - median(listlengthDF$yday),
                  nuSpecies = log(listlengthDF$nuSpecies) - log(median(listlengthDF$nuSpecies)),
                  singleList = listlengthDF$singleList,
                  shortList = listlengthDF$shortList,
                  nuRecs = log(listlengthDF$nuRecords) - log(median(listlengthDF$nuRecords)),
                  nuSS = log(listlengthDF$samplingSites) - log(median(listlengthDF$samplingSites)),# up to 3
                  expertise = log(listlengthDF$expertise)-median(log(listlengthDF$expertise)),
                  RpS = log(listlengthDF$RpS) - median(log(listlengthDF$RpS)),
                  y = as.numeric(occMatrix[,myspecies]),
                  siteRaums = siteRaums,
                  nsiteRaum = raumSiteNu)
listlengthDF$Species <- bugs.data$y

all(row.names(occMatrix)==listlengthDF$visit)

##set prior to zero if species never recorded in that state??
#temp <- ddply(listlengthDF,.(cnIndex),summarise,species=sum(Species))
#bugs.data$priorS <- ifelse(temp$species>0,0.99999,0.05)

#######################################################################################

# #set up the spline info
siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]

summary(siteInfo$x)
summary(siteInfo$y)

#centre and scale them?
siteInfo$x <- (siteInfo$x - median(siteInfo$x))/10000
siteInfo$y <- (siteInfo$y - median(siteInfo$y))/10000
  
#save(siteInfo,file="siteInfo.RData") 

# #get info on whether species was seen at each site
speciesSite <- ddply(listlengthDF,.(siteIndex),summarise,PA=max(Species,na.rm=T))
siteInfo$obs <- speciesSite$PA[match(siteInfo$siteIndex,speciesSite$siteIndex)]

# #fit normal gam
# library(mgcv)
# gam1 <- gam(obs ~ 1 + s(x,y), data=siteInfo, family=binomial)
# siteInfo$fits <- gam1$fitted.values
# library(ggplot2)
# qplot(x,y,data=siteInfo,colour=fits)+
#   scale_colour_gradient(low="blue",high="red")
# 
# ##########################################################################################
# 
# #set up space-time-spline info
# siteInfo <- unique(listlengthDF[,c("MTB_Q","stateIndex","siteIndex","boxIndex","yearIndex")])
# siteInfo$x <- mtbqsDF$x[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
# siteInfo$y <- mtbqsDF$y[match(siteInfo$MTB_Q,mtbqsDF$MTB_Q)]
# speciesSite <- ddply(listlengthDF,.(siteIndex,yearIndex),summarise,PA=max(Species,na.rm=T))
# siteInfo$obs <- speciesSite$PA[match(interaction(siteInfo$siteIndex,siteInfo$yearIndex),
#                                      interaction(speciesSite$siteIndex,speciesSite$yearIndex))]
# siteInfo <- subset(siteInfo,!is.na(obs))
# siteInfo <- subset(siteInfo,!is.na(x))
# 
# #fit normal gam
# library(mgcv)
# gam1 <- gam(obs ~ s(yearIndex) + s(x,y,yearIndex), data=siteInfo, family=binomial)
# siteInfo$fits <- gam1$fitted.values
# library(ggplot2)
# qplot(x,y,data=siteInfo,colour=fits)+
#   scale_colour_gradient(low="blue",high="red")+
#   facet_wrap(~yearIndex)

#######################################################################################

#get BUGS code using jagam
library(mgcv)
jags.ready <- jagam(obs ~ 1 + s(x,y,k=10), data=siteInfo, family=binomial,file="jagam.txt")


#extract bits for the model
bugs.data$X = jags.ready$jags.data$X#[,-1]
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero#[-1]
bugs.data$nspline = length(bugs.data$zero)
bugs.data$nspline1 = length(bugs.data$zero)+1
bugs.data$nspline2 = (length(bugs.data$zero))*2


#expand to include NA values ##########################################################

siteInfo_NAs <- mtbqsDF

#remove without county or baden-wuttenberg
siteInfo_NAs <- subset(siteInfo_NAs,!is.na(Counties))
siteInfo_NAs <- subset(siteInfo_NAs,Counties!="Baden-Württemberg")

siteInfo_NAs$obs <- siteInfo$obs[match(siteInfo_NAs$MTB_Q,siteInfo$MTB_Q)]
table(siteInfo_NAs$obs)
sum(is.na(siteInfo_NAs))
nrow(siteInfo_NAs)

#randomly insect 0 or 1 into the missing values
siteInfo_NAs$obs[is.na(siteInfo_NAs$obs)] <- sample(c(0,1),size=sum(is.na(siteInfo_NAs$obs)),replace=T)
saveRDS(siteInfo_NAs,file="splines/siteInfo_NAs.rds")

library(mgcv)
jags.ready <- jagam(obs ~ 1 + s(x,y,k=10), 
                    data=siteInfo_NAs, family=binomial,file="jagam.txt")
saveRDS(jags.ready,file="splines/jags.ready_NAs.rds")

########################################################################################

#specify initial values
library(reshape2)
zst <- acast(listlengthDF, siteIndex~yearIndex, value.var="Species",fun=max)
zst [is.infinite(zst)] <- 0
#inits <- function(){list(z = zst)}

inits <- function(){list(z = zst,
                         #state.a = runif(bugs.data$nstate,0.0001,0.01),
                         #lphi = runif(bugs.data$nstate,0,0.01),
                         #lgam = runif(bugs.data$nstate,0,0.01),
                         effort.p = runif(1,0,0.01),
                         mu.phenol = runif(1,-0.01,0.01),
                         mu.phenol2 = runif(1,-0.01,0.01))}

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
#n.cores = 3

###########################################################################################

#modelfile="R/BUGS_dynamic_nation_spline.txt"
modelfile="/data/idiv_ess/Odonata/BUGS_dynamic_nation_spline_year.txt"

effort = "nuSpecies"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
#params <- c("persistMean","colonizeMean")
params <- c("z")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=20,
            n.chains=n.cores, n.burnin=niterations/4,n.iter=niterations,parallel=T)
Sys.time()

#save as output file
saveRDS(out$summary,file=paste0("outSummary_dynamicspline_",stage,"_", myspecies,".rds"))

########################################################################################

#get nationwide boxes
load("mtbqsDF.RData")

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

#model that includes NAs
out <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/splines/outSummary_static_spline_NAs.rds")
siteInfo_NAs <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/splines/siteInfo_NAs.rds")
siteInfo_NAs$preds <- out[1:10882,1]

#add on x and y of the mtbq
siteInfo_NAs$x <- mtbqsDF$x[match(siteInfo_NAs$MTB_Q,mtbqsDF$MTB_Q)]
siteInfo_NAs$y <- mtbqsDF$y[match(siteInfo_NAs$MTB_Q,mtbqsDF$MTB_Q)]

library(ggplot2)
ggplot(siteInfo_NAs)+
  geom_point(aes(x=x,y=y,color=preds))+
  scale_color_viridis_c()
#looks nice!!!

########################################################################################