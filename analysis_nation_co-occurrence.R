########################################################################################

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
#get task id
#task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 
#get species for this task
#myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 
myspecies="Aeshna cyanea"
myspecies="Aeshna caerulea"#id2

#set stage
stage="adult"

#number of MCMC samples
niterations = 100000

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
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)
adultData <- subset(adultData, !MTB_Q %in% c("51561","50561","49563","55484","63012"))

###################################################################################
#filter to 1980 onwards

df <- subset(adultData, Year>=1980  & Year<2017)

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
#table(df$State2,df$State)

#####################################################################################

#remove any oversampled plots???
out <- ddply(df,.(MTB_Q,Year),summarise,nuDates = length(unique(Date)))
#out <- arrange(out,desc(nuDates))
summary(out$nuDates)

#subset to at most 20 dates per year
nrow(df)
df <- ddply(df, .(Year,MTB_Q),function(x){
  mydates <- ifelse(length(unique(x$Date))>50,
                    sample(unique(x$Date),50),unique(x$Date))
  subset(x, Date %in% mydates)
})
nrow(df)

######################################################################################

#define a visit
#df$visit <- paste(df$MTB_Q,df$Date,sep="_")
df$visit <- paste(df$MTB_Q,df$Date,df$Beobachter,sep="_")

#get occurence matrix  - detection and non-detection
getOccurrenceMatrix<-function(df){
  require(reshape2)
  out <- acast(df,Species~visit,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
  out[out>0] <- 1
  return(out)
}
occMatrix <- getOccurrenceMatrix(df)

###############################################################################

#Co-occurrence
#probably of genus appearing on a visit
occMatrix[1:5,1:5]

library(cooccur)
cooccur.finches <- cooccur(mat = occMatrix, type = "spp_site",thresh = TRUE, spp_names = TRUE)
summary(cooccur.finches)
temp <- prob.table(cooccur.finches)
plot(cooccur.finches)

