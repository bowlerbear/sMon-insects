########################################################################################

#get all libraries we need
#suppressMessages()
suppressMessages(library(rjags))
suppressMessages(library(R2WinBUGS))
suppressMessages(library(jagsUI))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))

#set stage
stage="adult"

#load in regional datasets
myfiles <- list.files("derived-data")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

#combine these files
adultData <- ldply(adultFiles,function(x){
  out<-readRDS(paste("derived-data",x,sep="/"))
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
all(df$MTB_Q%in% mtbqsDF$MTB_Q)
df$State <- mtbqsDF$Counties[match(df$MTB_Q,mtbqsDF$MTB_Q)]

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

saveRDS(df,file="df_sparta.rds")

######################################################################################

#apply sparta
library(devtools)

# Some users have reported issues with devtools not correctly installing
# dependencies. Run the following lines to avoid these issues
list.of.packages <- c("minqa", "lme4", "gtools", "gtable", "scales",
                      "assertthat", "magrittr", "tibble", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Now install sparta
install_github('BiologicalRecordsCentre/sparta')

# Load sparta
library(sparta)

?occDetModel
results <- occDetModel(taxa = df$Species,
                       site = df$MTB_Q,
                       survey = df$Date,
                       species_list = 'Aeshna cyanea',
                       write_results = TRUE,
                       n_iterations = 10000,
                       burnin = 100,
                       thinning = 2,
                       seed = 125, 
                       modeltype = c("ranwalk", "halfcauchy",
                                     "catlistlength","jul_date"))

