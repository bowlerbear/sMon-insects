
#get population data
library(plyr)
############################################################################################

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
library(lubridate)
adultData$Date<-as.Date(adultData$Date)
adultData$Year <- year(adultData$Date)
adultData$yday <- yday(adultData$Date)
adultData$week <- week(adultData$Date)
adultData$Month <- month(adultData$Date)
adultData$Day <- day(adultData$Date)
adultData <- subset(adultData,!is.na(Date))
adultData <- subset(adultData,yday!=1)

df <- subset(adultData, Year>=1975)

##########################################################################################

#get german county boundaries
germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")

#get MTBQs
library(rgdal)
mtbqs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                 layer="MTBQ_25833")

##################################################################################

#read in data from ralf schafer
tdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/WaterQuality/Data_monitoring/Data_monitoring"
samples <- read.csv(paste(tdir,"samples.csv",sep="/"))
sites <- read.csv(paste(tdir,"sites.csv",sep="/"))
head(samples)
head(sites)

#plot the sites
library(sp)
coordinates(sites) <- c("easting","northing")
proj4string(sites) <- CRS("+init=epsg:31467")
#sites <- spTransform(sites, CRS(proj4string(germanAdmin)))
#plot(sites,add=T)

#overlay the coordinates with the mtqb map
sites <- spTransform(sites, CRS(proj4string(mtbqs)))
plot(mtbqs)
plot(sites,add=T,col="red")

#get the mtbq of each site
sitesMTB_Q <- over(sites,mtbqs)
sitesMTB_Q$Q <- NA
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="NW")]<-1
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="NO")]<-2
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="SW")]<-3
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="SO")]<-4
sitesMTB_Q$MTB_Q <- paste0(as.character(sitesMTB_Q$Value),
                        as.character(sitesMTB_Q$Q))
head(sitesMTB_Q)

#########################################################################################

#add the MTBQ to the site file
sites <- read.csv(paste(tdir,"sites.csv",sep="/"))
sites$MTB_Q <- sitesMTB_Q$MTB_Q
samples$MTB_Q <- sites$MTB_Q[match(samples$site_id,sites$site_id)]

#merge with samples
library(lubridate)
samples$Date <- as.Date(samples$date)
samples$Year <- year(samples$Date)
samples$yday <- yday(samples$Date)
samples$week <- week(samples$Date)
samples$Month <- month(samples$Date)
samples <- subset(samples,!is.na(date))
samples <- subset(samples,Year>=1975 & Year<2019)

summary(samples$Year)

########################################################################################

#get info on which MTBS have water qualities samples (least 5 years)
waterSummary <- ddply(samples,.(MTB_Q),summarise,
                      nuYears=length(unique(Year)),
                      TimeSpan=max(Year)-min(Year))
waterSummary <- subset(waterSummary, nuYears>=5 & TimeSpan>=10)

#get info on which MTBS
dfSummary <- ddply(df, .(MTB_Q),summarise,
                   nuYears=length(unique(Year)),
                   TimeSpan=max(Year)-min(Year))
dfSummary <- subset(dfSummary, nuYears>=5 & TimeSpan>=10)

#get samples overlapping
samples <- subset(samples, MTB_Q %in% waterSummary$MTB_Q)
samples <- subset(samples, MTB_Q %in% dfSummary$MTB_Q)
nrow(samples)

#nrow(samples)
#[1] 7658027
# length(unique(samples$MTB_Q))
#[1] 1014
write.table(samples,file="samplesforDragonflies.txt",sep="\t",row.names=FALSE)
########################################################################################