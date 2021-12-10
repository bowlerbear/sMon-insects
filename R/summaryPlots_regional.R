library(rgdal)
library(ggplot2)
library(tidyverse)
library(raster)
library(maptools)
library(sp)
library(lubridate)
library(tmap)
library(sf)

#helper functions
source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

### Species list ##############################################################################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

### German maps #############################################################################################

#MTBQ - utm
mtbqMap <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                   layer="MTBQ_25833")
mtbqMap <- addMTBQ(mtbqMap)
mtbqMap_sf <- st_as_sf(mtbqMap)

#convert to raster - grid in lat/lon
mtbqMapR <- spTransform(mtbqMap,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
mtbqDF <- data.frame(mtbqMapR@data,x=coordinates(mtbqMapR)[,1],y=coordinates(mtbqMapR)[,2])
pixels <- SpatialPixelsDataFrame(points=mtbqDF[,c('x','y')], 
                                 data=mtbqDF[,c('Quadrant','MTB_Q')],
                                 tolerance = 0.916421)
raster <- raster(pixels[,'MTB_Q'])
plot(raster)

#in utm
load("mtbqsDF.RData")

#MTB
mtbMap <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",layer="MTB_25832")
mtbMap <- spTransform(mtbMap,CRS(proj4string(mtbqMap)))#utm
mtbMap_sf <- st_as_sf(mtbMap)

#get map of Germany
germanyMap <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
germanyMap <- spTransform(germanyMap,crs(mtbqMap))

###Adult records##############################################################################################

myfiles <- list.files("derived-data")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

#when available take updated data files
adultFiles <- adultFiles[-c(3,5,7)]

#combine these files
adultData <- plyr::ldply(adultFiles,function(x){
  out<-readRDS(paste("derived-data",x,sep="/"))
  out$File <- x
  return(out)
})


#extract state from file name
adultData$State <- sapply(adultData$File,function(x)strsplit(x,"\\.rds")[[1]][1])
adultData$State <- sapply(adultData$State,function(x)strsplit(x,"_")[[1]][3])
nrow(adultData)#1147558

#where do we have/ in the MTBQs
#only in the Bavaria dataset... fix at some point
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)

#add iNaturalist data to fill gaps
gbifdata <- readRDS("derived-data/datafile_iNaturalist.rds")
#combine the two
adultData <- rbind(adultData,gbifdata)
nrow(adultData)#1198708

#change state labels
adultData$State <- as.factor(adultData$State)
levels(adultData$State) <- c("Bavaria","Brandenberg","Baden-Wuttemberg","Hesse",
                             "Mecklenburg-Vorpommern",
                             "North Rhine-Westphalia","Lower Saxony",
                             "Rheinland-Pfalz","Saarland","Saxony-Anhalt",
                             "Saxony","Schleswig Holstein","Thuringia")

#format data
adultData$Date <- as.Date(adultData$Date)
adultData$Year <- year(adultData$Date)
adultData$yday <- yday(adultData$Date)
adultData$week <- week(adultData$Date)
adultData$Month <- month(adultData$Date)
adultData <- subset(adultData,!is.na(Date))
adultData <- subset(adultData,yday!=1)

#check  names
species<- read.delim("specieslist_odonata.txt")
adultSpecies <- sort(unique(adultData$Species))
adultSpecies[!adultSpecies%in%species$Species]#should be empty

#get MTB
adultData$MTB <- sapply(as.character(adultData$MTB_Q),function(x){
  substr(x,1,(nchar(x)-1))
})

#state to MTBQ mapping
stateDir <- unique(adultData[,c("MTB","MTB_Q","State")])


#add on ecoregion
df$Naturraum <- mtbqsDF$MTB_CoarseNatur[match(df$MTB,mtbqsDF$Value)]
df$Naturraum <- sapply(df$Naturraum,function(x)
  strsplit(x,"/")[[1]][2])

### define df ########################################################################

df <- subset(adultData, Year>=1990  & Year<=2016)

### remove sites visited once ########################################################

siteSummary <- plyr::ddply(df,"MTB_Q",summarise,nuYears=length(unique(Year)))
df <- subset(df, MTB_Q %in% siteSummary$MTB_Q[siteSummary$nuYears>1])
nrow(df)

#### effort calc #################################################################

#define a visit
df$visit <- paste(df$MTB_Q,df$Date,df$Beobachter,sep="_")

#summarize number of visits per year
annualVisits <- df %>%
                group_by(Year,MTB_Q) %>%
                summarise(nuVisits = length(unique(visit)), 
                          nuDates = length(unique(Date)),
                          nuRecs = length(Date)) %>%
                ungroup() %>%
                complete(Year,MTB_Q) %>%
                replace(is.na(.),0)


annualVisits$State <- stateDir$State[match(annualVisits$MTB_Q,stateDir$MTB_Q)]

#### state effort ts ##########################################################

annualVisits_byState <- df %>%
  group_by(Year,State) %>%
  summarise(nuVisits = length(unique(visit)), 
            nuDates = length(unique(Date)),
            nuRecs = length(Date)) %>%
  ungroup() %>%
  complete(Year,State) %>%
  replace(is.na(.),0)

ggplot(annualVisits_byState)+
  geom_line(aes(x = Year, y = nuDates))+
  facet_wrap(~State)

ggplot(annualVisits_byState)+
  geom_line(aes(x = Year, y = (nuDates+1))) +
  facet_wrap(~State) +
  scale_y_log10()
  
#### ecoregions effort ts #####################################################

annualVisits_byEcoregion <- df %>%
  group_by(Year,Naturraum) %>%
  summarise(nuVisits = length(unique(visit)), 
            nuDates = length(unique(Date)),
            nuRecs = length(Date)) %>%
  ungroup() %>%
  complete(Year,Naturraum) %>%
  replace(is.na(.),0)

ggplot(annualVisits_byEcoregion)+
  geom_line(aes(x = Year, y = nuVisits))+
  facet_wrap(~Naturraum)

ggplot(annualVisits_byEcoregion)+
  geom_line(aes(x = Year, y = (nuVisits+1))) +
  facet_wrap(~Naturraum) +
  scale_y_log10()


#### effort trend map #########################################################

#get trends in effort
annualTrends <- annualVisits %>%
  group_by(MTB_Q) %>%
  do(model = broom::tidy(glm(nuDates ~ Year, data = ., family=poisson))) %>% 
  unnest(model) %>%
  filter(term=="Year") %>%
  rename(Trend="estimate")

#limit outliers
summary(annualTrends$Trend)
annualTrends$Trend[annualTrends$Trend>quantile(annualTrends$Trend,0.95)] <- as.numeric(quantile(annualTrends$Trend,0.95))
annualTrends$Trend[annualTrends$Trend<quantile(annualTrends$Trend,0.05)] <- as.numeric(quantile(annualTrends$Trend,0.05))

#plot the trend
mtbqMap_sf$Trend <- annualTrends$Trend[match(mtbqMap_sf$MTB_Q,annualTrends$MTB_Q)]
tm_shape(mtbqMap_sf) + tm_fill("Trend")

#### decadal maps #############################################################

annualVisits$Decade <- plyr::round_any(annualVisits$Year,5)

decadeVisits <- annualVisits %>%
                group_by(Decade,MTB_Q,State) %>%
                summarise(nuDates = sum(nuDates),
                          nuYears = length(unique(Year)))
  
mtbqMap_sf_decade <- inner_join(mtbqMap_sf,decadeVisits,by="MTB_Q")
mtbqMap_sf_decade$nuDates[mtbqMap_sf_decade$nuDates==0] <- NA
tm_shape(mtbqMap_sf_decade) +
  tm_fill("nuDates", style="quantile") +
  tm_facets(by="Decade")

#### MTB effort trends #######################################################

annualVisits <- df %>%
  group_by(Year,MTB) %>%
  summarise(nuVisits = length(unique(visit)), 
            nuDates = length(unique(Date)),
            nuRecs = length(Date)) %>%
  ungroup() %>%
  complete(Year,MTB) %>%
  replace(is.na(.),0)
  
annualTrends <- annualVisits %>% 
  group_by(MTB) %>%
  do(model = broom::tidy(glm(nuDates ~ Year, data = ., family=poisson))) %>% 
  unnest(model) %>%
  filter(term=="Year") %>%
  rename(Trend="estimate")

#limit outliers
summary(annualTrends$Trend)
annualTrends$Trend[annualTrends$Trend>quantile(annualTrends$Trend,0.95)] <- as.numeric(quantile(annualTrends$Trend,0.95))
annualTrends$Trend[annualTrends$Trend<quantile(annualTrends$Trend,0.05)] <- as.numeric(quantile(annualTrends$Trend,0.05))

#plot
mtbMap_sf$Trend <- annualTrends$Trend[match(mtbMap_sf$Value,annualTrends$MTB)]
tm_shape(mtbMap_sf)+
  tm_fill("Trend")

#### Bavaria #################################################################

#are there really fewer records in Bavaria in the last years?
annualTrends$State <- stateDir$State[match(annualTrends$MTB,stateDir$MTB)]
#1990 was a big recording year!

#get MTBs with negative trends
annualTrends_Bavaria <- subset(annualTrends,State=="Bavaria")
hist(annualTrends_Bavaria$Trend)
mean(annualTrends_Bavaria$Trend<0)#65%???

decliningMTBs <- annualTrends_Bavaria$MTB[annualTrends_Bavaria$Trend<(-0.1)]
print(subset(annualVisits,MTB == decliningMTBs[1]),n=50)
#5527 in 1990/2009/2016
print(subset(annualVisits,MTB == decliningMTBs[2]),n=50)
#5527 in 1990/2009/2016
print(subset(annualVisits,MTB == decliningMTBs[4]),n=50)
#5527 in 1990/2009/2016

#yes, looks like there has been fewer

#### end ######################################################################