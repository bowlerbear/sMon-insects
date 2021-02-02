#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)
source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

###Species list#################################################################################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

###German maps#############################################################################################

library(maptools)
library(sp)

#get map of Germany
germanyMap <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")

#MTBQ
mtbqMap <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                   layer="MTBQ_25833")

#convert to raster
library(raster)
mtbqMapR <- spTransform(mtbqMap,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
mtbqDF <- data.frame(mtbqMapR@data,x=coordinates(mtbqMapR)[,1],y=coordinates(mtbqMapR)[,2])
mtbqDF$Q <- NA
mtbqDF$Q[which(mtbqDF$Quadrant=="NW")]<-1
mtbqDF$Q[which(mtbqDF$Quadrant=="NO")]<-2
mtbqDF$Q[which(mtbqDF$Quadrant=="SW")]<-3
mtbqDF$Q[which(mtbqDF$Quadrant=="SO")]<-4
mtbqDF$MTB_Q <- paste0(mtbqDF$Value,mtbqDF$Q)
pixels <- SpatialPixelsDataFrame(points=mtbqDF[,c('x','y')], 
                                 data=mtbqDF[,c('Q','MTB_Q')],
                                 tolerance = 0.916421)
raster <- raster(pixels[,'MTB_Q'])
plot(raster)

#in utm
mtbqDF <- data.frame(mtbqMap@data,x=coordinates(mtbqMap)[,1],y=coordinates(mtbqMap)[,2])
mtbqDF$Q <- NA
mtbqDF$Q[which(mtbqDF$Quadrant=="NW")]<-1
mtbqDF$Q[which(mtbqDF$Quadrant=="NO")]<-2
mtbqDF$Q[which(mtbqDF$Quadrant=="SW")]<-3
mtbqDF$Q[which(mtbqDF$Quadrant=="SO")]<-4
mtbqDF$MTB_Q <- paste0(mtbqDF$Value,mtbqDF$Q)

#MTB
mtbMap <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                  layer="MTB_25832")

#Get crs
proj4string(mtbqMap)
proj4string(mtbMap)
mtbMap <- spTransform(mtbMap,CRS(proj4string(mtbqMap)))

#Overlay
plot(mtbMap)
plot(mtbqMap,add=T,col="red")
#both in "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"

###Adult records##############################################################################################

myfiles <- list.files("derived-data")

#read in and combine all adult files
adultFiles <- myfiles[grepl("adult_datafile",myfiles)]
adultFiles <- adultFiles[grepl("rds",adultFiles)]

#when available take updated data files
adultFiles <- adultFiles[-c(3,5,7)]

#combine these files
adultData <- ldply(adultFiles,function(x){
  out<-readRDS(paste("derived-data",x,sep="/"))
  out$File <- x
  return(out)
})


#extract state from file name
adultData$State <- sapply(adultData$File,function(x)strsplit(x,"\\.rds")[[1]][1])
adultData$State <- sapply(adultData$State,function(x)strsplit(x,"_")[[1]][3])
nrow(adultData)#1147558

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

library(lubridate)
adultData$Date <- as.Date(adultData$Date)
adultData$Year <- year(adultData$Date)
adultData$yday <- yday(adultData$Date)
adultData$week <- week(adultData$Date)
adultData$Month <- month(adultData$Date)
adultData <- subset(adultData,!is.na(Date))
adultData$MTB_Q <- gsub("/","",adultData$MTB_Q)

#fix names
species<- read.delim("specieslist_odonata.txt")
adultSpecies <- sort(unique(adultData$Species))
adultSpecies[!adultSpecies%in%species$Species]

####Site visitation#########################################################################################

#how often is each site visited
out <- ddply(adultData,.(Year,MTB_Q),summarise,nu = length(unique(yday)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   3.107   3.000 125.000

#how many sites is a species usually present at
out <- ddply(adultData,.(Species),summarise,nu = length(unique(MTB_Q)))
length(unique(adultData$MTB_Q))#8231
summary(out$nu/8231)

#on how many visits is a species usually detected
out <- ddply(adultData,.(MTB_Q,Year),summarise,nuVisits=length(unique(yday)))
out2 <- ddply(adultData,.(Species,MTB_Q,Year),summarise,nuDets=length(unique(yday)))
out <- merge(out,out2,by=c("MTB_Q","Year"))
out$p <- out$nuDets/out$nuVisits
summary(out$p)

out <- ddply(out,.(Species),summarise,meanp=mean(p))
summary(out$meanp)

####Flight period###########################################################################################

#flight period check:

#for each species, in each state:
#each lower and upper 10% week of observations

adultDataS <- subset(adultData,Month%in%4:10)
adultDataS <- subset(adultDataS, Anzahl_min>0 & !is.na(Anzahl_min))

out <- ddply(adultDataS,.(Species,State),summarise,
      lower=quantile(week,0.05),
      upper=quantile(week,0.95),
      med=(lower+upper)/2)


#plot
library(ggplot2)
length(unique(out$Species))#79
ggplot(subset(out,Species%in%unique(out$Species)[1:40]))+
  geom_crossbar(aes(x=Species,y=med,ymin=lower,ymax=upper))+
  coord_flip()+
  facet_wrap(~State,nrow=1)
  
#remove genus
ggplot(subset(out,Species%in%unique(out$Species)[1:39]))+
  geom_crossbar(aes(x=Species,y=med,ymin=lower,ymax=upper))+
  coord_flip()+
  facet_wrap(~State,nrow=1)

####Serial dependence###########################################################################################

#serial dependence of observation?

#plot date of first observation for single-list
#plot date of first observation for all lists

#convert dat to julian date

#take Berd Trockur data
adultData_BT <- subset(adultData, Beobachter=="Bernd Trockur")

#for each species, in each year, get the first julian date
listSummary_BT <- ddply(adultData_BT,.(Year,State,Beobachter,Date),summarise,
                       nuList = length(unique(Species)),
                       nuRecs = length(Species))

adultData_BT <- merge(adultData_BT,listSummary_BT,merge=c("Year","State","Beobachter","Date"))

#single list or not
adultData_BT$ListType <- ifelse(adultData_BT$nuRec==1,"single list","multiple list")

#compare julian date each year for single and multiple lists
ggplot(adultData_BT)+
  geom_boxplot(aes(x=as.factor(Species),y=yday,colour=ListType))+
  coord_flip()

#same with all observers
listSummary <- ddply(adultData,.(Year,State,Beobachter,Date),summarise,
                        nuList = length(unique(Species)),
                        nuRecs = length(Species))
adultData <- merge(adultData,listSummary,merge=c("Year","State","Beobachter","Date"))
adultData$ListType <- ifelse(adultData$nuRec==1,"single list","multiple list")

ggplot(subset(adultData,Year>2010 & Year<2016))+
  geom_boxplot(aes(x=as.factor(Species),y=yday,colour=ListType))+
  coord_flip()+
  facet_wrap(~Year,nrow=1)

#look at total richness at sites and when they are studied
richnessSummary <- ddply(adultData,.(State,MTB_Q),summarise,
                         richness = length(unique(Species)))

adultData <- merge(adultData,richnessSummary,
                            merge=c("State","MTB_Q"))
head(adultData)

#look at average richness of sites sampled in each year
ggplot(subset(adultData,Year>1979 & Year <2016))+
  geom_boxplot(aes(x=as.factor(Year),y=Richness),outlier.shape=NA)+
  facet_wrap(~State)+
  coord_flip()
  
####Record map###########################################################################################

#plot number of records per MTB

adultDataS <- subset(adultData,Year>=1980 & Year <=2016)

#extract MTB from Q
adultDataS$MTB_Q <- gsub("/","",adultDataS$MTB_Q)

#if MTB is 5 characters, take first 4 as MTB and last character as Q
adultDataS$Q<- unlist(sapply(adultDataS$MTB, function(x){
  if(nchar(x)==5){
    substr(x,5,6)}
  else 
    x
}))

adultDataS$MTB<- sapply(adultDataS$MTB, function(x){
                                            if(nchar(x)==5){
                                              substr(x,1,4)
                                            }
                                            else x
})

#look at data without a quadrant
temp <- subset(adultDataS, as.numeric(Q)>4 |is.na(Q))
table(temp$Stat)#just 29 in SH
#these probably are supposed to have a zero at the start
adultDataS$MTB[which(adultDataS$MTB_Q==9162)]<-"916"
adultDataS$Q[which(adultDataS$MTB_Q==9162)]<-"2"
adultDataS$MTB[which(adultDataS$MTB_Q==9163)]<-"916"
adultDataS$Q[which(adultDataS$MTB_Q==9163)]<-"3"
adultDataS$MTB[which(adultDataS$MTB_Q==9164)]<-"916"
adultDataS$Q[which(adultDataS$MTB_Q==9164)]<-"4"


#work out records per quadrant and quadrant
nuRecs <- ddply(adultDataS,.(MTB,Q),
                summarise,nuRecs = length(Species),
                          nuYears = length(unique(Year)))


#plot
allDF <- merge(mtbqDF,nuRecs,by.x=c("Value","Q"),by.y=c("MTB","Q"))
nrow(mtbqDF)
nrow(allDF)
#data from 9992/12024

#plot map of germany
library(ggplot2)
germanyMap <- spTransform(germanyMap,proj4string(mtbMap))
AG <- fortify(germanyMap)

ggplot()+ geom_polygon(data=AG, aes(long, lat, group = group), colour = "grey", fill=NA)+
  geom_point(data=allDF,aes(x=x,y=y,colour=nuRecs),
             shape=15,size=rel(0.45))+
  scale_colour_viridis_c("Number of records",trans="log",
                         breaks=c(1,20,400,8000),
                         labels=c(1,20,400,8000),
                         option="A",
                         direction=-1)+
  xlab("X")+ylab("Y")+
  coord_equal()+
  theme_void()

ggsave(filename="plots/Adult_map_effort.png")
ggsave(filename="plots/Juv_map_effort.png",width=8,height=7)

#exclude sites visited once
gMap <- ggplot()+ geom_polygon(data=AG, aes(long, lat, group = group), 
                       colour = "grey", fill=NA)+
  geom_point(data=subset(allDF,nuYears>1),
             aes(x=x,y=y,colour=nuRecs),
             shape=15,size=rel(1.2))+
  scale_colour_viridis_c("Number of records",trans="log",
                         breaks=c(1,20,400,8000),
                         labels=c(1,20,400,8000),
                         option="A",
                         direction=-1)+
  xlab("X")+ylab("Y")+
  coord_equal()+
  theme_void()+
  theme(legend.position = "top")

gMap
ggsave(filename="plots/Adult_map_effort_sitesubset.png",width=8,height=7)
sum(nuRecs$nuRecs[nuRecs$nuYears>1])

#do we have data for all naturraums?
load("mtbqsDF.RData")
mtbqsDF$Data <- ifelse(mtbqsDF$MTB_Q %in% allDF$MTB_Q,1,0)

#coarse
ddply(mtbqsDF,.(CoarseNatur),summarise,nuRecs=sum(Data))

#fine
temp <- ddply(mtbqsDF,.(Natur),summarise,nuRecs=sum(Data))
subset(temp, nuRecs<10)#only one with zero!!!

####Record time series############################################################################################

#Time-series for whole Germany
timeSummary <- ddply(adultData,.(Year),summarise,
                     nuRecs=length(Species),
                     nuSpecies=length(unique(Species)),
                     nuPlots=length(unique(MTB_Q)))

g1 <- ggplot(subset(timeSummary,Year>1950 & Year <2017))+
  geom_bar(aes(x=Year,y=nuRecs),stat="identity",width=0.75)+
  theme_classic()+
  xlab("Year")+
  scale_x_continuous(breaks=c(1950,1960,1970,1980,1990,2000,2010),
                     labels=c(1950,1960,1970,1980,1990,2000,2010))+
  geom_vline(xintercept=1980,colour="red",linetype="dashed")+
  ylab("Total number of records")

g2 <- ggplot(subset(timeSummary,Year>1950 & Year <2017))+
  geom_bar(aes(x=Year,y=nuPlots),stat="identity",width=0.75)+
  theme_classic()+
  xlab("Year")+
  scale_x_continuous(breaks=c(1950,1960,1970,1980,1990,2000,2010),
                     labels=c(1950,1960,1970,1980,1990,2000,2010))+
  geom_vline(xintercept=1980,colour="red",linetype="dashed")+
  ylab("Total number of survey quadrants")

library(cowplot)
grid2 <- plot_grid(g1,g2,ncol=1)
plot_grid(gMap,grid2,
          labels = c("A","B"),
          rel_heights = c(2,1))

ggsave(filename="plots/Fig1.png",width=6,height=5)

#Time series over time for each state

timeSummary <- ddply(adultData,.(Year,State),summarise,
                     nuRecs=length(Species),
                     nuSpecies=length(unique(Species)),
                     nuPlots=length(unique(MTB_Q)),
                     nuVisits=length(unique(Beobachter)))

#number of records
library(ggplot2)
q1 <- qplot(Year,nuRecs,data=subset(timeSummary,Year>1979),colour=State)+
  geom_line()+
  theme_bw()+
  facet_wrap(~State,scales="free")+
  theme(legend.position="none")+
  ylab("Number of records")
#ggsave(filename="plots/Adult_timeseries_effort_q1.png",width=5,height=4)

#number of species seen per year
q2 <- qplot(Year,nuSpecies,data=subset(timeSummary,Year>1979),colour=State)+
  geom_line()+
  theme_bw()+
  theme(legend.position="none")+
  ylab("Number of species")

#number of visits
q3 <- qplot(Year,nuVisits,data=subset(timeSummary,Year>1979),colour=State)+
  geom_line()+
  theme_bw()+
  theme(legend.position="none")+
  scale_y_log10()+
  ylab("Number of observers")

#number of sampling points
q4 <- qplot(Year,nuPlots,data=subset(timeSummary,Year>1979),colour=State)+
  geom_line()+
  theme_bw()+
  theme(legend.position="none")+
  scale_y_log10()+
  ylab("Number of MTBQs")

library(cowplot)
plot_grid(q1,q2,q3,q4)
ggsave(filename="plots/Adult_timeseries_effort.png",width=12,height=7)
ggsave(filename="plots/Juv_timeseries_effort.png",width=12,height=7)

#examine relationship among these variables
library(GGally)
timeSummary[,3:6]<-sapply(timeSummary[,3:6],log)
ggpairs(timeSummary[,3:6])

### species summary statistics ##############################################################################

#Coenagrion hylas, Gomphus simillimus, Lestes macrostigma, Onychogomphus uncatus were excluded

adultDataS <- subset(adultDataS,!Species %in% c("Coenagrion hylas", "Gomphus simillimus", "Lestes macrostigma", 
                                                "Onychogomphus uncatus"))

nuRecs$MTB_Q <- paste0(nuRecs$MTB,nuRecs$Q)
adultDataS <- subset(adultDataS, MTB_Q %in% nuRecs$MTB_Q[nuRecs$nuYears>1])

speciesSummary <- ddply(adultDataS,
                        .(Species),summarise,
                        nuRecs=length(Date),
                        nuGrids=length(unique(MTB_Q)),
                        nuYears=length(unique(Year)))

speciesSummary <- arrange(speciesSummary,desc(nuRecs))
speciesSummary

subset(speciesSummary,nuYears<20)
median(speciesSummary$nuYears)
summary(speciesSummary$nuYears)
write.csv(speciesSummary,file="derived-data/SpeciesSummaries.csv",row.names=FALSE)

speciesDetectionYears <- ddply(subset(adultData,Year>1980 & Year <2017),.(Species,Year),
                               summarise,
                               nuRecs=length(Date))

#expand to all years and species
newgrid <- expand.grid(Species = unique(adultData$Species),
                       Year = 1980:2016)

speciesDetectionYears <- merge(speciesDetectionYears,newgrid,all=T)

### geographic statistics ###################################################################################

#number of records each month, year and state


### revisit statistics ######################################################################################

#how many dates are there per MTBQ and per Year

reVisits <- ddply(subset(adultDataS,Year>1980 & Year <2017),.(MTB_Q,Year),
                               summarise,nuVisits=length(unique(Date)))

#how many sites are revisied at least once
(temp <- ddply(reVisits,.(Year),summarise,mean(nuVisits>1)))
summary(temp$..1)


temp <- ddply(reVisits,.(Year),summarise,
      mean(nuVisits[!is.na(nuVisits & nuVisits>1)]))

summary(temp$..1)

###species check##############################################################################################

df <- readRDS("df_sparta.rds")

table(df$Year)

#check the following
"Crocothemis erythraea"
"Ophiogomphus cecilia"
"Anax ephippiger"
"Sympetrum meridionale"

ddply(subset(df,Species=="Crocothemis erythraea"),.(Year),summarise,nuRecs=length(Species))#nothing until 1986
ddply(subset(df,Species=="Ophiogomphus cecilia"),.(Year),summarise,nuRecs=length(Species))#increasing
ddply(subset(df,Species=="Aeshna subarctica"),.(Year),summarise,nuRecs=length(Species))#increasing...
ddply(subset(df,Species=="Gomphus flavipes"),.(Year),summarise,nuRecs=length(Species))#increasing...
ddply(subset(df,Species=="Gomphus vulgatissimus"),.(Year),summarise,nuRecs=length(Species))#increasing then decreasing in recent years
ddply(subset(df,Species=="Leucorrhinia albifrons"),.(Year),summarise,nuRecs=length(Species))#increasing...
ddply(subset(df,Species=="Leucorrhinia pectoralis"),.(Year),summarise,nuRecs=length(Species))#increasing...
ddply(subset(df,Species=="Leucorrhinia rubicunda"),.(Year),summarise,nuRecs=length(Species))#maybe increasing...
ddply(subset(df,Species=="Leucorrhinia pectoralis"),.(Year),summarise,nuRecs=length(Species))#maybe increase

#delete Anax ephippiger and Sympetrum meridionale from the model
ddply(subset(df,Species=="Anax ephippiger"),.(Year),summarise,nuRecs=length(Species))#consistently rare. seen in only 7 years
ddply(subset(df,Species=="Sympetrum meridionale"),.(Year),summarise,nuRecs=length(Species))#consistently rare, few obs per year

#when there are few observations, the predictions are pulled upwards...


###below is old###################################################################
###State models############################################################################################

#get annual indices:

modelFiles <- list.files("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_modelSummary_flightperiod/5098409")
modelFiles <- modelFiles[grepl("modelSummary",modelFiles)]

#read in each one
library(plyr)
modelSummary <- ldply(modelFiles, function(x){
  temp <- readRDS(paste("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_modelSummary_flightperiod/5098409",x,sep="/"))
  temp$Folder <- x
  return(temp)
})

#get new files(baseed on effort) - for MV and NS
modelFiles <- c("modelSummary_effort_tests_Odonata_adult_MV.rds",
                "modelSummary_effort_tests_Odonata_adult_NS.rds")

#read in each one
library(plyr)
modelSummary2 <- ldply(modelFiles, function(x){
  temp <- readRDS(paste("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs",x,sep="/"))
  temp$Folder <- x
  return(temp)
})

#get ones based on number of species for effort
modelSummary2 <- modelSummary2[grepl("nuSpecies",modelSummary2$File),]
modelSummary2$Folder <- gsub("effort_tests_","",modelSummary2$Folder)
modelSummary <- rbind(modelSummary,modelSummary2)

#extract state and stage
modelSummary$Stage <- sapply(modelSummary$Folder,function(x)strsplit(as.character(x),"_")[[1]][3])
modelSummary$State <- sapply(modelSummary$Folder,function(x)strsplit(as.character(x),"_")[[1]][4])
modelSummary$State <- sapply(modelSummary$State,function(x)strsplit(as.character(x),"\\.")[[1]][1])

#extract species information
modelSummary$Species <- gsub("\\.rds","",modelSummary$File)
modelSummary$Species <- gsub("out_nuSpecies_","",modelSummary$Species)
modelSummary$Species <- gsub("outSummary_nuSpecies_","",modelSummary$Species)
modelSummary$Species <- gsub("out_flightperiod_nuSpecies_","",modelSummary$Species)
#modelSummary$Species <- gsub("juv_","",modelSummary$Species)
modelSummary$Species <- gsub("adult_","",modelSummary$Species)
modelSummary$Species <- sapply(modelSummary$Species,function(x){
  if(grepl("_",x)){
  strsplit(as.character(x),"_")[[1]][2]
  }else {
    x
  }
})

#in total 73 species
unique(modelSummary$Species)
#unique(modelSummary$Species[modelSummary$Stage=="juv"])
head(modelSummary)

#format state information
modelSummary$kState <- modelSummary$State
modelSummary$State <- as.factor(modelSummary$State)
levels(modelSummary$State)<-c("Bavaria","Brandenburg","Mecklenburg-Vorpommern",
                              "North Rhine-Westphalia","Lower Saxony",
                              "Rheinland-Pfalz","Saarland",
                              "Saxony-Anhalt","Saxony",
                              "Schleswig Holstein","Thuringia")
#check Rhat

out <- subset(modelSummary,Rhat > 1.1)
table(out$State)
table(out$Species)
table(out$Param)

#standardize species names

species<- read.delim("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/specieslist_odonata.txt")
modelSpecies <- unique(modelSummary$Species)
modelSpecies[!modelSpecies%in%species$Species]

###Habitat directive#############################################################################################

#how many habitat directive species

#vulnerable or near threatened as per EU Red List
notLC<-c("Aeshna viridis","Coenagrion armatum","Coenagrion hylas","Coenagrion mercuriale",
         "Coenagrion ornatum","Cordulegaster bidentata","Gomphus simillimus","Lestes macrostigma",        
         "Leucorrhinia albifrons","Leucorrhinia caudalis","Oxygastra curtisii",
         "Sympetrum depressiusculum")
hd <- c("Aeshna viridis","Coenagrion hylas","Coenagrion mercuriale","Coenagrion ornatum","Gomphus flavipes",
        "Leucorrhinia albifrons","Leucorrhinia caudalis","Leucorrhinia pectoralis","Ophiogomphus cecilia",
        "Oxygastra curtisii","Sympecma paedisca") 
all<-c(notLC,hd)

#Restrict to occcurence indicies:

modelSummary <- modelSummary[grepl("psi.fs",modelSummary$Param),]
modelSummary$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", modelSummary$Param))

#get start year for each state
modelFiles <- list.files("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects")
modelFiles <- modelFiles[grepl("yearDF_adult",modelFiles)]
yearDF <- ldply(modelFiles, function(x){
  temp <- read.delim(paste("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects",x,sep="/"))
  temp$Folder <- x
  return(temp)
})
modelSummary$StartYear <- min(yearDF$Year[match(modelSummary$kState,yearDF$State)])
modelSummary$Year <- modelSummary$StartYear + modelSummary$ParamNu - 1

###Species code###############################################################################################

#simplify species names to 3-letter code:

modelSummary$Genus <- sapply(modelSummary$Species, function(x)strsplit(x," ")[[1]][1])
modelSummary$Genus <- sapply(modelSummary$Genus, function(x) substr(x,1,3))
modelSummary$Spec <- sapply(modelSummary$Species, function(x)strsplit(x," ")[[1]][2])
modelSummary$Spec <- sapply(modelSummary$Spec, function(x) substr(x,1,3))
modelSummary$Code <- paste(modelSummary$Genus,modelSummary$Spec,sep="_")
  
####State-species#############################################################################################

#plot them at the species level per state:

library(ggplot2)

#for each stage and state??
ggplot(subset(modelSummary,kState=="BB" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("Brandenberg - Adults")

ggplot(subset(modelSummary,kState=="Bav" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("Bavaria - Adults")

ggplot(subset(modelSummary,kState=="RLP" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("RLP - Adults")
 
ggplot(subset(modelSummary,kState=="SH" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("SH - Adults")

ggplot(subset(modelSummary,kState=="Sax" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("Saxony - Adults")

ggplot(subset(modelSummary,kState=="Thu" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("Thuringia - Adults")

ggplot(subset(modelSummary,kState=="SAnhalt" & Stage=="adult"))+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")+
  ggtitle("Sachsen Anhalt - Adults")

###State trends###########################################################################################

#how many are changing in each lander
load("derived-data/trendEstimates.RData")
library(reshape2)
library(ggplot2)
library(plyr)

trendEstimates$state <- NA
trendEstimates$state[which(trendEstimates$lowerCI>0 & trendEstimates$upperCI>0)] <- "increase"
trendEstimates$state[which(trendEstimates$lowerCI<0 & trendEstimates$upperCI<0)] <- "decrease"
trendEstimates$state[is.na(trendEstimates$state)]<-"stable"

trendSummary <- ddply(trendEstimates,.(Stage,State),summarise,
                      increase = sum(state=="increase"),
                      decrease = sum(state=="decrease"),
                      stable = sum(state=="stable"))
trendSummary<-melt(trendSummary,id=c("Stage","State"))
names(trendSummary)[which(names(trendSummary)=="variable")]<-"Trend"
trendSummary$Trend <- factor(trendSummary$Trend, levels=c("decrease","increase","stable"))

ggplot(trendSummary)+
  geom_bar(aes(x=State,y=value,fill=Trend),stat="identity")+
  facet_grid(~Stage)+
  ylab("number of species")+
  theme_bw()+
  coord_flip()

#Change the labels to German
levels(trendSummary$State)<-c("Bayern","Nordrhein-Westfalen","Rheinland-Pfalz" ,"Saarland",
                              "Sachsen-Anhalt","Saschen","Schleswig-Holstein","Thuringia")

levels(trendSummary$Trend) <- c("Abnahme","Zunahme","stabil")

ggplot(subset(trendSummary,Stage=="adult"))+
  geom_bar(aes(x=State,y=value,fill=Trend),stat="identity")+
  ylab("Anzahl Libellen Arten")+
  xlab("")+
  theme_bw()+
  coord_flip()

ggsave(file="German_Trends.tiff",dpi=300,width=5,height=3)

summary(lm(trend~Stage+State,data=trendEstimates))

####State comparison###########################################################################################

#pairwise comparison across lander

load("derived-data/trendEstimates.RData")

sortTrends <- dcast(trendEstimates,Species+Stage~State,value.var="trend")

#restrict to species seen in all lander
sortTrends <- subset(sortTrends,complete.cases(sortTrends))
sortTrends <- subset(sortTrends,Stage=="adult")
nrow(sortTrends)#35

#want to find out consistentl declining and consistently increasing
sortTrends$medianOccu <- apply(sortTrends[,3:7],1,median)
sortTrends$sdOccu <- apply(sortTrends[,3:7],1,sd)
sortTrends$maxOccu <- apply(sortTrends[,3:7],1,max)
sortTrends$minOccu <- apply(sortTrends[,3:7],1,min)
sortTrends$covOccu <- sortTrends$medianOccu/sortTrends$sdOccu
sortTrends$rangeOcc <- sortTrends$maxOccu-sortTrends$minOccu

out <- arrange(sortTrends,sdOccu)

#plot all trends on the same graphs

#order species by their sd of trend
trendEstimatesSD <- subset(trendEstimates,Species %in% sortTrends$Species & Stage=="adult")
trendEstimatesSD$Species <- factor(trendEstimatesSD$Species, level=rev(out$Species))

ggplot(trendEstimatesSD)+
  geom_point(aes(x=Species,y=trend,colour=State,size=1/trend_sd),alpha=0.5)+
  coord_flip()+
  geom_hline(yintercept = 0, colour="black", linetype="dashed")+
  theme_bw()
  
####Adult vs Juv############################################################################################

#relationship between adult and juvenile trends

load("trendEstimates.RData")

sortTrends <- dcast(trendEstimates,Species+State~Stage,value.var="trend")
sortTrends2 <- dcast(trendEstimates,Species+State~Stage,value.var="trend_sd")
names(sortTrends2)[3:4]<-c("adult_sd","juv_sd")
sortTrends <- cbind(sortTrends,sortTrends2[,3:4])
sortTrends$sd <- sortTrends$adult_sd + sortTrends$juv_sd

ggplot(sortTrends,aes(x=adult,y=juv))+
  geom_point()+
  facet_wrap(~State,scales="free")+
  stat_smooth(method="lm")+
  theme_bw()

####Community matrix##########################################################################################

#community samples

#drawn 1000 possible communities for each year

myMatrix <- matrix(data=NA,nrow=nrow(annualDF),ncol=1000) 
for(j in 1:ncol(myMatrix)){
  for(i in 1:nrow(myMatrix)){
  myp <- rnorm(n=1,mean=annualDF$mean[i],sd=annualDF$sd[i]) 
  myp[myp>1]<-1#to avoid boundary effects
  myp[myp<0]<-0#to avoid boundary effects
  #myMatrix[i,j] <- rbinom(1,1,myp)
  myMatrix[i,j] <- myp
  }
}
randomMatrix<-cbind(annualDF[,c("Year","Species")],myMatrix)

save(randomMatrix,file="randomMatrix.RData")
                                            
###Euro trends#########################################################################################

#plot european trends

euroTrends <- read.delim("euroTrends.txt",as.is=T)
library(ggplot2)
ggplot(euroTrends)+
  geom_bar(aes(x=country,y=nuSpecies,fill=Trend),stat="identity")+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("number of species")

###Thinning#########################################################################################

#thinning records
z <- coda::as.mcmc.list(out$samples)
z2 <- window(z, start=601, end=1000)
summary(z2)

gelman.diag(z2,multivariate=FALSE)
#Rhat is the potential scale reduction factor (at convergence, Rhat=1).

####Spline########################################################################################

#spline analysis
splineData <- readRDS("model-outputs/outSummary_dynamicspline_adult_Crocothemis erythraea2.rds")
summary(splineData$Rhat)
splineData$Param <- row.names(splineData)
splineData$ParamNu <- sub(".*\\[([^][]+)].*", "\\1", splineData$Param)
splineData$Site <- as.numeric(sapply(splineData$ParamNu,function(x)strsplit(x,",")[[1]][1]))
splineData$Year <- as.numeric(sapply(splineData$ParamNu,function(x)strsplit(x,",")[[1]][2]))
splineData <- subset(splineData,Param!="deviance")
summary(splineData$mean)

#get MTBQ data
load("siteInfo.RData")
splineData$MTBQ <- siteInfo$MTB_Q[match(splineData$Site,siteInfo$siteIndex)]

fill.na <- function(x, i=13) {
  if( is.na(x)[i] ) {
    return( median(x, na.rm=TRUE) )
  } else {
    return( x[i] )
  }
}

library(ggplot2)
#library(rasterVis)
library(viridis)

#for last year
for(i in 1:37){
splineData20 <- subset(splineData,Year==i)
summary(splineData20$mean)
mean(splineData20$mean==1)
mtbqDF$Occupancy <- splineData20$mean[match(mtbqDF$MTB_Q,splineData20$MTBQ)]

#as raster
pixels <- SpatialPixelsDataFrame(points=mtbqDF[,c('x','y')], 
                                 data=mtbqDF[,c('Q','MTB_Q','Occupancy')],
                                 tolerance = 0.916421)
r <- raster(pixels[,'Occupancy'])
#plot(r)

#interpolate over missing values
#https://gis.stackexchange.com/questions/181011/fill-the-gaps-using-nearest-neighbors
r2 <- focal(r, w = matrix(1,5,5), fun = fill.na, 
            pad = TRUE, na.rm = FALSE )

#remove 1
r2[r2==1] <- NA
r3 <- focal(r2, w = matrix(1,5,5), fun = fill.na, 
            pad = TRUE, na.rm = FALSE )

#mask cells outside germany
r3 <- mask(r3,germanyMap)

r2F <- as.data.frame(r3,xy=T)
names(r2F)[3] <- "Occupancy"
ggplot(r2F) +
  geom_tile(aes(x,y,fill=Occupancy)) +
  scale_fill_viridis(limits=c(0,1),na.value="white")+
  theme_void()+
  ggtitle(paste(i+1979))

ggsave(paste0("gifs/year",i,".tiff"),width=5,height=5)

}

#https://stackoverflow.com/questions/1298100/creating-a-movie-from-a-series-of-plots-in-r
#https://ryanpeek.Odonata_Github.io/2016-10-19-animated-gif_maps_in_R/

###GIF#########################################################################################
#make an animated gifs
library(magick)
library(purrr) 

list.files("gifs/") %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write("thermo.gif") # write to current dir

system("/opt/local/bin/convert -delay 80 *.png example_1.gif")

animation <- image_animate(img, fps = 2)
write.gif()

