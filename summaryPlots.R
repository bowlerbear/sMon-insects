#summaryPlots
library(rgdal)
library(ggplot2)

###################################################################################################

#species being analysed
mySpecies <- read.delim("speciesTaskID_adult.txt",as.is=T)$Species
################################################################################################

library(maptools)

#get map of Germany
germanyMap <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")

#MTBQ
mtbqMap <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                   layer="MTBQ_25833")
mtbqDF <- data.frame(mtbqMap@data,x=coordinates(mtbqMap)[,1],y=coordinates(mtbqMap)[,2])
mtbqDF$Q <- NA
mtbqDF$Q[which(mtbqDF$Quadrant=="NW")]<-1
mtbqDF$Q[which(mtbqDF$Quadrant=="NO")]<-2
mtbqDF$Q[which(mtbqDF$Quadrant=="SW")]<-3
mtbqDF$Q[which(mtbqDF$Quadrant=="SO")]<-4

#MTB
mtbMap <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                  layer="MTB_25832")

#Get crs
proj4string(mtbqMap)
proj4string(mtbMap)
mtbMap <- spTransform(mtbMap,CRS(proj4string(mtbqMap)))

#Overlay
plot(mtbMap)
plot(mtbqMap,add=T,col="red")
#both in "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"

#################################################################################################

#get adult records
myfiles <- list.files("derived-data")

#read in and combine all adult files
#adultFiles<-myfiles[grepl("juv",myfiles)]
adultFiles<-myfiles[grepl("adult",myfiles)]
adultFiles<-adultFiles[grepl("RData",adultFiles)]

#combine these files
library(plyr)
adultData <- ldply(adultFiles,function(x){
  load(paste("derived-data",x,sep="/"))
  adult_datafile$File <- x
  return(adult_datafile)
})

#adultData <- ldply(adultFiles,function(x){
#  load(paste("derived-data",x,sep="/"))
#  juv_datafile$File <- x
#  return(juv_datafile)
#})

#extract state from file name
adultData$State <- sapply(adultData$File,function(x)strsplit(x,"\\.RData")[[1]][1])
adultData$State <- sapply(adultData$State,function(x)strsplit(x,"_")[[1]][3])


library(lubridate)
adultData$Date<-as.Date(adultData$Date)
adultData$Year <- year(adultData$Date)

#change state labels
adultData$State <- as.factor(adultData$State)
levels(adultData$State) <- c("Bavaria","Brandenberg","Baden-Wuttemberg","Hesse",
                             "Mecklenburg-Vorpommern",
                             "North Rhine-Westphalia","Lower Saxony",
                             "Rheinland-Pfalz","Saarland","Saxony-Anhalt",
                             "Saxony","Schleswig Holstein","Thuringia")

library(lubridate)
adultData$Date <- as.Date(adultData$Date)
adultData$yday <- yday(adultData$Date)
adultData$week <- week(adultData$Date)
adultData$Month <- month(adultData$Date)
adultData <- subset(adultData,!is.na(Date))

#######################################################################################

#fix names
species<- read.delim("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/specieslist_odonata.txt")
adultSpecies <- sort(unique(adultData$Species))
adultSpecies[!adultSpecies%in%species$Species]

#############################################################################################

#how often is each site visited
out <- ddply(adultData,.(Year,MTB_Q),summarise,nu = length(unique(yday)))
summary(out$nu)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   2.958   3.000 106.000

#how many sites is a species usually present at
out <- ddply(adultData,.(Species),summarise,nu = length(unique(MTB_Q)))
length(unique(adultData$MTB_Q))#8231
summary(out$nu/8231)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0001215 0.0297048 0.1132305 0.1974859 0.3422427 0.7307739

#on how many visits is a species usually detected
out <- ddply(adultData,.(MTB_Q,Year),summarise,nuVisits=length(unique(yday)))
out2 <- ddply(adultData,.(Species,MTB_Q,Year),summarise,nuDets=length(unique(yday)))
out <- merge(out,out2,by=c("MTB_Q","Year"))
out$p <- out$nuDets/out$nuVisits
summary(out$p)

#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.009434 0.250000 0.500000 0.522878 1.000000 1.000000

out <- ddply(out,.(Species),summarise,meanp=mean(p))
summary(out$meanp)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3789  0.4509  0.4941  0.5140  0.5519  0.7560

###############################################################################################

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

###############################################################################################

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
  
###############################################################################################

#plot number of records per MTB

adultDataS <- adultData

#extract MTB from Q
adultDataS$MTB <- gsub("/","",adultDataS$MTB_Q)

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
nuRecs <- ddply(subset(adultDataS,Year>1979),.(MTB,Q),
                summarise,nuRecs = length(Species))


#plot
allDF <- merge(mtbqDF,nuRecs,by.x=c("Value","Q"),by.y=c("MTB","Q"))
nrow(mtbqDF)
nrow(allDF)
#data from 4186/12024

#plot map of germany
library(ggplot2)
germanyMap <- spTransform(germanyMap,proj4string(mtbMap))
AG <- fortify(germanyMap)
ggplot()+ geom_polygon(data=AG, aes(long, lat, group = group), colour = "grey", fill=NA)+
  geom_point(data=allDF,aes(x=x,y=y,colour=log(nuRecs)))+
  scale_colour_gradient(low="pink",high="darkred")+
  xlab("X")+ylab("Y")+
  theme_void()
ggsave(filename="plots/Adult_map_effort.png",width=8,height=7)
ggsave(filename="plots/Juv_map_effort.png",width=8,height=7)

################################################################################################

#Time series over time

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
  scale_y_log10()+
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

###############################################################################################

#examine relationship among these variables
library(GGally)
timeSummary[,3:6]<-sapply(timeSummary[,3:6],log)
ggpairs(timeSummary[,3:6])

###############################################################################################

#get annual indices:

modelFiles <- list.files("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/Odonata_modelSummary_flightperiod/5098409")
modelFiles <- modelFiles[grepl("modelSummary",modelFiles)]

#read in each one
library(plyr)
modelSummary <- ldply(modelFiles, function(x){
  temp <- readRDS(paste("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/Odonata_modelSummary_flightperiod/5098409",x,sep="/"))
  temp$Folder <- x
  return(temp)
})

#############################################################################################

#get new files(baseed on effort) - for MV and NS
modelFiles <- c("modelSummary_effort_tests_Odonata_adult_MV.rds",
                "modelSummary_effort_tests_Odonata_adult_NS.rds")

#read in each one
library(plyr)
modelSummary2 <- ldply(modelFiles, function(x){
  temp <- readRDS(paste("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs",x,sep="/"))
  temp$Folder <- x
  return(temp)
})

#get ones based on number of species for effort
modelSummary2 <- modelSummary2[grepl("nuSpecies",modelSummary2$File),]
modelSummary2$Folder <- gsub("effort_tests_","",modelSummary2$Folder)
modelSummary <- rbind(modelSummary,modelSummary2)

#################################################################################################

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

##############################################################################################

#format state information
modelSummary$kState <- modelSummary$State
modelSummary$State <- as.factor(modelSummary$State)
levels(modelSummary$State)<-c("Bavaria","Brandenburg","Mecklenburg-Vorpommern",
                              "North Rhine-Westphalia","Lower Saxony",
                              "Rheinland-Pfalz","Saarland",
                              "Saxony-Anhalt","Saxony",
                              "Schleswig Holstein","Thuringia")

#############################################################################################
#check Rhat

out <- subset(modelSummary,Rhat > 1.1)
table(out$State)
table(out$Species)
table(out$Param)

###############################################################################################

#standardize species names

species<- read.delim("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/specieslist_odonata.txt")
modelSpecies <- unique(modelSummary$Species)
modelSpecies[!modelSpecies%in%species$Species]

################################################################################################

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

#################################################################################################

#Restrict to occcurence indicies:

modelSummary <- modelSummary[grepl("psi.fs",modelSummary$Param),]
modelSummary$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", modelSummary$Param))

#get start year for each state
modelFiles <- list.files("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects")
modelFiles <- modelFiles[grepl("yearDF_adult",modelFiles)]
yearDF <- ldply(modelFiles, function(x){
  temp <- read.delim(paste("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects",x,sep="/"))
  temp$Folder <- x
  return(temp)
})
modelSummary$StartYear <- min(yearDF$Year[match(modelSummary$kState,yearDF$State)])
modelSummary$Year <- modelSummary$StartYear + modelSummary$ParamNu - 1

##################################################################################################

#simplify species names to 3-letter code:

modelSummary$Genus <- sapply(modelSummary$Species, function(x)strsplit(x," ")[[1]][1])
modelSummary$Genus <- sapply(modelSummary$Genus, function(x) substr(x,1,3))
modelSummary$Spec <- sapply(modelSummary$Species, function(x)strsplit(x," ")[[1]][2])
modelSummary$Spec <- sapply(modelSummary$Spec, function(x) substr(x,1,3))
modelSummary$Code <- paste(modelSummary$Genus,modelSummary$Spec,sep="_")
  
#################################################################################################

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

################################################################################################

#get average across states

nationSummary <- ddply(subset(modelSummary, Stage=="adult"),.(Code,Species,Year),summarise,
                       meanMean = median(mean,na.rm=T),
                       meanSD = mean(sd))
nationSummary$lowerCI <- nationSummary$meanMean - 1.96*nationSummary$meanSD
nationSummary$upperCI <- nationSummary$meanMean + 1.96*nationSummary$meanSD

ggplot(nationSummary)+
  geom_line(aes(x=Year,y=meanMean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+ ylab("Occupancy")

################################################################################################

#cluster them into groups of species with similar dyanamics

modelSummary2 <- nationSummary
library(reshape2)

#######################
#cluster on occurences#
#######################

mydata <- acast(modelSummary2,Species~Year,value.var="meanMean")
mydata[is.na(mydata)]<-0
#apply cluster analysis
d <- dist(mydata, method = "euclidean") 
fit <- hclust(d, method="ward.D2")
plot(fit) # display dendogram

wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

#pam
library(cluster)
fit <- kmeans(mydata, 6) 

# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster) 
mydata$Species <- row.names(mydata)

#add to dataframe 
modelSummary2$cluster <- mydata$fit.cluster[match(modelSummary2$Species,mydata$Species)]

#order
modelSummary2$cluster <- as.factor(modelSummary2$cluster)
modelSummary2$cluster <- factor(modelSummary2$cluster, levels=c("3","2","6","1","4","5"))
levels(modelSummary2$cluster) <- c("1","2","3","4","5","6")
ddply(modelSummary2,.(cluster),summarise,nuS = length(unique(Species)))
levels(modelSummary2$cluster) <- c("16 sp","13 sp","11 sp","6 sp","10 sp","17 sp")

#plot
ggplot(modelSummary2)+
  geom_line(aes(x=Year,y=meanMean,colour=Species))+
  facet_wrap(~cluster)+
  theme_bw()+ ylab("Occupancy")+
  theme(legend.position="none")


##########################
#cluster on growth rates##
##########################

mydata <- ddply(modelSummary2,.(Species),function(x){
  len = length(x$Year)
  growth = (x$meanMean[2:len]/(1-x$meanMean[2:len]))/(x$meanMean[1:(len-1)]/(1-x$meanMean[1:(len-1)]))
  data.frame(Species=rep(unique(x$Species),(len-1)),growth,Year=2:len)
})
mydata <- acast(mydata,Species~Year,value.var="growth")
mydata[is.na(mydata)]<-0

#apply cluster analysis
d <- dist(mydata, method = "euclidean") 
fit <- hclust(d, method="ward.D2")
plot(fit) # display dendogram

wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

#choose cluster
fit <- kmeans(mydata, 4) 

# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster) 
mydata$Species <- row.names(mydata)

#add to dataframe 
modelSummary2$cluster <- mydata$fit.cluster[match(modelSummary2$Species,mydata$Species)]

#plot
ggplot(modelSummary2)+
  geom_line(aes(x=Year,y=meanMean,colour=Species))+
  facet_wrap(~cluster)+
  theme_bw()+ ylab("Occupancy")+
  theme(legend.position="none")

################################################################################################

#get population trends per region and species
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/sparta_wrapper_functions.R')

trendEstimates <- ddply(modelSummary,.(Species,State,Stage),
                        function(x){
                          fitTrends(x)
                        })
save(trendEstimates,file="derived-data/trendEstimates.RData")

trendEstimates <- ddply(subset(modelSummary,Stage=="adult"),.(Species,State),
                        function(x){
                          fitTrendsBeta(x)
                        })

save(trendEstimates,file="trendEstimatesBeta.RData")

##############################################################################################

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

#############################################################################################

#national trends

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/sparta_wrapper_functions.R')
names(nationSummary)[4:5]<-c("mean","sd")

trendEstimates <- ddply(nationSummary,.(Species),
                        function(x){
                          fitTrends(x)
                        })

save(trendEstimates,file="derived-data/trendEstimatesNational.RData")

#plot as histogram
head(trendEstimates)
trendEstimates$state <- NA
trendEstimates$state[which(trendEstimates$lowerCI>0 & trendEstimates$upperCI>0)] <- "sig. increase"
trendEstimates$state[which(trendEstimates$lowerCI<0 & trendEstimates$upperCI<0)] <- "sig. decrease"
trendEstimates$state[is.na(trendEstimates$state)]<-"non sig"
trendEstimates$change <- factor(trendEstimates$state,levels=c("sig. decrease","sig. increase","non sig"))

ggplot(trendEstimates)+
  geom_histogram(aes(trend,fill=change))+
  theme_bw()+
  geom_vline(xintercept=0,colour="black",linetype="dashed")

###############################################################################################

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
  
################################################################################################

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

##############################################################################################

#community samples

#drawn 1000 possible communities for each year

modelSummary_Ad <- subset(modelSummary, Stage=="adult")

#get matrix of possible communities for each year and species and lander
myMatrix <- matrix(data=NA,nrow=nrow(modelSummary_Ad),ncol=1000) 
for(j in 1:ncol(myMatrix)){
  for(i in 1:nrow(myMatrix)){
  myp <- dnorm(1,modelSummary_Ad$mean[i],modelSummary_Ad$sd[i]) 
  myp[myp>1]<-0.9999
  myp[myp<0]<-0.0001
  myMatrix[i,j] <- rbinom(1,1,myp)
  }
}
randomMatrix<-cbind(modelSummary_Ad[,c("State","Year","Species")],myMatrix)
save(randomMatrix,file="randomMatrix.RData")
                                                 
#get species richness for each community
out <- ddply(randomMatrix,.(State,Year),function(x){
  numcolwise(sum)(x)})

#get mean and 95% CI across communities
out$meanRichness<-apply(out[,3:1002],1,median)
out$lowerRichness<-apply(out[,3:1002],1,function(x)quantile(x,0.025))
out$upperRichness<-apply(out[,3:1002],1,function(x)quantile(x,0.975))

out$State <- as.factor(out$State)
levels(out$State)<-c("Bavaria","North Rhine-Westphalia","Saarland","Saxony","Schleswig Holstein")

ggplot(out)+
  geom_line(aes(x=Year,y=meanRichness,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Average species richness")+
  facet_grid(~State)

#plot disimmilarity
library(vegan)
randomMatrixM<-melt(randomMatrix,id=c("State","Year","Species"))
randomMatrixM<-arrange(randomMatrixM,Species,Year,State)

out <- ddply(randomMatrixM,.(State,variable),function(x){
  moo2 <- acast(x,Year~Species,value.var="value")
  diss <- as.numeric(as.matrix(vegdist(moo2,binary=TRUE))[,1])
  year <- names(as.matrix(vegdist(moo2,binary=TRUE))[,1])
  data.frame(year,diss)
})

out <- subset(out,year!=1981)

#take average across all
out <- ddply(out,.(State,year),summarise,
             meanD = mean(diss),
             lowerCI = quantile(diss,0.025),
             upperCI = quantile(diss,0.975))

out$State <- as.factor(out$State)
levels(out$State)<-c("Bavaria","North Rhine-Westphalia","Saarland","Saxony","Schleswig Holstein")
out$year <-as.numeric(as.character(out$year))
ggplot(out)+
  geom_line(aes(x=year,y=meanD,colour=State))+
  geom_ribbon(aes(x=year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Average dissimilarity")+
  facet_grid(~State)

############################################################################################

#plot european trends

euroTrends <- read.delim("euroTrends.txt",as.is=T)
library(ggplot2)
ggplot(euroTrends)+
  geom_bar(aes(x=country,y=nuSpecies,fill=Trend),stat="identity")+
  coord_flip()+
  theme_bw()+
  xlab("")+ylab("number of species")

############################################################################################

#thinning records
z <- coda::as.mcmc.list(out$samples)
z2 <- window(z, start=601, end=1000)
summary(z2)

gelman.diag(z2,multivariate=FALSE)
#Rhat is the potential scale reduction factor (at convergence, Rhat=1).

############################################################################################