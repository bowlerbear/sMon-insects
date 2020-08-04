#stageEffects

#get data
load("derived-data/juv_datafile.RData")
load("derived-data/adult_datafile.RData")
adult_datafile$Anzahl_min <- as.numeric(as.character(adult_datafile$Anzahl_min))

#combine
juv_datafile$Stage <- "Juvenile"
adult_datafile$Stage <- "Adult"
datafile <- rbind(juv_datafile,adult_datafile)


#look at how often both were seen in the same visit
datafile <- subset(datafile,!is.na(Date) & !is.na(MTB_Q))
datafile$visit <- paste0(datafile$Date,datafile$MTB_Q)

#summaruse it
library(plyr)
visitSummary <- ddply(datafile,.(visit),summarise,
                      Juv= sum(Anzahl_min[Stage=="Juvenile"]>0,na.rm=T),
                      Adult= sum(Anzahl_min[Stage=="Adult"]>0,na.rm=T), 
                      Both= sum(Anzahl_min>0,na.rm=T))
visitSummary<-subset(visitSummary,Both>0)

#get detection probs for juvenile or adult
visitSummary$Juv_B <- ifelse(visitSummary$Juv>0,1,0) 
visitSummary$Adult_B <- ifelse(visitSummary$Adult>0,1,0) 

hist(visitSummary$Jprop)
hist(visitSummary$Adult)
hist(visitSummary$Juv)

hist(visitSummary$Juv_B/(visitSummary$Juv_B + visitSummary$Adult_B))
#distribution is not symmetrical

#either/or....
#juvenile visits shouldnt be treated as adult visits

Look at availability of each data type


#total records per year
#total species per year
#total grids per year
#total observers per year
#total visits per observer per year

library(plyr)
library(lubridate)

load("derived-data/juv_datafile.RData")
juv_datafile$Date<-as.Date(juv_datafile$Date)
juv_datafile$Year<-year(juv_datafile$Date)
juv_datafile
juv_datafile<-subset(juv_datafile,Anzahl_min>0)
juvSummary <- ddply(subset(juv_datafile,Year>1979),.(Year),summarise,
                    nuRecords=length(Species),
                    nuSpecies=length(unique(Species)),
                    nuGrids=length(unique(MTB_Q)),
                    nuObservers=length(unique(Expertise)),
                    nuVisits=length(unique(interaction(Date,MTB_Q,Expertise))))
juvSummary$Stage <- "Juvenile"


load("derived-data/adult_datafile.RData")
adult_datafile$Date<-as.Date(adult_datafile$Date)
adult_datafile$Year<-year(adult_datafile$Date)
adult_datafile$Anzahl_min <- as.numeric(as.character(adult_datafile$Anzahl_min))
adult_datafile<-subset(adult_datafile,Anzahl_min>0)
adultSummary <- ddply(subset(adult_datafile,Year>1979),.(Year),summarise,
                      nuRecords=length(Species),
                      nuSpecies=length(unique(Species)),
                      nuGrids=length(unique(MTB_Q)),
                      nuObservers=length(unique(Expertise)),
                      nuVisits=length(unique(interaction(Date,MTB_Q,Expertise))))
adultSummary$Stage<- "Adult"

#combine all
summaryData<-rbind(juvSummary,adultSummary)


g1<-ggplot(data=summaryData,aes(x=Year,y=nuRecords,group=Stage))+geom_line(aes(colour=Stage))+ylab("total records")
g2<-ggplot(data=summaryData,aes(x=Year,y=nuSpecies,group=Stage))+geom_line(aes(colour=Stage))+ylab("total species")
g3<-ggplot(data=summaryData,aes(x=Year,y=nuGrids,group=Stage))+
  geom_line(aes(colour=Stage))+ylab("total grids")
g4<-ggplot(data=summaryData,aes(x=Year,y=nuObservers,group=Stage))+geom_line(aes(colour=Stage))+ylab("total observers")

g5<-ggplot(data=summaryData,aes(x=Year,y=nuVisits,group=Stage))+geom_line(aes(colour=Stage))+ylab("total observer visits")

library(cowplot)

plot_grid(g1,g5,g2,g3)


```




