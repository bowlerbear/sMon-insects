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


