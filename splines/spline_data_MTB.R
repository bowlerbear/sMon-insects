library(tidyverse)

### load in the full dataset ####

adultData <- readRDS(paste("splines","adultData.rds",sep="/"))
load(paste(getwd(),"mtbqsDF.RData",sep="/"))
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

# subset  

df <- subset(adultData, Year>=1990  & Year<2017)
#table(df$Month)
df <- subset(df, Month %in% 4:10)

#get species

speciesTaskID <- read.delim(paste("model-auxfiles","speciesTaskID_adult.txt",sep="/"),as.is=T)
allspecies <- speciesTaskID$Species

#get MTB
df$MTB <- sapply(as.character(df$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,(len-1))})

### function to get data for a species ####

all(df$MTB %in% mtbqsDF$MTB)

df$State <- mtbsDF$Counties[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$State))

#check missing naturruam data
mtbsDF$Natur[is.na(mtbsDF$Natur)] <- mtbsDF$MTB_Natur[is.na(mtbsDF$Natur)]
df$Natur <- mtbsDF$Natur[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$Natur))

mtbsDF$CoarseNatur[is.na(mtbsDF$CoarseNatur)] <- mtbsDF$MTB_CoarseNatur[is.na(mtbsDF$CoarseNatur)]
df$CoarseNatur <- mtbsDF$CoarseNatur[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$CoarseNatur))

mtbsDF$MTB_MidNatur <- gdata::trim(mtbsDF$MTB_MidNatur)
df$MidNatur <- mtbsDF$MTB_MidNatur[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$MidNatur))

#define a visit
df$visit <- paste(df$MTB,df$Date,df$Beobachter,sep="_")

#get occurence matrix  - detection and non-detection
getOccurrenceMatrix<-function(df){
  out <- reshape2::acast(df,visit~Species,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
  out[out>0]<-1
  return(out)
}
occMatrix <- getOccurrenceMatrix(df)

#get list length
getListLength<-function(df){
  df %>% 
    group_by(visit,Date,MTB) %>%
    summarise(nuSpecies=length(unique(Species)),
                    nuRecords=length(Species)) %>%
    arrange(.,visit)
}
listlengthDF <- getListLength(df)

#check rows of occuMatrix match visits
all(listlengthDF$visit==row.names(occMatrix))

#add on some indices
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- lubridate::year(listlengthDF$Date)
listlengthDF$yday <- lubridate::yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
subset(listlengthDF,is.na(Naturraum))

listlengthDF$MidNaturraum <- mtbqsDF$MTB_MidNatur[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$mnIndex <- as.numeric(factor(listlengthDF$MidNaturraum))
subset(listlengthDF,is.na(MidNaturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)


### function for each species ####

getSpData <- function(myspecies){
  
#add on species data to the listlength object
listlengthDF$Species <- occMatrix[,myspecies]

listlengthDF %>%
  group_by(MTB,Year) %>%
  summarise(Detected = max(Species, na.rm=T),
            nuVisits = length(Species)) %>%
  ungroup %>%
  add_column(Species = myspecies)

}

spData <- lapply(allspecies, getSpData)
spData <- do.call(rbind,spData)
spData$Detected <- ifelse(spData$Detected==1,"Yes","No") 

#add on species coordinates
spData$x_MTB <- mtbsDF$x_MTB[match(spData$MTB, mtbsDF$MTB)]
spData$y_MTB <- mtbsDF$y_MTB[match(spData$MTB, mtbsDF$MTB)]

# plot maps

myYears <- sort(unique(spData$Year))


library(scales)

for(s in allspecies){
  for(i in 1:length(myYears)){
    
    ggplot(subset(spData, Species == s & Year == myYears[i]))+
      geom_point(aes( x=x_MTB, y = y_MTB, colour = nuVisits),size = 2.5, shape=15)+
      scale_color_viridis_c("Number of Visits",option = "C", direction = -1,
                            trans = log2_trans(), limits=c(1,324))+
      geom_point(data = subset(spData, Detected =="Yes" & Species == s & Year == myYears[i]),
                 aes(x=x_MTB, y = y_MTB), size = 3.5, shape=4)+
      theme_void()
    
    ggsave(file=paste0("plots/species/Data_",myYears[i],"_",s,".png"), width=5.5, height=6)    
    
  }
}
