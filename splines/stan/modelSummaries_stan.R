library(tidyverse)
library(ggthemes)

load("splines/mtbqsDF.RData")
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

#### choose model ####

#get list of stan models
#original using default spline properties
modelDirectory <- "model-outputs/Odonata_stan_spline/v1"

#with 2/1 dimension on spline and k = c(5,8)
modelDirectory <- "model-outputs/Odonata_stan_spline/v2"

#with 2/1 dimension on spline and k = c(5,5)
modelDirectory <- "model-outputs/Odonata_stan_spline/v3"

#with 2/1 dimension on spline and k = c(10,10) and bs/cs
modelDirectory <- "model-outputs/Odonata_stan_spline/v4"

#with 2/1 dimension on spline and k = c(15,10) and bs/cs
modelDirectory <- "model-outputs/Odonata_stan_spline/v5"

#with 2/1 dimension on spline and k = c(8,5)
modelDirectory <- "model-outputs/Odonata_stan_spline/v6"

stanFiles <- list.files(modelDirectory) %>% str_subset("m_fit")

### model fits ####

readStanModel <- function(file, get="fits"){
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl(get,temp$Param))
  
  #get species name
  temp$File <- file
  temp$Species <- gsub("m_fit_summary_spacetime_","",temp$File)
  temp$Species <- gsub(".rds","",temp$Species)
  
  return(temp)
}

modelSummaries <- stanFiles %>%
  map_df(readStanModel) %>%
  mutate(siteIndex = parse_number(Param))

#get site and year information
fitDF <- readRDS("splines/fitDF.rds")
fitDF$siteIndex <- 1:nrow(fitDF)
modelSummaries <- left_join(modelSummaries,fitDF, by="siteIndex") 

#add coordinates
modelSummaries <- left_join(modelSummaries,mtbsDF, by="MTB")
nuMTBs <- length(unique(modelSummaries$MTB))

### national predictions ####

#function to apply to each file
readStanModel <- function(file, get="psi"){
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl(get,temp$Param))
  
  #get species name
  temp$File <- file
  temp$Species <- gsub("m_fit_summary_spacetime_","",temp$File)
  temp$Species <- gsub(".rds","",temp$Species)
  
  return(temp)
}

modelSummaries <- stanFiles %>%
  map_df(readStanModel) %>%
  mutate(siteIndex = parse_number(Param))

#get site and year information
siteInfo_NAs <- readRDS("splines/siteInfo_NAs.rds") %>%
  select(!c(Species,SpeciesOrig))
modelSummaries <- left_join(modelSummaries,siteInfo_NAs, by="siteIndex")

nuMTBs <- length(unique(siteInfo_NAs$MTB))

#### time-series ####

#analyse nationwide time series for each species

myspecies <- sort(unique(modelSummaries$Species))
annualTS <- modelSummaries %>%
            dplyr::group_by(Species,Year) %>%
            dplyr::summarise(total = sum(mean)/nuMTBs)

ggplot(annualTS)+
  geom_line(aes(x=Year, y=total))+
  facet_wrap(~Species)+
  theme_bw()

#very smooth time-series with v1... maybe try with more knots
#more wiggly with v2 - only ran for 20 species
#very smooth again with v3
#more wiggly with v4
#more wiggly with v5

#### spatial maps ####
allspecies <- sort(unique(modelSummaries$Species))

#plot all years for one species
selectspecies <- allspecies[69]
ggplot(filter(modelSummaries, Species==selectspecies))+
  geom_point(aes(x=x_MTB, y=y_MTB, colour=mean))+
  facet_wrap(~Year) +
  scale_color_viridis_c()

#plot spatial maps for start and end year

#just compare the start and the end for each species
modelSummaries_Limits <- filter(modelSummaries, Year %in% c(1990,2016))
ggplot(filter(modelSummaries_Limits, Species %in% allspecies[1:10]))+
  geom_point(aes(x=x_MTB, y=y_MTB, colour=mean))+
    facet_grid(Year~Species)+
  scale_color_viridis_c()

ggplot(filter(modelSummaries_Limits, Species %in% allspecies[11:20]))+
  geom_point(aes(x=x_MTB, y=y_MTB, colour=mean))+
  facet_grid(Year~Species)+
  scale_color_viridis_c()

ggplot(filter(modelSummaries_Limits, Species %in% allspecies[21:30]))+
  geom_point(aes(x=x_MTB, y=y_MTB, colour=mean))+
  facet_grid(Year~Species)+
  scale_color_viridis_c()

#v1 - not stripey
#v2
#v3 - very smooth - not stripey
#v4 - stripey for some
#v5 - stripey for some

#plot each year for each species
myYears <- 1991:2016

#mean
for(s in unique(modelSummaries$Species)){
  
  myMax <- max(modelSummaries$mean[modelSummaries$Species==s])
  
for(i in 1:length(myYears)){
  
    ggplot(subset(modelSummaries, Species == s & Year == myYears[i]))+
      geom_point(aes( x=x_MTB, y = y_MTB, colour = mean), size = 3)+
      scale_color_viridis_c("Occupancy",option = "A", direction = -1, limits=c(0,myMax))+
      theme_void()
    
ggsave(file=paste0("plots/species/spatial_maps/Map_",myYears[i],"_",s,".png"), width=5.5, height=6)    
       
  }
}

#and error maps
for(s in unique(modelSummaries$Species)){
  
  myMax <- max(modelSummaries$sd[modelSummaries$Species==s])
  myMin <- min(modelSummaries$sd[modelSummaries$Species==s])
  
  for(i in 1:length(myYears)){
    
    ggplot(subset(modelSummaries, Species == s & Year == myYears[i]))+
      geom_point(aes( x=x_MTB, y = y_MTB, colour = sd), size = 3)+
      scale_color_viridis_c("Occupancy_error",option = "A", direction = -1, limits=c(myMin,myMax))+
      theme_void()
    
    ggsave(file=paste0("plots/species/spatial_maps/error/Map_error_",myYears[i],"_",s,".png"), width=5.5, height=6)    
    
  }
}