library(tidyverse)
library(ggthemes)

load("splines/mtbqsDF.RData")
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

#### choose model directory ####

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

#spatial only and k = 10
modelDirectory <- "model-outputs/Odonata_stan_spline/v7"

#spatial-temporal and k = 8
modelDirectory <- "model-outputs/Odonata_stan_spline/v8"

#spatial-temporal and k = 12
modelDirectory <- "model-outputs/Odonata_stan_spline/v9"

#spatial-temporal and k = 8 + expanded grid
modelDirectory <- "model-outputs/Odonata_stan_spline/v10"

#spatial-temporal and k = 7 + expanded grid
modelDirectory <- "model-outputs/Odonata_stan_spline/v11"

### get list of models ####

stanFiles <- list.files(modelDirectory) %>% str_subset("m_fit")

### model fits ####

# readStanModel <- function(file, get="fits"){
#   
#   temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
#   
#   #get rows for the parameter of interest
#   temp$Param <- row.names(temp)
#   temp <- subset(temp, grepl(get,temp$Param))
#   
#   #get species name
#   temp$File <- file
#   temp$Species <- gsub("m_fit_summary_spacetime_","",temp$File)
#   temp$Species <- gsub(".rds","",temp$Species)
#   
#   return(temp)
# }
# 
# modelSummaries <- stanFiles %>%
#   map_df(readStanModel) %>%
#   mutate(siteIndex = parse_number(Param))
# 
# #get site and year information
# fitDF <- readRDS("splines/fitDF.rds")
# fitDF$siteIndex <- 1:nrow(fitDF)
# modelSummaries <- left_join(modelSummaries,fitDF, by="siteIndex") 
# 
# #add coordinates
# modelSummaries <- left_join(modelSummaries,mtbsDF, by="MTB")
# nuMTBs <- length(unique(modelSummaries$MTB))

### national predictions ####

#function to apply to each file
readStanModel <- function(file, get="psi"){
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl(get,temp$Param))
  temp <- subset(temp, !grepl("beta_psi",temp$Param))
                 
  #get species name
  temp$File <- file
  temp$Species <- gsub("m_fit_summary_spacetime_","",temp$File)
  temp$Species <- gsub(".rds","",temp$Species)
  
  return(temp)
}

modelSummaries <- stanFiles %>%
  map_df(readStanModel) %>%
  mutate(siteIndex = parse_number(Param))

#get site and year information - using only MTBQ as the siteInfo
#siteInfo_NAs <- readRDS("splines/siteInfo_NAs.rds") %>%
#  select(!c(Species,SpeciesOrig))
siteInfo_NAs <- readRDS("splines/siteInfo_NAs.rds") %>% #using the MTBQ data frame
    dplyr::select(!c(Species,SpeciesOrig)) %>%
    dplyr::filter(type!="extension")

#merge
modelSummaries <- inner_join(modelSummaries,siteInfo_NAs, by="siteIndex")
nuMTBs <- length(unique(siteInfo_NAs$MTB))
#modelSummaries$x_MTB <- modelSummaries$x
#modelSummaries$y_MTB <- modelSummaries$y

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

#v1 - very smooth time-series with v1... maybe try with more knots
#v2 - more wiggly with v2 - only ran for 20 species
#v3 - very smooth again with v3
#v4 - more wiggly with v4
#v5 - more wiggly with v5
#v6 - very smooth
#v7 - constant - spatial only model
#v8 - looks good
#v9 - wiggly and look good
#v10 - looks good
#v11 - smooth but ok

#### spatial maps ####

allspecies <- sort(unique(modelSummaries$Species))

#plot all years for one species
selectspecies <- allspecies[11]

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
#v2 - some stripey
#v3 - very smooth - stripey
#v4 - stripey for some
#v5 - stripey for some
#v6 - stripey
#v7 - spatial only

#plot each year for each species
myYears <- 1990:2016

#mean
for(s in unique(modelSummaries$Species)){
  
  myMax <- max(modelSummaries$mean[modelSummaries$Species==s])
  
for(i in 1:length(myYears)){
  
    ggplot(subset(modelSummaries, Species == s & Year == myYears[i]))+
      geom_point(aes( x=x_MTB, y = y_MTB, colour = mean), size = 3)+
      scale_color_viridis_c("Occupancy",option = "A", direction = -1, limits=c(0,myMax))+
      theme_void()
    
ggsave(file=paste0("plots/species/spatial_maps/v11/Map_",s,"_",myYears[i],".png"), width=5.5, height=6)    
       
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