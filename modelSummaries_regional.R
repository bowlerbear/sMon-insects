#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyverse)
library(sf)
library(tmap)
library(boot)

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

plotTS_regional <- function(x){
  require(ggplot2)
  g1 <- ggplot(data=x,aes(x=Year,y=mean,group=factor(Site)))+
    geom_line(aes(x=Year,y=mean,colour=factor(Site)))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.,fill=factor(Site)),alpha=0.5)+
    facet_wrap(~Species)
  print(g1)
}

###Species list#############################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

allSpecies <- read.delim("specieslist_odonata.txt",as.is=T)$Species

### get ecoregions ########################################

nr <- st_read("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/Narurraeume/naturraeume_polygon.shp")

plot(nr)

#simplify coarse naturraum column
nr$CoarseNaturraum <- sapply(as.character(nr$FolderPath),function(x)strsplit(x,"/")[[1]][2])

#dissolve 
nr_dissolved <-
  nr %>%
  group_by(CoarseNaturraum) %>%
  summarise()

#plot
tm_shape(nr_dissolved)+
  tm_polygons("CoarseNaturraum")


#1:7 are:
CoarseRegions <- c("Alpen","Alpenvorland","Nordostdeutsches Tiefland","Nordwestdeutsches Tiefland", 
                   "Oestliches Mittelgebirge","Suedwestdeutsches Mittelgebirge ",
                   "Westliches Mittelgebirge")

### ecoregion 1 (coarse) ##################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_sparta_regional/6869992"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#missing these species:
#[1] "Calopteryx splendens"     "Coenagrion hastulatum"    "Erythromma viridulum"     "Ischnura pumilio"    
#[5] "Libellula fulva"          "Libellula quadrimaculata" "Orthetrum albistylum"     "Platycnemis pennipes"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_sparta_regional_nation_naturraum_adult_")

#annual time series - psi.fs
annualDF <- getBUGSfitsII(modelDF,param="psi.fs")
annualDF$Year <- annualDF$Year + 1979
plotTS_regional(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#1281 15960

#plot just for two species
ggplot(data=subset(annualDF,Species=="Anax imperator"))+
      geom_line(aes(x=Year,y=mean,))+
      geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
      facet_wrap(~Site)
#same for all of them...

ggplot(data=subset(annualDF,Species=="Crocothemis erythraea"))+
  geom_line(aes(x=Year,y=mean,))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Site)
#same for all of them

#plot cr.a instead
annualDF <- getBUGSfitsII(modelDF,param="cr.a")
annualDF$Year <- annualDF$Year + 1979

#run same plots as above

#better!!
#need to be back transformed...

#do we have 7 regions for all species??
annualDF %>%
  group_by(Species) %>%
  summarise(nuSites = n_distinct(Site))
#all 7

#get mean per ecoregion:

ecoregionMeans <- annualDF %>%
  group_by(Site) %>%
  summarise(meanOcc  = mean(inv.logit(mean)))
ecoregionMeans$CoarseNaturraum <- CoarseRegions

#plot it
nr_dissolved  <- left_join(nr_dissolved,ecoregionMeans,by="CoarseNaturraum")

tm_shape(nr_dissolved)+
  tm_polygons("meanOcc")


#mean change per lab
ecoregionMeans <- annualDF %>%
  group_by(Site) %>%
  summarise(occChange  = mean(inv.logit(mean[Year==37]) - inv.logit(mean[Year==1])))
ecoregionMeans$CoarseNaturraum <- CoarseRegions

nr_dissolved  <- left_join(nr_dissolved,ecoregionMeans,by="CoarseNaturraum")

tm_shape(nr_dissolved)+
  tm_polygons("occChange")

#change for one species
ecoregionMeans <- annualDF %>%
  filter (Species=="Anax imperator") %>%
  group_by(Site) %>%
  summarise(spoccChange  = inv.logit(mean[Year==37]) - inv.logit(mean[Year==1]))
ecoregionMeans$CoarseNaturraum <- CoarseRegions

nr_dissolved  <- left_join(nr_dissolved,ecoregionMeans,by="CoarseNaturraum")

tm_shape(nr_dissolved)+
  tm_polygons("spoccChange")

#### ecoregion 3 (finest) #################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_regional_nation_raumunit_sparta/7180407"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_sparta_regional_nation_raumunit_adult_")

#annual time series - psi.fs
annual_z_DF <- getBUGSfitsII(modelDF,param="muZ")

#plot time-series
ggplot(annual_z_DF)+
  geom_line(aes(x=Year,y=mean,colour=factor(Site)))+
  facet_wrap(~Code)+
  theme(legend.position = "none")
#there is change over time...????

#match with raums

#get siteInfo from HPC script
annual_z_DF  <- merge(annual_z_DF,siteInfo,by.x="Site",by.y="siteIndex")
#fix German accents
annual_z_DF$Name <- sapply(annual_z_DF$Natur,function(x)iconv(x,from="UTF-8",to="latin1"))

#what species do we have data for
unique(annual_z_DF$Code)

#pull out data for 1 species
annual_z_DF_species1 <- subset(annual_z_DF,Code=="Aes_vir")

#map to naturraum polygons
nr_sf <- st_as_sf(nr)
nr_sf <- full_join(nr_sf,annual_z_DF_species1,by="Name")
nr_sf <- subset(nr_sf, !is.na(mean))

tm_shape(nr_sf)+
    tm_polygons("mean")+
    tm_facets(by="Year")

#looks like very little change over time...

### ecoregion 2 (mid) ##################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_regional_midnaturraumtrend_sparta/7474533"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_sparta_regional_nation_midnaturraumtrends_adult_")

#what parameters do we have?
unique(modelDF$Param)

#annual time series - psi.fs
annualtrendsDF <- getBUGSfits(modelDF,param="mraumtrend")

#find out which indices correspond to which
siteInfo <- readRDS("siteInfo_midnaturraum.rds")
annualtrendsDF$MidNaturraum <- siteInfo$MidNaturraum[
  match(annualtrendsDF$ParamNu,siteInfo$mnIndex)]

#get % of increasing species per naturraum
increasingProp <- annualtrendsDF %>%
                    group_by(MidNaturraum) %>%
                    summarise(incr_prop = mean(mean>0))

#add also to the nr shapefile
load("mtbqsDF.RData")
nr$MidNaturraum <- gdata::trim(mtbqsDF$MTB_MidNatur[match(nr$Name,mtbqsDF$MTB_Natur)])

#dissolve to midnaturraum
nr_dissolved <-
  nr %>%
  group_by(MidNaturraum) %>%
  summarise()
plot(nr_dissolved)

#merge summary data and shapefile
all(increasingProp$MidNaturraum %in% nr$MidNaturraum)#TRUE
nr_dissolved  <- left_join(nr_dissolved,increasingProp,by="MidNaturraum")

tm_shape(nr_dissolved)+
  tm_polygons("incr_prop")

#no obvious pattern
#no data from some midnaturraums

#plot for a specific species
speciesChange <- subset(annualtrendsDF,Species=="Crocothemis erythraea")
nr_dissolved_sp  <- left_join(nr_dissolved,speciesChange,by="MidNaturraum")

tm_shape(nr_dissolved_sp)+
  tm_polygons("mean")

speciesChange <- subset(annualtrendsDF,Species=="Sympetrum danae")
nr_dissolved_sp  <- left_join(nr_dissolved,speciesChange,by="MidNaturraum")

tm_shape(nr_dissolved_sp)+
  tm_polygons("mean")

### spline model ##########################################



###end######################################################
