#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyverse)
library(sf)
library(tmap)
library(boot)

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

plotTS_regional <- function(x){
  require(ggplot2)
  g1 <- ggplot(data=x,aes(x=Year,y=mean,group=factor(Site)))+
    geom_line(aes(x=Year,y=mean,colour=factor(Site)))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.,fill=factor(Site)),alpha=0.5)+
    facet_wrap(~Species)
  print(g1)
}

###Species list#################################################################################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

allSpecies <- read.delim("specieslist_odonata.txt",as.is=T)$Species

### get ecoregions #########################################################################################

nr <- st_read("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/Narurraeume/naturraeume_polygon.shp")

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

###Model summaries#########################################################################################

### ecoregion 1 ############################################################################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_sparta_regional/6869992"

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

#### ecoregion 2 model ###################################################################



### spline model #########################################################################



