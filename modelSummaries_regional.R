#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyverse)
library(sf)
library(tmap)
library(boot)
library(wesanderson)

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

###Species list#############################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

allSpecies <- read.delim("specieslist_odonata.txt",as.is=T)$Species

### get ecoregions ########################################

nr <- st_read("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/Narurraeume/naturraeume_polygon.shp")

plot(nr)

#simplify coarse naturraum column
nr$CoarseNaturraum <- gdata::trim(sapply(as.character(nr$FolderPath),function(x)strsplit(x,"/")[[1]][2]))

#dissolve 
nr_dissolved <-
  nr %>%
  group_by(CoarseNaturraum) %>%
  summarise()

#plot
tm_shape(nr_dissolved)+
  tm_polygons("CoarseNaturraum")


#1:7 are:
nr_dissolved$CoarseNaturraum
CoarseRegions <- c("Alpen","Alpenvorland","Nordostdeutsches Tiefland","Nordwestdeutsches Tiefland", 
                   "Oestliches Mittelgebirge","Suedwestdeutsches Mittelgebirge",
                   "Westliches Mittelgebirge")

### ecoregion 1 (coarse) ##################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_regional_naturraum/76676"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModels(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_sparta_regional_nation_naturraum_adult_")

#annual time series for each coarse naturruam
#pull out muZ.cra
annualDF <- getBUGSfitsII(modelDF,param="muZ.cra")
annualDF$Year <- annualDF$Year + 1989
plotTS_regional(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#633 14458

#match region indices to actual names
coarseRaums <- data.frame(Index=1:7,
                          Names=c("Alpen",
                 "Alpenvorland",
                 "Nordostdeutsches Tiefland",
                 "Nordwestdeutsches Tiefland",
                 "Oestliches Mittelgebirge",
                 "Suedwestdeutsches Mittelgebirge",
                 "Westliches Mittelgebirge"))
annualDF$Naturraum <- coarseRaums$Names[match(annualDF$Site,coarseRaums$Index)]
plotTS_regional(annualDF)

#plot just for two species
ggplot(data=subset(annualDF,Species=="Anax imperator"))+
      geom_line(aes(x=Year,y=mean,))+
      geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
      facet_wrap(~Naturraum)

ggplot(data=subset(annualDF,Species=="Crocothemis erythraea"))+
  geom_line(aes(x=Year,y=mean,))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Naturraum)

#do we have 7 regions for all species??
annualDF %>%
  group_by(Species) %>%
  summarise(nuSites = n_distinct(Site))
#all 7

#density ridge plots
library(ggridges)
theme_set(theme_ridges())

#change between first and last year
occChange <- annualDF %>%
  group_by(Species,Naturraum) %>%
  summarise(occChange  = median(mean[Year==max(annualDF$Year)]/mean[Year==min(annualDF$Year)]))

#extreme value
summary(occChange$occChange)
occChange$occChange[occChange$occChange>2] <- 2

ggplot(occChange, aes(x = occChange, y = Naturraum)) +
  geom_density_ridges(aes(fill = Naturraum)) +
  theme(legend.position = "none")+
  geom_vline(xintercept=1, linetype="dashed")
#two humps???


#trends - see if there is a regres.psi in the modelDF??
occTrends <- annualDF %>%
  group_by(Species,Naturraum) %>%
  do(model = broom::tidy(lm(mean ~ Year, data = .))) %>% 
  unnest(model) %>%
  filter(term=="Year") %>%
  rename(trend="estimate",trend_se="std.error")

#order naturaum by proportion of trends bigger than zero

ggplot(occTrends, aes(x = trend, y = reorder(Naturraum,desc(Naturraum)))) +
  geom_density_ridges(aes(fill = Naturraum)) +
  theme(legend.position = "none")+
  geom_vline(xintercept=0, linetype="dashed")+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 7,type = "continuous"))+
  ylab("")

### maps ####

#plot mean per ecoregion:
ecoregionMeans <- annualDF %>%
  group_by(Naturraum) %>%
  summarise(meanOcc  = sum(mean)) %>%
  rename(CoarseNaturraum=Naturraum)
nr_dissolved  <- left_join(nr_dissolved,ecoregionMeans,by="CoarseNaturraum")
tm_shape(nr_dissolved)+tm_polygons("meanOcc")

#mean change per ecoregion
ecoregionMeans <- annualDF %>%
  group_by(Naturraum) %>%
  summarise(occChange  = median(mean[Year==2017]/mean[Year==1990]))%>%
  rename(CoarseNaturraum=Naturraum)
nr_dissolved  <- left_join(nr_dissolved,ecoregionMeans,by="CoarseNaturraum")
tm_shape(nr_dissolved)+tm_polygons("occChange")

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
speciesChange <- subset(annualtrendsDF,Species=="Anax imperator")
nr_dissolved_sp  <- left_join(nr_dissolved,speciesChange,by="MidNaturraum")

tm_shape(nr_dissolved_sp)+
  tm_polygons("mean")

speciesChange <- subset(annualtrendsDF,Species=="Coenagrion hastulatum")
nr_dissolved_sp  <- left_join(nr_dissolved,speciesChange,by="MidNaturraum")

tm_shape(nr_dissolved_sp)+
  tm_polygons("mean")


#### state model ##########################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_state/992948"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModels(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="outSummary_sparta_regional_nation_state_adult_")

#what parameters do we have?
unique(modelDF$Param)

#annual time series - psi.fs
annualtrendsDF <- getBUGSfitsII(modelDF,param="muZ.state")

#relabels sites as states
myStates <- data.frame(ID=1:13,
                       State=c("Baden-Württemberg","Bayern","Brandenburg","Hessen","Mecklenburg-Vorpommern", "Niedersachsen","Nordrhein-Westfalen","Rheinland-Pfalz","Saarland","Sachsen","Sachsen-Anhalt", "Schleswig-Holstein", "Thüringen")) 
annualtrendsDF$State <- myStates$State[match(annualtrendsDF$Site,myStates$ID)]

#format year
annualtrendsDF$Year <- annualtrendsDF$Year+1989

#plot data
ggplot(annualtrendsDF)+
  geom_line(aes(x=Year, y=mean, colour=State))+
  facet_wrap(~Species)
  
#long-term trends
annualtrendsDF <- getBUGSfits(modelDF,param="regres.psi")
annualtrendsDF$State <- myStates$State[match(annualtrendsDF$ParamNu,myStates$ID)]

#ridges
library(ggridges)
theme_set(theme_ridges())
ggplot(annualtrendsDF, aes(x = mean, y = State)) +
  geom_density_ridges(aes(fill = State)) +
  theme(legend.position = "none")+
  geom_vline(xintercept=0, linetype="dashed")

#order by proportion of species with negative trends
annualTrendsDF_Prop <- annualtrendsDF %>%
  dplyr::group_by(State) %>%
  dplyr::summarise(prop = mean(mean<0))%>%
  arrange(.,desc(prop))
annualtrendsDF$State <- factor(annualtrendsDF$State, 
                               levels=annualTrendsDF_Prop$State)

ggplot(annualtrendsDF, aes(x = mean, y = State)) +
  geom_density_ridges(aes(fill = State)) +
  theme(legend.position = "none")+
  geom_vline(xintercept=0, linetype="dashed")

n### spline model ##########################################


###end######################################################
