#species richness patterns

setwd("~/winhome/dbowler/sMon-insects/traits")

#distribution data
distr <- read.csv("10750_2017_3495_MOESM2_ESM.csv",dec=",",as.is=T)
distr <- subset(distr, Country=="Germany")

#traits data
traits <- read.csv("10750_2017_3495_MOESM1_ESM.csv")
traits$species <- traits$Species
traits$Species <- paste(traits$Genus, traits$Species)
traits <- subset(traits, Species %in% distr$Species.name)

#change coordinate system to 
#UTM grid
library(mgrs)
library(plyr)
out <- ldply(distr$MGRS.WGS84, function(x) mgrs_to_utm(x))
distr<- data.frame(distr, out)

#majority are in zone 32
table(distr$zone)

#convert all into zone 31 to 32
distr31 <- subset(distr, zone == 31)
coordinates(distr31)<-c("easting","northing")
proj4string(distr31) <- CRS("+proj=utm +zone=31 +datum=WGS84")
distr31 <- spTransform(distr31,CRS("+proj=utm +zone=32 +datum=WGS84"))

#convert all into zone 33 to 32
distr33 <- subset(distr, zone == 33)
coordinates(distr33)<-c("easting","northing")
proj4string(distr33) <- CRS("+proj=utm +zone=33 +datum=WGS84")
distr33 <- spTransform(distr33,CRS("+proj=utm +zone=32 +datum=WGS84"))

#add to original dataset
distr$easting[distr$zone==31] <- distr31@coords[,1]
distr$northing[distr$zone==31] <- distr31@coords[,2]
distr$easting[distr$zone==33] <- distr33@coords[,1]
distr$northing[distr$zone==33] <- distr33@coords[,2]

#round coordinates
distr$easting <- signif(distr$easting, digits = 3)
distr$northing <- signif(distr$northing, digits = 3)

#overlay with the german states...
setwd("~/winhome/dbowler/sMon-insects/Spatial_data")
library(sp)
germany <- readRDS("gadm36_DEU_1_sp.rds")
germany <- spTransform(germany,"+proj=utm +zone=32 +datum=WGS84")
germanyF <- fortify(germany)
ggplot() + geom_path(data=germanyF, aes(long, lat, group = group))
             
#plotting
library(ggplot2)

#total species
distrSummary1 <- ddply(distr,.(easting,northing),summarise,nuSpecies=length(unique(Species.name)))
g1 <- ggplot(subset(distrSummary1,nuSpecies>30))+
  geom_point(aes(x=easting, y=northing, colour=nuSpecies),shape=15,size=5)+
  scale_colour_gradient2(low="steelblue",mid="white",high="red",midpoint=mean(distrSummary1$nuSpecies))+
  geom_path(data=germanyF, aes(long, lat, group = group))+
  ggtitle("total species")

summary(distrSummary1$nuSpecies)
  
#number of red list
table(traits$EU.Red.List.Europe)
#Critically Endangered        Data Deficient            Endangered         Least Concern       Near Threatened        Not Applicable 
#0                     0                     0                    69                     9                     0 
#Not Evaluated           Not present            Vulnerable 
#0                     0                     3 

table(traits$Habitat.Directive)
#No Yes 
#70  11

#total species
distrSummary2 <- ddply(subset(distr, Species.name %in% traits$Species[traits$Habitat.Directive=="Yes"]),.(easting,northing),summarise,nuSpecies=length(unique(Species.name)))
g2<- ggplot(distrSummary2)+
  geom_point(aes(x=easting, y=northing, colour=nuSpecies),shape=15,size=5)+
  scale_colour_gradient2(low="steelblue",mid="white",high="red",midpoint=mean(distrSummary2$nuSpecies))+
  geom_path(data=germanyF, aes(long, lat, group = group))+
  ggtitle("habitat directive species")

summary(distrSummary2$nuSpecies)

#combine
library(cowplot)
plot_grid(g1,g2)




