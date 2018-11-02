#formatting traits
library(plyr)

###########################################################################################

#get total specieslist
species <- read.delim("traits/specieslist_odonata.txt",as.is=T)

###########################################################################################

#get new atlas traits data data
traits1<- read.csv("traits/10750_2017_3495_MOESM1_ESM.csv",as.is=T)
traits1$species <- traits1$Species
traits1$Species <- paste(traits1$Genus,traits1$Species)

#get new atlas distribution data
atlas <- read.csv("traits/10750_2017_3495_MOESM2_ESM.csv",as.is=T)
#make sure first letter is a capital
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
atlas$Species.name <- sapply(atlas$Species.name,function(x)firstup(x))
atlas$Species.name[which(atlas$Species.name=="Leucorrhinia Dubia")]<-"Leucorrhinia dubia"

atlas <- ddply(atlas,.(Species.name),summarise,
               lowerLat=min(Latitude),
               upperLat=max(Latitude),
               minLong = min(Longitude),
               maxLong = max(Longitude),
               nuEuroGrids = length(MGRS.WGS84))

#check names match up
traits1$Species[!traits1$Species %in% atlas$Species.name]
#all there
atlas$Species.name[!atlas$Species.name %in% traits1$Species]
traits1 <- merge(traits1,atlas,by.x="Species",by.y="Species.name")
all(species$Species %in% traits1$Species) 

#get my trait compilation based off the atlas species list
traits2 <- read.csv("traits/odonata_Germany_list.csv",as.is=T)
traits2 <- merge(traits2,traits1,by="Species")
nrow(traits2)#should be 81

#####################################################################################

#get oliver schweigers

traits3 <- read.csv("traits/from Oliver Schweiger/dragonflies_STI_larval development time_withNotes.csv",as.is=T)
traits3$species %in% species$Species
traits3$species[!traits3$species %in% species$Species]
traits3$species[which(traits3$species=="Lestes viridis")]<-"Chalcolestes viridis"
all(traits3$species %in% species$Species)
traits2 <- merge(traits2,traits3,by.x="Species",by.y="species",all.x=T)

#####################################################################################

#get christian hofs

#species temperature preferences
traits4 <- read.delim("traits/from Christian Hof/Dragonfly_temperature_niche_6.txt",as.is=T)
traits4$Species <- sapply(traits4$Species,function(x)strsplit(x,"\\.")[[1]][1])

#merge with other data
traits6 <- read.csv("traits/from Christian Hof/dragons_for_Diana_habitats.csv",as.is=T)
traits6$all_spec <- sapply(traits6$all_spec,function(x)strsplit(x,"\\.")[[1]][1])
traits4 <- merge(traits4,traits6,by.x=("Species"),by.y="all_spec",all=T)

#fill in missing names
traits4$Species[traits4$species_name==""]
traits4$species_name[traits4$species_name==""]<-
  c("Coenagrion puella","Coenagrion pulchellum","Cordulegaster heros",
    "Cordulegaster insignis","Cordulegaster picta","Enallagma cyathigerum",
    "Gomphus flavipes","Gomphus schneiderii","Anax ephippiger",
    "Ischnura elegans","Ischnura genei","Lestes Macrostigma","Lestes viridis",
    "Onychogomphus costae","Orthetrum chrysostigma","Pantala flavescens","Platycnemis pennipes",
    "Pyrrhosoma nymphula","Somatochlora sahlbergi","Zygonyx torridus")

#get other traits data from CH
traits5 <- read.csv("traits/from Christian Hof/dragons_for_Diana.csv",as.is=T)
traits5$species <- gsub("_"," ",traits5$species)
traits4 <- merge(traits4, traits5, by.x="species_name",by.y="species",all=T)
traits4 <- traits4[,-2]

traits2$Species[!traits2$Species %in% traits4$species_name]
#"Chalcolestes viridis" "Lestes macrostigma"

#rename christian's names to our
traits4$species_name[which(traits4$species_name=="Lestes viridis")]<-"Chalcolestes viridis"
traits4$species_name[which(traits4$species_name=="Lestes Macrostigma")]<-"Lestes macrostigma"

#combine
alltraits <- merge(traits2,traits4,by.x="Species",by.y="species_name",all.x=T)

#####################################################################################

#get new body size and wing length information
bodysize <- read.csv("traits/from Christian Hof/DragonflyTraits_bodysize.csv")
bodysize$medAb <- (bodysize$Ab_Lower + bodysize$Ab_Upper)/2
bodysize$medHw <- (bodysize$Hw_Lower + bodysize$Hw_Upper)/2
species$Species [!species$Species %in% bodysize$Species]

#"Aeshna caerulea"       "Boyeria irene"         "Coenagrion hylas"     
#[4] "Coenagrion scitulum"   "Epitheca bimaculata"   "Gomphus simillimus"   
#[7] "Lestes macrostigma"    "Nehalennia speciosa"   "Onychogomphus uncatus"
#[10] "Orthetrum albistylum"  "Oxygastra curtisii" 

alltraits$medAb <- bodysize$medAb[match(alltraits$Species,bodysize$Species)] 
alltraits$medHw <- bodysize$medHw[match(alltraits$Species,bodysize$Species)] 

save(alltraits, file="alltraits.RData")

######################################################################################

#check out missing data
library(reshape2)
alltraits <- alltraits[,c(1,5:6,28:48)]
alltraitsM <- melt(alltraits,id=c("Species","total","since"))
alltraitsM <- subset(alltraitsM, is.na(value))
alltraitsM
#not too bad
  
######################################################################################

#look at histograms
library(ggplot2)
ggplot(data=alltraits,aes(x=total))+
  geom_histogram()+
  theme_bw()+
  xlab("Range size")+
  ylab("Number of species")


ggplot(data=alltraits,aes(x=Habitat.y))+
  geom_histogram(stat="count")+
  theme_bw()+
  xlab("Habitat preference")+
  ylab("Number of species")

