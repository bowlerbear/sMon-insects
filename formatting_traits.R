#formatting traits
library(plyr)
library(ggplot2)

###########################################################################################

#get total specieslist
species <- read.delim("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/specieslist_odonata.txt",as.is=T)

###########################################################################################

#get new atlas traits data data
traits1<- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/10750_2017_3495_MOESM1_ESM.csv",as.is=T)
traits1$species <- traits1$Species
traits1$Species <- paste(traits1$Genus,traits1$Species)

#get new atlas distribution data
atlas <- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/10750_2017_3495_MOESM2_ESM.csv",as.is=T)
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
traits2 <- read.delim("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/specieslist_Germany_atlastraits.txt",as.is=T)
traits2 <- merge(traits2,traits1,by="Species")
nrow(traits2)#should be 81
head(traits2)

# names(traits2)
# [1] "Species"                             "German_name"                        
# [3] "Family"                              "Suborder"                           
# [5] "total"                               "since"                              
# [7] "Europe_trend"                        "Flight_start"                       
# [9] "Flight_end"                          "Water"                              
# [11] "Habitat.x"                           "Comment"                            
# [13] "Taxon"                               "Genus"                              
# [15] "Habitat.y"                           "Endemic.to.Europe"                  
# [17] "EU.Red.List.Europe"                  "Habitat.Directive"                  
# [19] "Threath.status.on.European.Red.List" "species"                            
# [21] "lowerLat"                            "upperLat"                           
# [23] "minLong"                             "maxLong"                            
# [25] "nuEuroGrids"

#####################################################################################

#get oliver schweigers - development tinme, voltinism, migration

traits3 <- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/from Oliver Schweiger/dragonflies_STI_larval development time_withNotes.csv",as.is=T)
traits3$species %in% species$Species
traits3$species[!traits3$species %in% species$Species]
traits3$species[which(traits3$species=="Lestes viridis")]<-"Chalcolestes viridis"

#all names correct?
all(traits3$species %in% species$Species)

#all species there?
all(traits2$Species %in% traits3$species)
traits2$Species[!traits2$Species %in% traits3$species]
#[1] "Aeshna caerulea"       "Aeshna viridis"        "Boyeria irene"        
#[4] "Ceriagrion tenellum"   "Coenagrion hylas"      "Coenagrion scitulum"  
#[7] "Gomphus pulchellus"    "Gomphus simillimus"    "Lestes macrostigma"   
#[10] "Onychogomphus uncatus" "Orthetrum albistylum"  "Oxygastra curtisii"

traits2 <- merge(traits2,traits3,by.x="Species",by.y="species",all.x=T)

#####################################################################################

#get christian hofs - not needed anymore

#species temperature preferences
# traits4 <- read.delim("traits/from Christian Hof/Dragonfly_temperature_niche_6.txt",as.is=T)
# traits4$Species <- sapply(traits4$Species,function(x)strsplit(x,"\\.")[[1]][1])
# 
# #merge with other data
# traits6 <- read.csv("traits/from Christian Hof/dragons_for_Diana_habitats.csv",as.is=T)
# traits6$all_spec <- sapply(traits6$all_spec,function(x)strsplit(x,"\\.")[[1]][1])
# traits4 <- merge(traits4,traits6,by.x=("Species"),by.y="all_spec",all=T)
# 
# #fill in missing names
# traits4$Species[traits4$species_name==""]
# traits4$species_name[traits4$species_name==""]<-
#   c("Coenagrion puella","Coenagrion pulchellum","Cordulegaster heros",
#     "Cordulegaster insignis","Cordulegaster picta","Enallagma cyathigerum",
#     "Gomphus flavipes","Gomphus schneiderii","Anax ephippiger",
#     "Ischnura elegans","Ischnura genei","Lestes Macrostigma","Lestes viridis",
#     "Onychogomphus costae","Orthetrum chrysostigma","Pantala flavescens","Platycnemis pennipes",
#     "Pyrrhosoma nymphula","Somatochlora sahlbergi","Zygonyx torridus")
# 
# #get other traits data from CH
# traits5 <- read.csv("traits/from Christian Hof/dragons_for_Diana.csv",as.is=T)
# traits5$species <- gsub("_"," ",traits5$species)
# traits4 <- merge(traits4, traits5, by.x="species_name",by.y="species",all=T)
# traits4 <- traits4[,-2]
# 
# traits2$Species[!traits2$Species %in% traits4$species_name]
# #"Chalcolestes viridis" "Lestes macrostigma"
# 
# #rename christian's names to our
# traits4$species_name[which(traits4$species_name=="Lestes viridis")]<-"Chalcolestes viridis"
# traits4$species_name[which(traits4$species_name=="Lestes Macrostigma")]<-"Lestes macrostigma"
# 
# #combine
# alltraits <- merge(traits2,traits4,by.x="Species",by.y="species_name",all.x=T)

####################################################################################

#get new temp preference data based on new atlas

tempPref <- read.delim("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/speciesNiches_germanDragonflies.txt")
all(tempPref$Species.name %in% traits2$Species) 

alltraits <- merge(traits2,tempPref,by.x="Species",by.y="Species.name",all=T)

#####################################################################################

#get new body size and wing length information

bodysize <- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/sizes/specieslist_odonata_size.csv")
bodysize$medAb <- (bodysize$Ab.lower.in.mm + bodysize$Ab.upper.in.mm)/2
bodysize$medHw <- (bodysize$Hindwing.lower.in.mm + bodysize$Hindwing.upper.in.mm)/2
species$Species [!species$Species %in% bodysize$Species]

#add to my data
alltraits$medAb <- bodysize$medAb[match(alltraits$Species,bodysize$Species)] 
alltraits$medHw <- bodysize$medHw[match(alltraits$Species,bodysize$Species)] 
alltraits$HabitatBreadth <- bodysize$HabitatSpecialist[match(alltraits$Species,bodysize$Species)] 

######################################################################################

#traits from Frank Suhling:
library(gdata)
setwd("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/from_Frank_Suhling")
ftraits <- read.xls("Traits_German_Odonata.xls",skip=2)
ftraits <- ftraits[,c(4,10:21)]
ftraits <- ftraits[,c(1:6,10:13)]
names(ftraits)[1] <- "Species"
species$Species [!species$Species %in% ftraits$Species]
#"Gomphus flavipes"   "Libellula depressa" "Libellula fulva"
#check for synonyms

#weighted mean of voltinism
ftraits$VoltinismProp <- apply(ftraits,1,function(x)(t(as.numeric(x[2:6])) %*% 5:1)/10)
#check meaning of columns

#overwintering strategy
subset(ftraits,Larva>5)
ftraits$Overwintering <- NA
ftraits$Overwintering[ftraits$Larva>5] <- "Larva"
ftraits$Overwintering[ftraits$Adult>5|ftraits$Not>5] <- "Other"
ftraits$Overwintering[ftraits$Egg>5] <- "Egg"
table(ftraits$Overwintering)

alltraits <- merge(alltraits,ftraits[,c("Species","VoltinismProp","Overwintering")],by="Species",all.x=T)

#######################################################################################

#from Domnik
setwd("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/from_Dominik_Jab")
zeigerwerte <- read.xls("Zeigerwerte_final.xls")

######################################################################################

#data from freshwater ecology
fwecol <- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/from FWecolol/tvtexport_data.csv")
#species name conversion
fwecolSpecies<- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/from FWecolol/specieslist_odonata.csv")
fwecol$Species <- fwecolSpecies$Species[match(fwecol$Fwecol_Species,fwecolSpecies$Fwecol_Species)]

#get mean szp for each species
szp <- fwecol[,grepl("szp",names(fwecol))]
szp[is.na(szp)] <- 0
SZP <- apply(szp,1,function(x)(as.numeric(x) %*% 1:9)/sum(x))
SZP[is.na(SZP)]<- NA

#get mean saprotic index
saprobit_Austria <- fwecol$saprobity...Austria_si

#also for the netherlands

#xeno - non existent
#oligo - normal
#beta -moderate
#alpha - strong
#poly - present locally
saprobit_Netherlands <- fwecol[,grepl("saprobity...Netherlands",names(fwecol))]
saprobit_Netherlands[is.na(saprobit_Netherlands)] <- 0
saprobit_Netherlands <- apply(saprobit_Netherlands,1,function(x)(as.numeric(x) %*% 1:5)/sum(x))
saprobit_Netherlands[is.na(saprobit_Netherlands)]<- NA

#combine all
fwecol <- data.frame(Species=fwecol$Species,
                     SZP=SZP,
                     saprobit_Austria=saprobit_Austria,
                     saprobit_Netherlands=saprobit_Netherlands)
pairs(fwecol[,2:4])
#higher the zone preference, the higher the saprobit index
alltraits <- merge(alltraits,fwecol,by="Species")

#####################################################################################

#rename some column headings
names(alltraits)[which(names(alltraits)=="total")] <- "germanRange"
names(alltraits)[which(names(alltraits)=="Habitat.y")] <- "Habitat"
names(alltraits)[which(names(alltraits)=="development.time..mean.")] <- "devTime"

######################################################################################

save(alltraits, file="alltraits.RData")

######################################################################################

#investigating correlations among traits

vars <- c("germanRange","Flight_start","Habitat","nuEuroGrids",
          "devTime","meanTemp","sdTemp","medHw","VoltinismProp")

alltraitsV <- alltraits[,vars]
unique(alltraitsV$Habitat)
alltraitsV$Habitat <- ifelse(alltraitsV$Habitat=="Standing",0,1)
str(alltraitsV)
names(alltraitsV) <- c("german range","start flight","running water","europe range",
                       "dev time","temp pref","temp range","wing length","voltinism")

#options
#hetcor() offers me the discrimination into polyserial and polychoric correlations, but no p-values.
#corrplot

library(corrplot)

M <- cor(alltraitsV,use="pairwise.complete.obs")

corrplot.mixed(M, lower = "ellipse", upper="number", 
               order = "FPC", 
               tl.cex = 0.8,tl.col="black")

#range sizes area relates
#voltinism and dev time

chisq.test(alltraits$HabitatBreadth,alltraits$Habitat)

######################################################################################
#fill missing data using random forests

taxon <- alltraits[,c("Family","Suborder")]
taxon$Genus <- sapply(alltraits$Species,function(x)strsplit(x," ")[[1]][1])
taxon[] <- lapply(taxon,factor)
head(taxon)

library(missForest)
alltraitsMF <- cbind(alltraitsV,taxon)
alltraitsMF <- missForest(alltraitsMF)

save(alltraitsMF,file="alltraitsMF.RData")

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

#########################################################################################

#distribution data analysis
atlas <- read.csv("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Insects/traits/10750_2017_3495_MOESM2_ESM.csv",
                  as.is=T,dec=",")
#save grid
atlasGrid <- unique(atlas[,c("Latitude","Longitude","MGRS.WGS84")])
write.csv(atlasGrid,file="atlasGrid.csv")

#make sure first letter is a capital
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
atlas$Species.name <- sapply(atlas$Species.name,function(x)firstup(x))
atlas$Species.name[which(atlas$Species.name=="Leucorrhinia Dubia")]<-"Leucorrhinia dubia"

#check all names match with our specieslist
species$Species[!species$Species %in% atlas$Species.name]

#subset atlas data to our species
atlas <- subset(atlas, Species.name %in% species$Species)
head(atlas)

#convert units to utm
#devtools::install_git("https://git.sr.ht/~hrbrmstr/mgrs")
library(mgrs)
atlas$MGRS.y <- mgrs_to_utm(atlas$MGRS.WGS84)[,"northing"]
atlas$MGRS.x <- mgrs_to_utm(atlas$MGRS.WGS84)[,"easting"]
temp <- unique(atlas[,c("Longitude","Latitude","MGRS.y","MGRS.x")])
coordinates(temp) <- c("Longitude","Latitude")

qplot(MGRS.x,MGRS.y,data=atlas)
qplot(Longitude,Latitude,data=atlas)

#create a raster with the mgrs plots
#exlude points outside mainland europe
library(rgdal)
euro <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/Euro",
                layer="europe")
euro <- subset(euro, !CNTRY_NAME %in% c("Russia","West Bank","Turkey","Tunisia","Syria","Saudi Arabia",
                                        "Morocco","Libya","Kazakhstan","Lebanon","Kuwait","Jordan",
                                        "Iraq","Iran","Israel","Gaza Strip","Iceland","Egypt","Azerbaijan",
                                        "Algeria","Georgia","Armenia"))
plot(euro)
plot(temp,col="red",add=T)
proj4string(temp) <- crs(euro)

#subset points to this region
pointsS <- temp[!is.na(over(as(temp,'SpatialPoints'),as(euro,"SpatialPolygons"),fn=NULL)),]
plot(euro)
plot(pointsS,col="red",add=T)

#now make into a raster
plot(pointsS@data$MGRS.x,pointsS@data$MGRS.y)
