library(rgdal)
library(maptools)
library(raster)
library(rgeos)

######################################################################################################
#get datafiles

#natur raume
setwd("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data")

nr <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/Narurraeume",
              layer="naturraeume_polygon")
#"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
nr_buffer <- gBuffer(nr,byid=TRUE,width=5000)

#get german county boundaries
germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
germany <- unionSpatialPolygons(germanAdmin,rep(1, length(germanAdmin)))

#get MTBQs
mtbqs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                layer="MTBQ_25833")
mtbqsDF <- data.frame(mtbqs@data)

#get MTBs
mtbs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                 layer="MTB_25832")

########################################################################################################

#get centroid of each mtbq and use that to determin the county it is assigned to

library(sp)
mtbqs_centroids <- data.frame(coordinates(mtbqs))
coordinates(mtbqs_centroids) <- c("X1","X2")
proj4string(mtbqs_centroids) <- proj4string(mtbqs)

#add on x and y
mtbqsDF$x <- mtbqs_centroids@coords[,1]
mtbqsDF$y <- mtbqs_centroids@coords[,2]

#add on county information
mtbqs_centroids <- spTransform(mtbqs_centroids,proj4string(germanAdmin))
mtbqs_counties <- over(mtbqs_centroids,germanAdmin)
mtbqsDF$Counties <- as.character(mtbqs_counties$NAME_1)

#add on lon/lat info
mtbqsDF$lon <- mtbqs_centroids@coords[,1]
mtbqsDF$lat <- mtbqs_centroids@coords[,2]

#######################################################################################################

#but many mtbqs have NA for county based on this method

#get nr to right project
germanAdmin <- spTransform(germanAdmin,proj4string(mtbqs))

#get 500 uniformly space points in each box
library(sp)
library(plyr)
#for each polygon for list of points
all.points <- llply(mtbqs@polygons,function(x){
  my.sample<-spsample(x,500,type="regular")
  proj4string(my.sample)<-proj4string(mtbqs)
  return(my.sample)
})

#overlay with the germanAdmin
all.nr <- ldply(1:length(all.points),function(x){
  temp <- over(all.points[[x]],germanAdmin)
  temp$MTB <- mtbqs$Value[x]
  temp$Q <- mtbqs$Quadrant[x]
  return(temp)
})

#get the mode Name for each mtbq
Mode <- function(x) {
  ux <- unique(x)
  x <- x[!is.na(x)]
  if(length(x)>0){
    return(ux[which.max(tabulate(match(x, ux)))])
  }
  else{
    return(NA)
    }
}

all.nr_Mode <- ddply(all.nr,.(MTB,Q),summarise,Name=Mode(NAME_1[!is.na(NAME_1)]))

#get MTB level
MTB_County <- unique(all.nr_Mode[,c("MTB","Name")])
MTB_County <- subset(MTB_County,!is.na(Name))
all.nr_Mode$MTB_County <- MTB_County$Name[match(all.nr_Mode$MTB,MTB_County$MTB)]
all.nr_Mode$Name[is.na(all.nr_Mode$Name)] <- all.nr_Mode$MTB_County[is.na(all.nr_Mode$Name)]

#add on to the mtqbs dataframe
mtbqsDF$Counties<-all.nr_Mode$Name[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
                                      interaction(all.nr_Mode$MTB,all.nr_Mode$Q))] 
########################################################################################################

#get Natur raum for each MTBQ

#get list of natur raums
sort(unique(nr$Name))#501

#get nr to right project
nr <- spTransform(nr,proj4string(mtbqs))

#get 100 uniformly space points in each box
library(sp)
library(plyr)
#for each polygon for list of points
all.points <- llply(mtbqs@polygons,function(x){
  my.sample<-spsample(x,500,type="regular")
  proj4string(my.sample)<-proj4string(mtbqs)
  return(my.sample)
})

#overlay with the nr
all.nr <- ldply(1:length(all.points),function(x){
  temp <- over(all.points[[x]],nr)
  temp$MTB <- mtbqs$Value[x]
  temp$Q <- mtbqs$Quadrant[x]
  return(temp)
})

#get the mode Name for each quadrant and mtqbs
all.nr_Mode <- ddply(all.nr,.(MTB,Q),summarise,
                     Name=Mode(Name[!is.na(Name)]), 
                     FolderPath=Mode(FolderPath[!is.na(FolderPath)]))

#add on to the mtqbs dataframe
mtbqsDF$Natur<-all.nr_Mode$Name[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
                                                  interaction(all.nr_Mode$MTB,all.nr_Mode$Q))] 

mtbqsDF$CoarseNatur<-all.nr_Mode$FolderPath[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
                                      interaction(all.nr_Mode$MTB,all.nr_Mode$Q))] 

#there are quite a few NAs

#####################################################################################################

#get Natur raum of each MTB

#get nr to right project
nr <- spTransform(nr,proj4string(mtbs))

#get 100 uniformly space points in each box
#for each polygon for list of points
all.points <- llply(mtbs@polygons,function(x){
  my.sample<-spsample(x,500,type="regular")
  proj4string(my.sample)<-proj4string(mtbs)
  return(my.sample)
})

#overlay with the nr
all.nr <- ldply(1:length(all.points),function(x){
  temp <- over(all.points[[x]],nr)
  temp$MTB <- mtbs$Value[x]
  return(temp)
})

#get the mode Name for each quadrant and mtqbs
all.nr_Mode <- ddply(all.nr,.(MTB),summarise,
                     Name=Mode(Name[!is.na(Name)]),
                     FolderPath=Mode(FolderPath[!is.na(FolderPath)]))

mtbqsDF$MTB_Natur<-all.nr_Mode$Name[match(mtbqsDF$Value,all.nr_Mode$MTB)]
mtbqsDF$MTB_CoarseNatur<-all.nr_Mode$FolderPath[match(mtbqsDF$Value,all.nr_Mode$MTB)]

###################################################################################################
 
#sort MTBQs
mtbqsDF$Q <- NA
mtbqsDF$Q[which(mtbqsDF$Quadrant=="NW")]<-1
mtbqsDF$Q[which(mtbqsDF$Quadrant=="NO")]<-2
mtbqsDF$Q[which(mtbqsDF$Quadrant=="SW")]<-3
mtbqsDF$Q[which(mtbqsDF$Quadrant=="SO")]<-4
mtbqsDF$MTB_Q <- paste0(as.character(mtbqsDF$Value),
                        as.character(mtbqsDF$Q))
library(ggplot2)
qplot(x,y,data=mtbqsDF,color=MTB_Natur)+
  theme(legend.position="none")

###################################################################################################

save(mtbqsDF,file="mtbqsDF.RData")

###################################################################################################

#work on higher-level grouping possibilities:

# #50 km grid
library(raster)
newRaster <- raster(extent(mtbqs))
res(newRaster) <- 50000
newRaster[] <- 1:ncell(newRaster)
projection(newRaster) <- proj4string(mtbqs)
plot(newRaster)
germanAdmin <- spTransform(germanAdmin,proj4string(mtbqs))
plot(germanAdmin,add=T)
#writeRaster(newRaster,file='km50grid.tif',format="GTiff",overwrite=T)

#add to mtbq to each one
mtbqs_centroids <- spTransform(mtbqs_centroids,proj4string(mtbqs))
mtbqsDF$km50 <- extract(newRaster,mtbqs_centroids)

#100 km grid
library(raster)
newRaster <- raster(extent(mtbs))
res(newRaster) <- 100000
newRaster[] <- 1:ncell(newRaster)
projection(newRaster) <- proj4string(mtbs)
plot(newRaster)
germanAdmin <- spTransform(germanAdmin,proj4string(mtbs))
plot(germanAdmin,add=T)
#writeRaster(newRaster,file='km100grid.tif',format="GTiff",overwrite=T)

#add to mtbq to each one
mtbqs_centroids <- spTransform(mtbqs_centroids,proj4string(mtbs))
mtbqsDF$km100 <- extract(newRaster,mtbqs_centroids)

length(unique(mtbqsDF$km50))#191
length(unique(mtbqsDF$km100))#54
# 
# #################################################################################################
# 
# #hexagonal grid
# #https://stackoverflow.com/questions/29374004/how-do-i-generate-a-hexagonal-grid-in-r
# #http://strimas.com/spatial/hexagonal-grids/
# require(sp)
# library(rgeos)
# 
# #making hexagonal polygons
# meuse.sr = germany
# meuse.sr = spTransform(germany,CRS(proj4string(mtbqs)))
# plot(meuse.sr)
# meuse.large = gBuffer(meuse.sr, width = 200000)
# HexPts <-spsample(meuse.large, type="hexagonal", cellsize=100000)
# HexPols <- HexPoints2SpatialPolygons(HexPts)
# plot(HexPols[meuse.sr,], add=TRUE)
# 
# #make data frame
# Spol <- HexPols
# #Creating a dataframe with Spol IDs
# Spol_df<- as.data.frame(sapply(slot(Spol, "polygons"), function(x) slot(x, "ID")))
# #Making the IDs row names 
# row.names(Spol_df) <- sapply(slot(Spol, "polygons"), function(x) slot(x, "ID"))
# # Making the spatial polygon data frame
# HexPols <- SpatialPolygonsDataFrame(Spol, data =Spol_df)
# 
# #write shape file
# writeOGR(HexPols, dsn = "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data",
#          layer = "HexPols",driver="ESRI Shapefile")
# 
# germanyDF <- as(germany,'SpatialPolygonsDataFrame')
# writeOGR(germanyDF, dsn = "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data",
#          layer = "Germany",driver="ESRI Shapefile")
# 
# #read in edited file by Volker
# HexPols <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/Hexagonal_polygons",
#                    layer="HexPolsEdit")
# germany <- spTransform(germany,CRS(proj4string(HexPols)))
# plot(germany)
# plot(HexPols,add=T)
# 
# #crop to germany
# HexPols <- gIntersection(HexPols,gUnaryUnion(mtbqs),byid=TRUE)
# plot(HexPols)
# #label the polygons
# nums <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "ID"))
# Spol <- HexPols
# Spol_df<- as.data.frame(sapply(slot(Spol, "polygons"), function(x) slot(x, "ID")))
# row.names(Spol_df) <- sapply(slot(Spol, "polygons"), function(x) slot(x, "ID"))
# HexPols <- SpatialPolygonsDataFrame(Spol, data =Spol_df)
# 
# #for each mtbq, get the polygon
# getmode <- function(v) {
#   v <- v[!is.na(v)]
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }
# library(plyr)
# mtbqsPolygons <- ldply(1:nrow(mtbqs),
#                        function(x){
#                          mtbqsCentres = spsample(mtbqs[x,],100,type="regular")
#                          Polygon = getmode(over(mtbqsCentres,HexPols))
#                          data.frame(Name = mtbqs@data$Name[x],
#                                     Value =mtbqs@data$Value[x],
#                                     Quadrant =mtbqs@data$Quadrant[x],
#                                     Polygon)
#                        })
# mtbqsDF$HexPol <- mtbqsPolygons$Polygon[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
#                                               interaction(mtbqsPolygons$Value,mtbqsPolygons$Quadrant))]
# 
# #get centroids of the polygons
# id <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "ID"))
# HexPol_x <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "labpt")[1])
# HexPol_y <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "labpt")[2])
# mtbqsDF$HexPol_x <- HexPol_x[match(mtbqsDF$HexPol,id)]
# mtbqsDF$HexPol_y <- HexPol_y[match(mtbqsDF$HexPol,id)]

###################################################################################################

#fill in remaining missing county data manually

out <- subset(mtbqsDF,is.na(Counties))
mtbqsDF$Counties[mtbqsDF$Value=="1226"]<-"Schleswig-Holstein"
mtbqsDF$Counties[mtbqsDF$Value=="1633"]<-"Schleswig-Holstein"
mtbqsDF$Counties[mtbqsDF$Value=="1815"]<-"Schleswig-Holstein"
mtbqsDF$Counties[mtbqsDF$Value=="1816"]<-"Schleswig-Holstein"
mtbqsDF$Counties[mtbqsDF$Value=="2452"]<-"Mecklenburg-Vorpommern"
mtbqsDF$Counties[mtbqsDF$Value=="4003"]<-"Nordrhein-Westfalen"
mtbqsDF$Counties[mtbqsDF$Value=="4455"]<-"Sachsen"
mtbqsDF$Counties[mtbqsDF$Value=="5155"]<-"Sachsen"
mtbqsDF$Counties[mtbqsDF$Value=="5702"]<-"Rheinland-Pfalz"
mtbqsDF$Counties[mtbqsDF$Value=="5840"]<-"Nordrhein-Westfalen"
mtbqsDF$Counties[mtbqsDF$Value=="7349"]<-"Bayern"
mtbqsDF$Counties[mtbqsDF$Value=="8439"]<-"Bayern"
mtbqsDF$Counties[mtbqsDF$Value=="8535"]<-"Bayern"
unique(mtbqsDF$Counties)

#city states to surrounding ones 
mtbqsDF$Counties[mtbqsDF$Counties=="Hamburg"]<-"Schleswig-Holstein"
mtbqsDF$Counties[mtbqsDF$Counties=="Bremen"]<-"Niedersachsen"
mtbqsDF$Counties[mtbqsDF$Counties=="Berlin"]<-"Brandenburg"

##################################################################################################

save(mtbqsDF,file="mtbqsDF.RData")

###################################################################################################

load("mtbqsDF.RData")
table(mtbqsDF$Counties)# - usually 500 ish per state

###################################################################################################

#check mtbqs for which we have data but no coordinates
#63012,50561,51561
subset(mtbqsDF,Value=="6301")
subset(mtbqsDF,Value=="5056")
subset(mtbqsDF,Value=="5156")

mtbsDF <- mtbs@data
subset(mtbsDF,Value=="6301")
subset(mtbsDF,Value=="5056")
subset(mtbsDF,Value=="5156")
#dont exist....

###################################################################################################

#formatting terrestrial ecoregion data
tdir <- "I:/ess/03_Projects_(incl_Data)/sMon/Gis-DataData/TEOW/official_teow/official"
ecoregions <- readOGR(dsn = tdir,
                      layer= "wwf_terr_ecos")
ecoregionsC <- crop(ecoregions,extent(germany))
plot(ecoregionsC)
head(ecoregionsC)
unique(ecoregionsC$ECO_NAME)#there are 5

#match with MTBQ
ecoregions <- spTransform(ecoregions,CRS(proj4string(mtbqs)))
mtbqsEco <- over(mtbqs,ecoregions)
mtbqsEco <- cbind(mtbqsEco,mtbqs@data)
mtbqsEco$ECO_NAME <- factor(mtbqsEco$ECO_NAME)
table(mtbqsEco$ECO_NAME)

#is there missing data
nrow(subset(mtbqsEco,is.na(ECO_NAME)))
head(subset(mtbqsEco,is.na(ECO_NAME)))

#get coords for these
mtbqsEco <- cbind(mtbqsEco,coordinates(mtbqs))
names(mtbqsEco)[25:26]<-c("x","y")

#plot
library(ggplot2)
ggplot(mtbqsEco)+
  geom_point(aes(x=x,y=y,colour=ECO_NAME))

###################################################################################################

#deal with missing values

#fill in MTB blanks with MTBQ data
lapply(mtbqsDF,function(x)sum(is.na(x)))
mtbqsDF$MTB_Natur[is.na(mtbqsDF$MTB_Natur)] <- mtbqsDF$Natur[is.na(mtbqsDF$MTB_Natur)]
mtbqsDF$MTB_CoarseNatur[is.na(mtbqsDF$MTB_CoarseNatur)] <- mtbqsDF$CoarseNatur[is.na(mtbqsDF$MTB_CoarseNatur)]
#1633 NW

#and vice versa
mtbqsDF$Natur[is.na(mtbqsDF$Natur)] <- mtbqsDF$MTB_Natur[is.na(mtbqsDF$Natur)]
mtbqsDF$CoarseNatur[is.na(mtbqsDF$CoarseNatur)] <- mtbqsDF$MTB_CoarseNatur[is.na(mtbqsDF$CoarseNatur)]

#NAs are all around the edges...
out <- subset(mtbqsDF,is.na(Natur))
qplot(x,y,data=out)
out <- subset(mtbqsDF,is.na(MTB_Natur))
qplot(x,y,data=out)
out <- subset(mtbqsDF,is.na(CoarseNatur))
qplot(x,y,data=out)
out <- subset(mtbqsDF,is.na(MTB_CoarseNatur))
qplot(x,y,data=out)

mtbqsDF <- subset(mtbqsDF,!is.na(Natur))
mtbqsDF <- subset(mtbqsDF,!is.na(MTB_Natur))
mtbqsDF <- subset(mtbqsDF,!is.na(CoarseNatur))
mtbqsDF <- subset(mtbqsDF,!is.na(MTB_CoarseNatur))

################################################################################

#group fine-scale naturraum to the higher level
load("mtbqsDF.RData")

library(sf)
nr <- st_as_sf(nr)
plot(nr)

#get next order file
nr_orders <- read.csv("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/Narurraeume/DWD_naturraum.csv")

#trim white space
library(gdata)
nr_orders$NR_order1 <- trim(nr_orders$NR_order1)
nr_orders$NR_order2 <- trim(nr_orders$NR_order2)

#change accents and fix mistakes
nr_orders$NR_order2 <- gsub("ä","ae",nr_orders$NR_order2)
nr_orders$NR_order2 <- gsub("ü","ue",nr_orders$NR_order2)
nr_orders$NR_order2 <- gsub("ö","oe",nr_orders$NR_order2)
nr_orders$NR_order2 <- gsub("Ö","oe",nr_orders$NR_order2)
nr_orders$NR_order2 <- gsub("Nordriesische","Nordfriesische",nr_orders$NR_order2)

#check that it finds the first fine-scale naturraum
"Allgaeuer Hochalpen"
nr_orders$NR_order1[grepl("Allgaeuer Hochalpen",nr_orders$NR_order2)]

myfun <- function(x){
  if(length(nr_orders$NR_order1[grepl(x,nr_orders$NR_order2)])==1){
    return(nr_orders$NR_order1[grepl(x,nr_orders$NR_order2)])
  }
  else{ NA
    }
}
myfun("Allgaeuer Hochalpen")
myfun("Halligen")

mtbqsDF$mtb_nr_order2 <- sapply(mtbqsDF$MTB_Natur,function(x)myfun(x))
unique(subset(mtbqsDF,is.na(mtb_nr_order2) & !is.na(mtbqsDF$MTB_Natur))[,"MTB_Natur"])
#many missing values

#match on first 3 and last 3 characters to see if we get more matches
x <- "Allgaeuer Hochalpen"
x <- "Nordfriesische Marsch"
firstChar <- substr(x,1,4)
nChar <- nchar(x)
endChar <- substr(x,nChar-4,nChar)
mySet <- grep(paste0(firstChar,".*",endChar),nr_orders$NR_order2, value=T)
nr_orders$NR_order1[grepl(mySet,nr_orders$NR_order2)]

myfun <- function(x){
  
  if(!is.na(x)){
    
  x <- as.character(x)
  firstChar <- substr(x,1,3)
  nChar <- nchar(x)
  endChar <- substr(x,nChar-2,nChar)
  mySet <- grep(paste0(firstChar,".*",endChar),
                nr_orders$NR_order2, value=T)
  
  
  mySet_len <- nchar(paste(mySet,collapse = ""))
   
  if(mySet_len< 10 ){
     mySet <- substr(mySet,1,mySet_len)
  }
   else{
     mySet <- substr(mySet,1,7)  
  }
  
  if(length(mySet)==1){
    out <- nr_orders$NR_order1[grepl(mySet,nr_orders$NR_order2)]
    if(length(out)==1){
    return(out)
    }else{
      return(NA)
    }
  }else{ 
    return(NA)
  }
  
  }
  else{
    return(NA)
  }
  
}

mtbqsDF$mtb_nr_order2_try2<- sapply(mtbqsDF$MTB_Natur,function(x)myfun(x))

unique(subset(mtbqsDF,is.na(mtb_nr_order2_try2) & !is.na(mtbqsDF$MTB_Natur))[,c("MTB_Natur","mtb_nr_order2")])

#some matched twice
mtbqsDF$mtb_nr_order2_try2[is.na(mtbqsDF$mtb_nr_order2_try2)] <- mtbqsDF$mtb_nr_order2[is.na(mtbqsDF$mtb_nr_order2_try2)]
unique(subset(mtbqsDF,is.na(mtb_nr_order2_try2) & !is.na(mtbqsDF$MTB_Natur))[,c("MTB_Natur")])#still 57 missing ones...

tempDF <- unique(mtbqsDF[,c("MTB_Natur","mtb_nr_order2_try2")])
write.csv(tempDF,file="nr_order2.csv",row.names=FALSE)
#these were fixed manually

#read back in
tempDF <- read.csv("nr_order2.csv")
nr$order2 <- tempDF$mtb_nr_order2_try2[match(nr$Name,tempDF$MTB_Natur)]
nr$order2 <- trim(nr$order2)
nr$area <- st_area(nr)

myOrders <- sort(unique(nr$order2))

library(tidyverse)
nr2 <-
  nr %>%
  group_by(order2) %>%
  summarise(area = sum(area))

plot(nr2)

#add to mtbqDF
mtbqsDF$MTB_MidNatur <- tempDF$mtb_nr_order2_try2[match(mtbqsDF$MTB_Natur,tempDF$MTB_Natur)]
nrow(subset(mtbqsDF,is.na(MTB_MidNatur)))
#0
nrow(subset(mtbqsDF,is.na(MTB_CoarseNatur)))
#0

save(mtbqsDF,file="mtbqsDF.RData")
#check 5101
subset(mtbqsDF,Value==5101)

####################################################################################################

## get MTB x and y
mtbs <- spTransform(mtbs,CRS(proj4string(mtbqs)))
mtbDF <- as.data.frame(mtbs)
mtbDF$x <- coordinates(mtbs)[,1]
mtbDF$y <- coordinates(mtbs)[,2]

#add it to the MTBQ data frame
mtbqsDF$x_MTB <- mtbDF$x[match(mtbqsDF$Value,mtbDF$Value)]
mtbqsDF$y_MTB <- mtbDF$y[match(mtbqsDF$Value,mtbDF$Value)]
save(mtbqsDF,file="mtbqsDF.RData")
