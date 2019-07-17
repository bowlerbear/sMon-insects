library(rgdal)
library(maptools)

######################################################################################################
#get datafiles

#natur raume
setwd("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data")

nr <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/Narurraeume",
              layer="naturraeume_polygon")
#"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#get german county boundaries
germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
germany <- unionSpatialPolygons(germanAdmin,rep(1, length(germanAdmin)))

#get MTBQs
mtbqs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                layer="MTBQ_25833")
mtbqsDF <- data.frame(mtbqs@data)

#get MTBQs
mtbs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
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
  my.sample<-spsample(x,100,type="regular")
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
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
all.nr_Mode <- ddply(all.nr,.(MTB,Q),summarise,Name=Mode(Name[!is.na(Name)]))

#add on to the mtqbs dataframe
mtbqsDF$Natur<-all.nr_Mode$Name[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
                                                  interaction(all.nr_Mode$MTB,all.nr_Mode$Q))] 

#there are quite a few NAs

#####################################################################################################

#get Natur raum of each MTB

#get nr to right project
nr <- spTransform(nr,proj4string(mtbs))

#get 100 uniformly space points in each box
#for each polygon for list of points
all.points <- llply(mtbs@polygons,function(x){
  my.sample<-spsample(x,100,type="regular")
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
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
all.nr_Mode <- ddply(all.nr,.(MTB),summarise,Name=Mode(Name[!is.na(Name)]))

mtbqsDF$Natur<-all.nr_Mode$Name[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
                                      interaction(all.nr_Mode$MTB,all.nr_Mode$Q))]

mtbqsDF$MTB_Natur<-all.nr_Mode$Name[match(mtbqsDF$Value,all.nr_Mode$MTB)]

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

#50 km grid
library(raster)
newRaster <- raster(extent(mtbqs))
res(newRaster) <- 50000
newRaster[] <- 1:ncell(newRaster)
projection(newRaster) <- proj4string(mtbqs)
plot(newRaster)
germanAdmin <- spTransform(germanAdmin,proj4string(mtbqs))
plot(germanAdmin,add=T)
writeRaster(newRaster,file='km50grid.tif',format="GTiff",overwrite=T)

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
writeRaster(newRaster,file='km100grid.tif',format="GTiff",overwrite=T)

#add to mtbq to each one
mtbqs_centroids <- spTransform(mtbqs_centroids,proj4string(mtbs))
mtbqsDF$km100 <- extract(newRaster,mtbqs_centroids)

length(unique(mtbqsDF$km50))#191
length(unique(mtbqsDF$km100))#54

#################################################################################################

#hexagonal grid
#https://stackoverflow.com/questions/29374004/how-do-i-generate-a-hexagonal-grid-in-r
#http://strimas.com/spatial/hexagonal-grids/
require(sp)
library(rgeos)

#making hexagonal polygons
meuse.sr = germany
meuse.sr = spTransform(germany,CRS(proj4string(mtbqs)))
plot(meuse.sr)
meuse.large = gBuffer(meuse.sr, width = 200000)
HexPts <-spsample(meuse.large, type="hexagonal", cellsize=100000)
HexPols <- HexPoints2SpatialPolygons(HexPts)
plot(HexPols[meuse.sr,], add=TRUE)

#make data frame
Spol <- HexPols
#Creating a dataframe with Spol IDs
Spol_df<- as.data.frame(sapply(slot(Spol, "polygons"), function(x) slot(x, "ID")))
#Making the IDs row names 
row.names(Spol_df) <- sapply(slot(Spol, "polygons"), function(x) slot(x, "ID"))
# Making the spatial polygon data frame
HexPols <- SpatialPolygonsDataFrame(Spol, data =Spol_df)

#write shape file
writeOGR(HexPols, dsn = "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data",
         layer = "HexPols",driver="ESRI Shapefile")

germanyDF <- as(germany,'SpatialPolygonsDataFrame')
writeOGR(germanyDF, dsn = "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data",
         layer = "Germany",driver="ESRI Shapefile")

#read in edited file by Volker
HexPols <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/Spatial_data",
                   layer="HexPolsEdit")
germany <- spTransform(germany,CRS(proj4string(HexPols)))
plot(germany)
plot(HexPols,add=T)

#crop to germany
HexPols <- gIntersection(HexPols,gUnaryUnion(mtbqs),byid=TRUE)
plot(HexPols)
#label the polygons
nums <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "ID"))
Spol <- HexPols
Spol_df<- as.data.frame(sapply(slot(Spol, "polygons"), function(x) slot(x, "ID")))
row.names(Spol_df) <- sapply(slot(Spol, "polygons"), function(x) slot(x, "ID"))
HexPols <- SpatialPolygonsDataFrame(Spol, data =Spol_df)

#for each mtbq, get the polygon
getmode <- function(v) {
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
library(plyr)
mtbqsPolygons <- ldply(1:nrow(mtbqs),
                       function(x){
                         mtbqsCentres = spsample(mtbqs[x,],100,type="regular")
                         Polygon = getmode(over(mtbqsCentres,HexPols))
                         data.frame(Name = mtbqs@data$Name[x],
                                    Value =mtbqs@data$Value[x],
                                    Quadrant =mtbqs@data$Quadrant[x],
                                    Polygon)
                       })
mtbqsDF$HexPol <- mtbqsPolygons$Polygon[match(interaction(mtbqsDF$Value,mtbqsDF$Quadrant),
                                              interaction(mtbqsPolygons$Value,mtbqsPolygons$Quadrant))]

#get centroids of the polygons
id <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "ID"))
HexPol_x <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "labpt")[1])
HexPol_y <- sapply(slot(HexPols, "polygons"), function(x) slot(x, "labpt")[2])
mtbqsDF$HexPol_x <- HexPol_x[match(mtbqsDF$HexPol,id)]
mtbqsDF$HexPol_y <- HexPol_y[match(mtbqsDF$HexPol,id)]

save(mtbqsDF,file="mtbqsDF.RData")

###################################################################################################
load("mtbqsDF.RData")
table(mtbqsDF$Counties)# - usually 500 ish per state
###################################################################################################