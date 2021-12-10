#Extend MTBQ buffer grid beyond edge to avoid buffer effects
library(tidyverse)
library(sf)
library(tmap)

### get data #### 

#mtbqs
mtbqs <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",layer="MTBQ_25833")

mtbqs$Q <- NA
mtbqs$Q[which(mtbqs$Quadrant=="NW")]<-1
mtbqs$Q[which(mtbqs$Quadrant=="NO")]<-2
mtbqs$Q[which(mtbqs$Quadrant=="SW")]<-3
mtbqs$Q[which(mtbqs$Quadrant=="SO")]<-4
mtbqs$MTB_Q <- paste0(as.character(mtbqs$Value),
                        as.character(mtbqs$Q))

#naturaum
nr <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Narurraeume",
              layer="naturraeume_polygon")


#remove small island in north sea
nr <- subset(nr,Name!="Helgoland")

### coarse naturraum ####

nr <- st_transform(nr, st_crs(mtbqs))

nr <- nr %>%
      dplyr::group_by(FolderPath) %>%
      dplyr::summarise() %>%
      rename(Coarsenaturraum = FolderPath)

st_crs(mtbqs)==st_crs(nr)

### nr of mtbqs points ####

mtbqs_points <- st_cast(mtbqs,"POLYGON") %>%
                  st_centroid()
# crs = 25833
nrow(mtbqs_points)

tm_shape(nr)+
  tm_polygons()+
tm_shape(mtbqs_points)+
  tm_dots()
#ones at the edge dont overlap

mtbqs_points1 <- st_join(mtbqs_points, nr, left = TRUE) %>% #overlap in utm
                dplyr::filter(!is.na(Coarsenaturraum))
mtbqs_points1 <- subset(mtbqs_points1,!Value %in% c(1815,1816))#north sea islands
nrow(mtbqs_points1)

### nr buffer of mtbqs points ####

nr_buffer15 <- st_buffer(nr, dist = 12500)

tm_shape(nr_buffer15)+
  tm_polygons()+
  tm_shape(mtbqs_points)+
  tm_dots()
#few dont overlap but all in the north

mtbqs_points2 <- st_join(mtbqs_points, nr_buffer15, left = TRUE)# overlap on utm scale

#remove the 4 blanks - islands in the north east
mtbqs_points2 <- subset(mtbqs_points2, !is.na(Coarsenaturraum))
mtbqs_points2 <- subset(mtbqs_points2,!Value %in% c(1815,1816))
mtbqs_points2 <- subset(mtbqs_points2,!duplicated(MTB_Q))
nrow(mtbqs_points2)

### extend points ####

#make raster to extend as regular grid 
library(raster)

mtbqs_points <- st_transform(mtbqs_points,crs = 4326)
mtbqs_raster <- SpatialPixelsDataFrame(points = st_coordinates(mtbqs_points),
                                       data = data.frame(MTB_Q = as.numeric(mtbqs_points$MTB_Q)),
                                       tolerance = 0.523346,
                                       proj4string = CRS("+init=epsg:4326"))
mtbqs_raster <- raster(mtbqs_raster)
plot(mtbqs_raster)
mtbqs_raster_buffer <- raster::buffer(mtbqs_raster, width = 50000)
plot(mtbqs_raster_buffer)

#see how new points compare with old points
tm_shape(mtbqs_raster_buffer)+
  tm_raster()+
  tm_shape(mtbqs_points)+
  tm_dots()
#looks good! 

#get new list of points
mtbqs_raster_buffer_points <- as.data.frame(mtbqs_raster_buffer, xy=T) %>%
                              filter(layer==1)
head(mtbqs_raster_buffer_points)

### nr of raster points ####

mtbqs_raster_buffer_pointsS <- st_as_sf(mtbqs_raster_buffer_points,
                                        coords =c("x","y"),
                                        crs = 4326)
                                
#check overlap again
tm_shape(mtbqs_raster_buffer_pointsS)+
  tm_dots()+
  tm_shape(mtbqs_points)+
  tm_dots(col="red")

#nr_buffer1 <- st_buffer(nr, dist = 1000)
#nr_buffer_dissolve <- dplyr::summarise(nr)
#nr_buffer50 <- nr %>%
#              st_buffer(., dist = 50000) %>% 
#              st_union(.,by_feature = TRUE) %>% 
#              st_difference(.,nr_buffer_dissolve)
#st_write(nr_buffer50,"nc.shp")

#following the above lines, a few tweaks were made by Volker in QGIS 
nr_buffer50 <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Narurraeume/buffer50km",
                       layer="nc_manipulated")

nr_buffer50$Coarsenaturraum <- NA
nr_buffer50$Coarsenaturraum[nr_buffer50$code==1] <- "Naturräume Deutschland/Alpen"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==3] <- "Naturräume Deutschland/Nordostdeutsches Tiefland"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==4] <- "aturräume Deutschland/Nordwestdeutsches Tiefland"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==5] <- "Naturräume Deutschland/Oestliches Mittelgebirge"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==6] <- "Naturräume Deutschland/Suedwestdeutsches Mittelgebirge / Stufenland"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==7] <- "Naturräume Deutschland/Westliches Mittelgebirge"
nr_buffer50 <- nr_buffer50 %>% group_by(Coarsenaturraum) %>% summarise()

#again check
tm_shape(nr)+
  tm_polygons(col="blue")+
tm_shape(nr_buffer50)+
  tm_polygons(col="red")

mtbqs_points3 <- st_transform(mtbqs_raster_buffer_pointsS,
                                            crs = st_crs(nr)) %>%
                                st_intersection(.,nr_buffer50) %>%
                filter(!is.na(Coarsenaturraum )) %>%
                filter(!duplicated(geometry))

#check overlap
tm_shape(nr_buffer50)+
  tm_polygons()+
  tm_shape(mtbqs_points3)+
  tm_dots()

head(mtbqs_points3)

### combine all ####

#nr of centroids
head(mtbqs_points1)
mtbqs_points1c <- mtbqs_points1 %>%
                  st_set_geometry(.,NULL) %>%
                  as_tibble() %>%
                  dplyr::select(c("Value","MTB_Q","Coarsenaturraum")) %>%
                  dplyr::mutate(x = st_coordinates(mtbqs_points1)[,1],
                                y = st_coordinates(mtbqs_points1)[,2])%>%
                  add_column(type = "original")

ggplot(mtbqs_points1c)+
  geom_point(aes(x=x,y=y, colour=Coarsenaturraum),size=0.5,alpha=0.5)+
  coord_fixed()

#nr for missing points
head(mtbqs_points2)
mtbqs_points2c <- mtbqs_points2 %>%
                st_set_geometry(.,NULL) %>%
                as_tibble() %>%
                dplyr::select(c("Value","MTB_Q","Coarsenaturraum")) %>%
                dplyr::mutate(x = st_coordinates(mtbqs_points2)[,1],
                              y = st_coordinates(mtbqs_points2)[,2])%>%
                dplyr::filter(!MTB_Q %in% mtbqs_points1c$MTB_Q)%>%
                add_column(type = "edges")

ggplot(mtbqs_points2c)+
  geom_point(aes(x=x,y=y, colour=Coarsenaturraum),size=0.5,alpha=0.5)+
  coord_fixed()

#nr for extended points
head(mtbqs_points3)
mtbqs_points3c <- mtbqs_points3 %>%
                  st_set_geometry(.,NULL) %>%
                  as_tibble() %>%
                  dplyr::select(c("Coarsenaturraum")) %>%
                  dplyr::mutate(x = st_coordinates(mtbqs_points3)[,1],
                                y = st_coordinates(mtbqs_points3)[,2]) %>%
                 add_column(type = "extension")

                
#combine these two
mtbqs_pointsc <- bind_rows(mtbqs_points1c, mtbqs_points2c,mtbqs_points3c)

ggplot(mtbqs_pointsc)+
  geom_point(aes(x=x,y=y, colour=type),size=0.5,alpha=0.5)+
  coord_fixed()

table(mtbqs_pointsc$type)
#looks good!!
