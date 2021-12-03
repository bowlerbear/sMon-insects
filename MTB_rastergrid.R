#Extend mtb buffer grid beyond edge to avoid buffer effects
library(tidyverse)
library(sf)
library(tmap)

### get data #### 

#mtbs
mtbs <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                layer="MTB_25832")

mtbs<- subset(mtbs,!Value %in% c(1815,1816))#north sea islands

#naturaum
nr <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Narurraeume",
              layer="naturraeume_polygon")


#remove small island in north sea
nr <- subset(nr,Name!="Helgoland")

### coarse naturraum ####

nr <- st_transform(nr, st_crs(mtbs))

nr <- nr %>%
      dplyr::group_by(FolderPath) %>%
      dplyr::summarise() %>%
      rename(Coarsenaturraum = FolderPath)

st_crs(mtbs)==st_crs(nr)

### nr of mtbs points ####

mtbs_points <- st_cast(mtbs,"POLYGON") %>%
                  st_centroid()
# crs = 25833
nrow(mtbs_points)

tm_shape(nr)+
  tm_polygons()+
tm_shape(mtbs_points)+
  tm_dots()
#ones at the edge dont overlap

mtbs_points1 <- st_join(mtbs_points, nr, left = TRUE) %>% #overlap in utm
                dplyr::filter(!is.na(Coarsenaturraum))
nrow(mtbs_points1)

### nr buffer of mtbs points ####

nr_buffer15 <- st_buffer(nr, dist = 8000)

tm_shape(nr_buffer15)+
  tm_polygons()+
  tm_shape(mtbs_points)+
  tm_dots()
#few dont overlap but all in the north

mtbs_points2 <- st_join(mtbs_points, nr_buffer15, left = TRUE)# overlap on utm scale

#remove the 4 blanks - islands in the north east
mtbs_points2 <- subset(mtbs_points2, !is.na(Coarsenaturraum))
mtbs_points2 <- subset(mtbs_points2,!duplicated(Value))
nrow(mtbs_points2)

### extend points ####

#make raster to extend as regular grid 

library(raster)

mtbs_points <- st_transform(mtbs_points,crs = 4326)#have to make it on the lat/lon scale
mtbs_raster <- SpatialPixelsDataFrame(points = st_coordinates(mtbs_points),
                                       data = data.frame(Value = as.numeric(mtbs_points$Value)),
                                       tolerance = 0.523346,
                                       proj4string = CRS("+init=epsg:4326"))
mtbs_raster <- raster(mtbs_raster)

#ideal extent of new raster
mtbs_raster_extend <- extend(mtbs_raster, 10, value=NA)
plot(mtbs_raster_extend)
mtbs_raster_buffer <- raster::buffer(mtbs_raster_extend, width = 150000)#just to be sure
plot(mtbs_raster_buffer)

#see how new points compare with old points
tm_shape(mtbs_raster_buffer)+
  tm_raster()+
  tm_shape(mtbs_points)+
  tm_dots()
#looks good! 

#get new list of points
mtbs_raster_buffer_points <- as.data.frame(mtbs_raster_buffer, xy=T) %>%
                              filter(layer==1)
head(mtbs_raster_buffer_points)

mtbs_raster_buffer_pointsS <- st_as_sf(mtbs_raster_buffer_points,
                                        coords =c("x","y"),
                                        crs = 4326)#lat lon scale
                                
#check overlap again
tm_shape(mtbs_raster_buffer_pointsS)+
  tm_dots()+
  tm_shape(mtbs_points)+
  tm_dots(col="red")

#now get nr buffer with 50 km extension

#nr_buffer_dissolve <- dplyr::summarise(nr)
#nr_buffer50 <- nr %>%
#              st_buffer(., dist = 50000, endCapStyle = "SQUARE") %>% 
#              st_union(.,by_feature = TRUE) %>% 
#              st_difference(.,nr_buffer_dissolve)

# tm_shape(nr)+
#   tm_polygons("Coarsenaturraum")+
# tm_shape(nr_buffer50)+
#   tm_polygons("Coarsenaturraum")

#st_write(nr_buffer50,"nc.shp")

#following the above lines, a few tweaks were made by Volker in QGIS 
nr_buffer50 <- st_read(dsn="C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Narurraeume/buffer50km",
                       layer="nc_manipulated")

nr_buffer50$Coarsenaturraum <- NA
nr_buffer50$Coarsenaturraum[nr_buffer50$code==1] <- "Naturräume Deutschland/Alpen"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==3] <- "Naturräume Deutschland/Nordostdeutsches Tiefland"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==4] <- "Naturräume Deutschland/Nordwestdeutsches Tiefland"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==5] <- "Naturräume Deutschland/Oestliches Mittelgebirge"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==6] <- "Naturräume Deutschland/Suedwestdeutsches Mittelgebirge / Stufenland"
nr_buffer50$Coarsenaturraum[nr_buffer50$code==7] <- "Naturräume Deutschland/Westliches Mittelgebirge"
nr_buffer50 <- nr_buffer50 %>% group_by(Coarsenaturraum) %>% summarise()

#there is still an overlap problem
nr_buffer50_Alpen <- subset(nr_buffer50,Coarsenaturraum=="Naturräume Deutschland/Alpen")
nr_buffer50_NO <- subset(nr_buffer50,Coarsenaturraum=="Naturräume Deutschland/Nordostdeutsches Tiefland")
nr_buffer50_NW <- subset(nr_buffer50,Coarsenaturraum=="Naturräume Deutschland/Nordwestdeutsches Tiefland")
nr_buffer50_O <- subset(nr_buffer50,Coarsenaturraum=="Naturräume Deutschland/Oestliches Mittelgebirge")
nr_buffer50_SW <- subset(nr_buffer50,Coarsenaturraum=="Naturräume Deutschland/Suedwestdeutsches Mittelgebirge / Stufenland")
nr_buffer50_W <- subset(nr_buffer50,Coarsenaturraum=="Naturräume Deutschland/Westliches Mittelgebirge")

#remove overlap from neighbours
nr_buffer50_Alpen <- nr_buffer50_Alpen %>% st_difference(.,nr_buffer50_SW) %>% st_difference(.,nr_buffer50_O)

nr_buffer50_O <- nr_buffer50_O %>% st_difference(.,nr_buffer50_Alpen) %>% st_difference(.,nr_buffer50_SW)

nr_buffer50_NO <- nr_buffer50_NO %>% st_difference(.,nr_buffer50_NW) %>% st_difference(.,nr_buffer50_O)

nr_buffer50_NW <- nr_buffer50_NW %>% st_difference(.,nr_buffer50_NO) %>% st_difference(.,nr_buffer50_W)

nr_buffer50_W <- nr_buffer50_W %>% st_difference(.,nr_buffer50_NW) %>% st_difference(.,nr_buffer50_SW)

nr_buffer50_SW <- nr_buffer50_SW %>% st_difference(.,nr_buffer50_W) %>% st_difference(.,nr_buffer50_Alpen)

#combine all
nr_buffer50_mod <- bind_rows(nr_buffer50_Alpen,nr_buffer50_O,nr_buffer50_NO,
                             nr_buffer50_NW,nr_buffer50_W,nr_buffer50_SW)

tm_shape(nr_buffer50_mod)+
  tm_polygons("Coarsenaturraum")

tm_shape(nr_buffer50_mod)+
  tm_polygons("Coarsenaturraum")+
  tm_facets(by="Coarsenaturraum")

#convert to utm crs
mtbs_raster_buffer_pointsS <- st_transform(mtbs_raster_buffer_pointsS, crs = st_crs(nr_buffer50_mod)) 
                             
#check overlap
tm_shape(nr_buffer50_mod)+
  tm_polygons()+
  tm_shape(mtbs_raster_buffer_pointsS)+
  tm_dots()

#get overlap i.e. nr of extended points
mtbs_points3 <- st_join(mtbs_raster_buffer_pointsS, nr_buffer50_mod)# overlap on utm scale
mtbs_points3 <- subset(mtbs_points3, !is.na(Coarsenaturraum))

head(mtbs_points3)

tm_shape(mtbs_points3)+
  tm_dots("Coarsenaturraum")
table(mtbs_points3$Coarsenaturraum)

### combine all ####

#nr of centroids
head(mtbs_points1)
mtbs_points1c <- mtbs_points1 %>%
                  st_set_geometry(.,NULL) %>%
                  as_tibble() %>%
                  dplyr::select(c("Value","Coarsenaturraum")) %>%
                  dplyr::mutate(x = st_coordinates(mtbs_points1)[,1],
                                y = st_coordinates(mtbs_points1)[,2])%>%
                  add_column(type = "original")

ggplot(mtbs_points1c)+
  geom_point(aes(x=x,y=y, colour=Coarsenaturraum),size=1,alpha=0.5)+
  coord_fixed()

#nr for missing points
head(mtbs_points2)
mtbs_points2c <- mtbs_points2 %>%
                st_set_geometry(.,NULL) %>%
                as_tibble() %>%
                dplyr::select(c("Value","Coarsenaturraum")) %>%
                dplyr::mutate(x = st_coordinates(mtbs_points2)[,1],
                              y = st_coordinates(mtbs_points2)[,2])%>%
                dplyr::filter(!Value %in% mtbs_points1c$Value)%>%
                add_column(type = "edges")

ggplot(mtbs_points2c)+
  geom_point(aes(x=x,y=y, colour=Coarsenaturraum),size=1,alpha=0.5)+
  coord_fixed()

#nr for extended points
head(mtbs_points3)
mtbs_points3c <- mtbs_points3 %>%
                  st_set_geometry(.,NULL) %>%
                  as_tibble() %>%
                  dplyr::select(c("Coarsenaturraum")) %>%
                  dplyr::mutate(x = st_coordinates(mtbs_points3)[,1],
                                y = st_coordinates(mtbs_points3)[,2]) %>%
                 add_column(type = "extension")

ggplot(mtbs_points3c)+
  geom_point(aes(x=x,y=y, colour=Coarsenaturraum),size=1,alpha=0.3)+
  coord_fixed()               
 
#combine these
mtbs_pointsc <- bind_rows(mtbs_points1c, mtbs_points2c,mtbs_points3c)

ggplot(mtbs_pointsc)+
  geom_point(aes(x=x,y=y, colour=type),size=1,alpha=0.5)+
  coord_fixed()

ggplot(mtbs_pointsc)+
  geom_point(aes(x=x,y=y, colour=Coarsenaturraum),size=1,alpha=0.5)+
  coord_fixed()

ggplot(mtbs_pointsc)+
  geom_point(aes(x=x,y=y),size=1,alpha=0.5)+
  facet_wrap(~Coarsenaturraum)+
  coord_fixed()

table(mtbs_pointsc$type)
#looks good!!

saveRDS(mtbs_pointsc,file="MTB_extendedpoints.rds")
