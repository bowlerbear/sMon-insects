#landuseData

### germany-wide corine ########################################

#data frame
load("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/mtbqsDF.RData")

#also shape file
mtbqs <- st_read("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile/MTBQ_25833.shp")

#germany
germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")
germanAdmin <- st_as_sf(germanAdmin)
germanAdmin <- st_transform(germanAdmin,crs=st_crs(mtbqs))
germanAdmin <- st_union(germanAdmin)
german_buffer <- st_buffer(germanAdmin,dist=10000)

#1990 data
r <- raster("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Corine/CLC_1990/u2000_clc1990_v2020_20u1_raster100m/DATA/U2000_CLC1990_V2020_20u1.tif")
german_buffer <- st_transform(german_buffer,crs=st_crs(r,asText=TRUE))
german_buffer <- as(german_buffer,'Spatial')
r <- crop(r,german_buffer)
plot(r)

#extract for each mtbq
mtbqsTrans <- mtbqs %>% 
  st_transform(crs(r)) %>% 
  as(.,'Spatial')

circlesLU <- raster::extract(r,mtbqsTrans,df=T)
names(circlesLU)[2] <- "Land.cover.class"
circlesLU_2 <- ddply(circlesLU,.(ID),
                     summarise,
                     urban = mean(Land.cover.class %in% c(1:6),na.rm=T),
                     water = mean(Land.cover.class %in% c(40:41),na.rm=T),
                     tree = mean(Land.cover.class %in% c(23:25),na.rm=T),
                     grass = mean(Land.cover.class %in% c(18,26:28),na.rm=T),
                     crop =  mean(Land.cover.class %in% c(12:17,19:21),na.rm=T),
                     tot = length(Land.cover.class))

circlesLU_2$Value <- mtbqs$Value
circlesLU_2$Quadrant <- mtbqs$Quadrant


#repeat for 2018
r <- raster("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Corine/CLC_2018/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")
