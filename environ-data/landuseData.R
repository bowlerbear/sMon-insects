### germany extent #############################################

library(sf)
library(raster)

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

### corine ########################################

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

#### esa ###############################################

#2016-2018

years <- c(2016,2017,2018)

#get urban land cover for each year and MTBQ
ldir <- "C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/landUsePerMTBQ/esaCciLandCover_2016_2018/dataset-satellite-land-cover-7ceb489b-98a0-4153-b5f6-809c82b288c1"

myfun <- function(year){
  myfile <- list.files(ldir)[grep(year,list.files(ldir))]
  r <- raster(paste(ldir,myfile,sep="/"))
  german_buffer <- st_transform(german_buffer,crs = st_crs(r,asTest=TRUE))
  Germany_sp <- as(german_buffer,'Spatial')
  r <- crop(r,Germany_sp)
  
  #extract for each mtbq
  mtbqsTrans <- mtbqs %>% 
    st_transform(crs(r)) %>% 
    as(.,'Spatial')
  
  circlesLU <- raster::extract(r,mtbqsTrans,df=T)
  circlesLU <- plyr::ddply(circlesLU,"ID",
                     summarise,
                     urban = mean(Land.cover.class.defined.in.LCCS==190),
                     water = mean(Land.cover.class.defined.in.LCCS==210),
                     wetland = mean(Land.cover.class.defined.in.LCCS==180),
                     tree = mean(Land.cover.class.defined.in.LCCS %in% c(50,60:62,70:72,80:82,90,100,160,170)),
                     grass = mean(Land.cover.class.defined.in.LCCS %in% c(110,130)),
                     crop =  mean(Land.cover.class.defined.in.LCCS %in% c(10:12,20,30,40)),
                     tot = length(Land.cover.class.defined.in.LCCS))
  
  #add year
  circlesLU$Year <- year
  
  #return output
  return(circlesLU)
  
}

#apply function to each year
library(tidyverse)
lateYearsDF <- years %>%
  map(myfun) %>%
  bind_rows()
saveRDS(lateYearsDF,
        file = "environ-data/lateYearsDF_environData.rds")

# 1992-2015 
#run on RStudio server
years <- c(1992:2015)

#get urban land cover for each year and MTBQ

ldir <- "I:/ess/03_Projects_(incl_Data)/sMon/Gis-DataData/landUsePerMTBQ/esaCciLandCover_1992_2015"

myfun <- function(year){
  myfile <- list.files(ldir)[grep(year,list.files(ldir))]
  r <- raster(paste(ldir,myfile,sep="/"))
  german_buffer <- st_transform(german_buffer,crs = st_crs(r,asTest=TRUE))
  Germany_sp <- as(german_buffer,'Spatial')
  r <- crop(r,Germany_sp)
  
  #extract for each mtbq
  mtbqsTrans <- mtbqs %>% 
    st_transform(crs(r)) %>% 
    as(.,'Spatial')
  circlesLU <- raster::extract(r,mtbqsTrans,df=T)
  names(circlesLU)[2] <- "Land.cover.class.defined.in.LCCS"
  circlesLU <- plyr::ddply(circlesLU,"ID",
                     summarise,
                     urban = mean(Land.cover.class.defined.in.LCCS==190),
                     water = mean(Land.cover.class.defined.in.LCCS==210),
                     tree = mean(Land.cover.class.defined.in.LCCS %in% c(50,60:62,70:72,80:82,90,100,160,170)),
                     grass = mean(Land.cover.class.defined.in.LCCS %in% c(110,130)),
                     crop =  mean(Land.cover.class.defined.in.LCCS %in% c(10:12,20,30,40)),
                     tot = length(Land.cover.class.defined.in.LCCS))
  
  #add year
  circlesLU$Year <- year
  
  #return output
  return(circlesLU)
  
}

#apply function to each year
firstYearsDF <- years %>%
  map(myfun) %>%
  bind_rows()


### combine land use data ####################################

#landUseDF <- rbind(firstYearsDF,lateYearsDF)
df1 <- readRDS("environ-data/firstYearsDF_90s.rds")
df2 <- readRDS("environ-data/firstYearsDF_00s.rds")
df3 <- readRDS("environ-data/firstYearsDF_10s.rds")
df4 <- readRDS("environ-data/lateYearsDF_environData.rds")[,-4]
df <- rbind(df1,df2,df3,df4)

#add mtbqs
df <- cbind(df,mtbqsTrans)
saveRDS(df,file="environ-data/esacci_MTBQ.rds")

#check
df$x <- coordinates(mtbqsTrans)[,1]
df$y <- coordinates(mtbqsTrans)[,2]
library(ggplot2)
qplot(x,y,data=subset(df,Year==1992),colour=urban)

###end ########################################################
