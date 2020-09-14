library(raster)
library(tidyverse)
library(tmap)

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

#get years of interest
myyears <- c(1980:2019)
myyearsC <- paste(myyears,collapse="|")

#define spring file
tfiles <- "I:/ess/03_Projects_(incl_Data)/sMon/Gis-DataData/Klima Deutschland/DWD/Spring_temp"
yearfiles <- list.files(tfiles)[sapply(list.files(tfiles),function(x)grepl(myyearsC,x))]

#for each file, read it in, and convert to a df
springTemps <- plyr::ldply(yearfiles,function(x){
  r <- raster(read.asc(paste(tfiles,x,sep="/"),gz=T))
  r <- raster::aggregate(r,fact=10,fun=median)
  projection(r) <- CRS("+init=epsg:31467")
  temp <- as.data.frame(r,xy=T)
  temp$Year <- x
  temp$Year <- gsub("grids_germany_seasonal_air_temp_mean_","",temp$Year)
  temp$Year <- as.numeric(gsub("13.asc.gz","",temp$Year))
  return(subset(temp,!is.na(layer)))
})

#get trend across all years
springTrends <- springTemps %>%
  group_by(x,y) %>%
  do(model = lm(layer ~ Year, data = .)) %>% 
  rowwise() %>% 
  broom::tidy(model) %>%
  filter(term=="Year")

#plot it
dfr <- rasterFromXYZ(springTrends[,c("x","y","estimate")]) 
plot(dfr,main="Spring temp trends")
projection(dfr) <- CRS("+init=epsg:31467")
tm_shape(dfr) +
  tm_raster(palette = "YlOrRd",style="order")


#plot 5 year means
springTemps$HalfDecade <- plyr::round_any(springTemps$Year,5,floor)

springMeans <- springTemps %>%
  group_by(x,y,HalfDecade) %>%
  summarize(temp=mean(layer,na.rm=T))%>%
  spread(HalfDecade,temp) 

#how many years in each one
springTemps %>%
  group_by(HalfDecade)%>%
  summarise(nuYears = n_distinct(Year))

#convert to a raster stack
meuse.sp = SpatialPixelsDataFrame(points = 
                                      springMeans[c("x", "y")], data = springMeans[,-c(1,2)], 
                                    proj4string = CRS("+init=epsg:31467"))
meuse.r <- as(meuse.sp, "RasterStack")
tm_shape(meuse.r) +
    tm_raster(palette = "YlOrRd",style="order",title="Spring temp",midpoint=NA)
  
#map to MTBQs
library(rgdal)
mtbqs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                 layer="MTBQ_25833")

#get years of interest
myyears <- c(1980:2019)
myyearsC <- paste(myyears,collapse="|")

#define spring files
tfiles <- "I:/ess/03_Projects_(incl_Data)/sMon/Gis-DataData/Klima Deutschland/DWD/Spring_temp"
yearfiles <- list.files(tfiles)[sapply(list.files(tfiles),function(x)grepl(myyearsC,x))]


#for each file, read it in, and convert to a df
springTemps <- plyr::ldply(yearfiles,function(x){
  r <- raster(read.asc(paste(tfiles,x,sep="/"),gz=T))
  projection(r) <- CRS("+init=epsg:31467")
  mtbqs <- spTransform(mtbqs,CRS(projection(r)))
  mtbqs$meanTemp <- raster::extract(r,mtbqs,fun=median,na.rm=T)
  mtbqs$Year <- x
  mtbqs$Year <- gsub("grids_germany_seasonal_air_temp_mean_","",mtbqs$Year)
  mtbqs$Year <- as.numeric(gsub("13.asc.gz","",mtbqs$Year))
  return(mtbqs@data)
})

