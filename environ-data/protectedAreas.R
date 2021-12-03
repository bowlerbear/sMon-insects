#protectedAreas
library(sf)
library(tidyverse)

#get proj
mtbqs <- st_read("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile","MTBQ_25833")

#add MTBQ
mtbq_fn <- function(x){
  
  x$MTB_Q <- NA  
  
  x$MTB_Q= ifelse(x$Quadrant == "NW" ,
                  paste(x$Value, 1, sep = ""),
                  ifelse(x$Quadrant == "NO",
                         paste(x$Value, 2, sep = ""),
                         ifelse(x$Quadrant == "SW", 
                                paste(x$Value, 3, sep = ""),
                                paste(x$Value, 4, sep = ""))))
  return(x)
}
mtbqs <- mtbq_fn(mtbqs)


#data from Bfn
paFiles <- list.files("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/ProtectedAreas_DE/BFN")
paFiles <- paFiles[grepl(".zip",paFiles)]
#paFiles 1 is gros naturraume
#paFiles 2 is naturraume - looks interesting - habitat info

#limit to schutz files
paFiles <- paFiles[grepl("SCH",paFiles)]

#limit to Naturschutzgebiet (NSGs) and National Parks (NLPs)
setwd("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/ProtectedAreas_DE/BFN")

#get NSG
temp <- tempfile()
unzip(zipfile = "BFN_SCH_NSG.zip", exdir = temp)
your_SHP_file<-list.files(temp, pattern = ".shp$",full.names=TRUE)
paShape <- sf::read_sf(your_SHP_file)
paShape <- st_transform(paShape,crs=st_crs(mtbqs))
paShape <- st_buffer(paShape,dist=0)
#plot(paShape)

#get NLP
temp <- tempfile()
unzip(zipfile = "BFN_SCH_NLP.zip", exdir = temp)
your_SHP_file<-list.files(temp, pattern = ".shp$",full.names=TRUE)
paShape2 <- sf::read_sf(your_SHP_file)
paShape2 <- st_transform(paShape2,crs=st_crs(mtbqs))
paShape2 <- st_buffer(paShape2,dist=0)
#plot(paShape2)
#pa_Union <- st_union(paShape,paShape2)#memory issues

#read in Netras object
pa_Union <- sf::read_sf(dsn = "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_bias/bias_template/MTBQ_data/ProtectedAreas_DE/pre-processing/protected_area",
                       layer = "pa_union")
pa_Union <- st_transform(pa_Union,crs=st_crs(mtbqs))
pa_Union$area <- st_area(pa_Union)
sum(pa_Union$area) #29348526878 [m^2]

#get area of each MTB_Q
mtbqs$mtbq_area <- st_area(mtbqs)
sum(mtbqs$mtbq_area) #390641923430 [m^2]
#7.5%

#get area per MTBQ
MTBQ_pa  <-  st_intersection(mtbqs[,"MTB_Q"],st_make_valid(pa_Union))
MTBQ_pa$pa_area <- st_area(MTBQ_pa)
sum(duplicated(MTBQ_pa$MTB_Q))
sum(MTBQ_pa$pa_area)
#sum per MTBQ
MTBQ_pa <- MTBQ_pa %>%
            dplyr::group_by(MTB_Q) %>%
            dplyr::summarise(pa_area = sum(pa_area))


#combine
mtbqs$pa_area <- MTBQ_pa$pa_area[match(mtbqs$MTB_Q,MTBQ_pa$MTB_Q)]
mtbqs$pa_area[is.na(mtbqs$pa_area)] <- 0


#put into a data frame
protected_area <- data.frame(MTB_Q = mtbqs$MTB_Q, PA_area = as.numeric(mtbqs$pa_area/mtbqs$mtbq_area))
saveRDS(protected_area,file="environ-data/protectedarea_MTBQ.rds")

#old

# setwd("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/ProtectedAreas_DE/BFN")
# #create a temp file
# for(i in 1:length(paFiles)){
#   temp <- tempfile()
#   #unzip the contents and save unzipped content in 'temp'
#   unzip(zipfile = paFiles[i], exdir = temp)
#   #finds the filepath of the shapefile (.shp) file in the temp2 unzip folder
#   #the $ at the end of ".shp$" ensures you are not also finding files such as .shp.xml 
#   your_SHP_file<-list.files(temp, pattern = ".shp$",full.names=TRUE)
#   
#   #read the shapefile. 
#   paShape <- sf::read_sf(your_SHP_file)
#   
#   #transform to utm
#   paShape <- st_transform(paShape,crs=st_crs(mtbqs))
#   
#   #clean in
#   paShape <- st_buffer(paShape,dist=0)
#   #head(paShape)
#   #plot(as(paShape,'Spatial'))
#   assign(paste0("paShape",i), paShape)
# }     
# 
# #merge them
# object.size(paShape1)
# object.size(paShape2)#large
# object.size(paShape3)
# object.size(paShape4)
# object.size(paShape5)#large
# object.size(paShape6)
# object.size(paShape7)
# object.size(paShape8)#large
# object.size(paShape9)
# object.size(paShape10)
# object.size(paShape11)
# object.size(paShape12)#large
# object.size(paShape13)
# 
# paShape1 <- st_make_valid(paShape1)
# paShape3 <- st_make_valid(paShape3)
# temp <- st_union(paShape1,paShape3,by_feature=TRUE)
# 
# paShape4 <- st_make_valid(paShape4)
# temp <- st_union(temp,paShape4,by_feature=TRUE)
# 
# paShape6 <- st_make_valid(paShape6)
# temp <- st_union(temp,paShape6,by_feature=TRUE)
# 
# paShape7 <- st_make_valid(paShape7)
# temp <- st_union(temp,paShape7,by_feature=TRUE)
