#the script asks whether species about:
#range filler or range expander, and 
#change in number of occupied sites versus change in areas of convex hull

#https://help.natureserve.org/biotics/content/record_management/Element_Files/Element_Ranking/ERANK_Definitions_of_Extent_of_Occurrence_and_Area_of_Occupancy.htm
#The range extent will be overestimated by using a minimum convex polygon (also called a convex hull) to calculate the area. In these cases, the α-hull is recommended. The α-hull can be estimated by making a Delauney triangulation of the points in a sample (connect the points with lines, constrained so that no lines intersect), and then deleting lines that are longer than two times the average line length. The range extent is then the sum of enclosed areas. 
#https://www.ala.org.au/spatial-portal-help/aoo/

#Joppa, L. N., Butchart, S. H. M., Hoffmann, M., Bachman, S. P., Akçakaya, H. R., Moat, J. F., Böhm, M., #Holland, R. A., Newton, A., Polidoro, B. and Hughes, A. (2016), Impact of alternative metrics on #estimates of extent of occurrence for extinction risk assessment. Conservation Biology, 30: 362–370. doi:10.1111/cobi.12591

#dealing with outliers
#https://search.r-project.org/CRAN/refmans/CoordinateCleaner/html/cc_outl.html


### get data ####

#run the modelSummaries_stan

### subset data #####

#just to test below code
modelSummaries$PA <- sapply(modelSummaries$mean, function(x) rbinom(1,1,x))

#are we subsetting?
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

allspecies <- sort(unique(modelSummaries$Species))

allYears <- 1990:2016

### functions for analysis #####

getRangeArea <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  
  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(nuGrids = length(unique(MTB))) %>%
    dplyr::ungroup() %>%
    tidyr::complete(Year = allYears) %>%
    dplyr::mutate(nuGrids = ifelse(is.na(nuGrids),0,nuGrids))
  
}

getRangeExtents <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]

  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(max_Y = max(y_MTB), 
                     min_Y = min(y_MTB)) %>%
    dplyr::ungroup()
  
}

getConvexHull <- function(species, modelSummaries_Limits){
  
  require(sp)
  require(rgdal)
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
                   
    rangeHull <- sapply(allYears,function(x){
      
      ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
      if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
        ch <- chull(ydat)
        coords <- ydat[c(ch, ch[1]), ] 
        plot(ydat, pch=19, main = species)
        lines(coords, col="red")
        sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)),
                                 proj4string=CRS("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
        rangeSize <- raster::area(sp_poly)
        return(rangeSize)
    }else {
          return(0)
        }})
    
    data.frame(Species = species, Year = allYears, rangeHull = rangeHull) 
}


getAlphaHull <- function(species, modelSummaries_Limits){
  
  require(sp)
  require(rgdal)
  require(alphahull)
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeAlpha <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      
      dist <- max(ydat$y_MTB)-min(ydat$y_MTB)
      ah <- ahull(ydat, alpha = dist)
      plot(ah, main = species)
      areaahull(ah)
      
    }else {
      return(0)
    }})
    
    data.frame(Species = species, Year = allYears, rangeAlpha = rangeAlpha) 
  
}


getConcaveMan <- function(species, modelSummaries_Limits){
  
  require(concaveman)
  require(sf)
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeMan <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      
      ydat <- st_as_sf(ydat, coords =c("x_MTB", "y_MTB"),
                       crs = 25833)
      sp_poly <- concaveman(ydat)
      plot(ydat, main = species)
      plot(sp_poly,add=T)
      rangeSize <- as.numeric(st_area(sp_poly))#m2 units
      return(rangeSize)
      
    }else {
      return(0)
    }})
  
    data.frame(Species = species, Year = allYears, rangeMan = rangeMan)
    
}

### range area (AOO) ####

applyRangeArea <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeArea(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  
  #summarise change
  if(summary == "change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year,values_from=nuGrids) %>%
    janitor::clean_names() %>%  
    tidyr::complete() %>%
    dplyr::mutate(x1990 = ifelse(is.na(x1990),1,x1990)) %>%
    dplyr::mutate(x2016 = ifelse(is.na(x2016),1,x2016)) %>%
    dplyr::mutate(change = log((x2016)/(x1990)))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(myMedian = median(change), 
                     lower = quantile(change, 0.25),
                     upper = quantile(change, 0.75))
    
  }else if (summary == "annual") {
    
  #summarise annual
   temp %>%
     dplyr::group_by(Species,Year) %>%
     dplyr::summarise(myMedian = median(nuGrids), 
                      lower = quantile(nuGrids, 0.25),
                      upper = quantile(nuGrids, 0.75))
    
  }
  
}

### latitudinal extents ####

applyRangeExtent <- function(species, modelSummaries_Limits, summary="change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeExtents(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary =="change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year,values_from=c(max_Y,min_Y)) %>%
    janitor::clean_names() %>%  
    tidyr::complete() %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(change_maxY = (max_y_2016 - max_y_1990), change_minY = (min_y_2016 - min_y_1990))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(myMedian_Max = median(change_maxY), 
                     lower_Max = quantile(change_maxY, 0.25),
                     upper_Max = quantile(change_maxY, 0.75),
                     myMedian_Min = median(change_minY), 
                     lower_Min = quantile(change_minY, 0.25),
                     upper_Min = quantile(change_minY, 0.75))
    
  } else if(summary =="annual"){
    
  #summarise annual
   temp %>%
     dplyr::group_by(Species,Year) %>%
     dplyr::summarise(myMedian_Max = median(max_Y), 
                      lower_Max = quantile(max_Y, 0.25),
                      upper_Max = quantile(max_Y, 0.75),
                      myMedian_Min = median(min_Y), 
                      lower_Min = quantile(min_Y, 0.25),
                      upper_Min = quantile(min_Y, 0.75))
  }
  
}

### extent of occupancy (EOCC) ####

applyConcaveMan <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getConcaveMan(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeMan) %>%
    janitor::clean_names() %>%  
    dplyr::mutate(change = log((x2016+1)/(x1990+1)))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(myMedian = median(change), 
                     lower = quantile(change, 0.25),
                     upper = quantile(change, 0.75))
    
  } else if(summary == "annual"){
    
    #summarise annual
     temp %>%
       dplyr::group_by(Species,Year) %>%
       dplyr::summarise(myMedian = median(rangeMan), 
                        lower = quantile(rangeMan, 0.25),
                        upper = quantile(rangeMan, 0.75))
    
  }
  
}

### application ####

#get realizations
PAs <- lapply(modelSummaries$mean, function(x) rbinom(10,1,x))
PA_matrix <- do.call(rbind,PAs)

#area
areaChanges <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits)})
areaChanges  <- do.call(rbind,areaChanges)

#order species by range changes
areaChanges <- arrange(areaChanges, myMedian)
areaChanges$Species <- factor(areaChanges$species, levels=areaChanges$species)

#plot
ggplot(areaChanges)+
  geom_pointrange(aes(x = Species, y = myMedian, ymin = lower, max = upper))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_few()

#lat extent
rangeExtents <- lapply(allspecies, function(x){
  applyRangeExtent(x,modelSummaries_Limits)})
rangeExtents  <- do.call(rbind,rangeExtents)

#order species by range changes
rangeExtents <- arrange(rangeExtents, myMedian_Max)
rangeExtents$Species <- factor(rangeExtents$species, levels=rangeExtents$species)

#plot
ggplot(rangeExtents)+
  geom_pointrange(aes(x = Species, y = myMedian_Max, ymin = lower_Max, max = upper_Max),
                  colour="blue")+
  geom_pointrange(aes(x = Species, y = myMedian_Min, ymin = lower_Min, max = upper_Min),
                  colour="red")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_few()

#hull
hullChanges <- lapply(allspecies, function(x){
  applyConcaveMan(x,modelSummaries)})
hullChanges  <- do.call(rbind,hullChanges)

#order species by range changes
hullChanges <- arrange(hullChanges, myMedian)
hullChanges$Species <- factor(hullChanges$species, levels=hullChanges$species)

#plot
ggplot(hullChanges)+
  geom_pointrange(aes(x = Species, y = myMedian, ymin = lower, max = upper))+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_few()

### species ####

myspecies <- "Sympetrum danae"
myspecies <- "Crocothemis erythraea"

hullSpecies <- applyConcaveMan(myspecies,modelSummaries, summary="annual")
g1 <- ggplot(hullSpecies)+
  geom_pointrange(aes(x = Year, y = myMedian, ymin = lower, ymax = upper))+
  ylab("EOCC") + ggtitle(myspecies) +
  theme_few()

areaSpecies <- applyRangeArea(myspecies,modelSummaries, summary="annual")
mult <- 31000
g2 <- ggplot(areaSpecies)+
  geom_pointrange(aes(x = Year, y = myMedian *mult, ymin = lower*mult, ymax = upper*mult))+
  ylab("AOCC") + ggtitle(myspecies) +
  theme_few()

cowplot::plot_grid(g1,g2)

#ratio
hullSpecies$packing <- hullSpecies$myMedian/areaSpecies$myMedian
ggplot(hullSpecies)+
  geom_line(aes(x = Year, y = packing))+
  ylab("Packing") + ggtitle(myspecies) +
  theme_few()
