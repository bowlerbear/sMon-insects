#run the modelSummaries_stan

#the script asks whether species about:
#range filler or range expander, and 
#change in number of occupied sites versus change in areas of convex hull

#https://help.natureserve.org/biotics/content/record_management/Element_Files/Element_Ranking/ERANK_Definitions_of_Extent_of_Occurrence_and_Area_of_Occupancy.htm
#The range extent will be overestimated by using a minimum convex polygon (also called a convex hull) to calculate the area. In these cases, the α-hull is recommended. The α-hull can be estimated by making a Delauney triangulation of the points in a sample (connect the points with lines, constrained so that no lines intersect), and then deleting lines that are longer than two times the average line length. The range extent is then the sum of enclosed areas. 
#https://www.ala.org.au/spatial-portal-help/aoo/


#function to get get range area and extent for each species
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))
allspecies <- sort(unique(modelSummaries_Limits$Species))

### functions for analysis #####

#take from https://babichmorrowc.github.io/post/2019-03-18-alpha-hull/
ashape2poly <- function(ashape){
  # Convert node numbers into characters
  ashape$edges[,1] <- as.character(ashape$edges[,1])
  ashape_graph <- graph_from_edgelist(ashape$edges[,1:2], directed = FALSE)
  if (!is.connected(ashape_graph)) {
    stop("Graph not connected")
  }
  if (any(degree(ashape_graph) != 2)) {
    stop("Graph not circular")
  }
  if (clusters(ashape_graph)$no > 1) {
    stop("Graph composed of more than one circle")
  }
  # Delete one edge to create a chain
  cut_graph <- ashape_graph - E(ashape_graph)[1]
  # Find chain end points
  ends = names(which(degree(cut_graph) == 1))
  path = get.shortest.paths(cut_graph, ends[1], ends[2])$vpath[[1]]
  # this is an index into the points
  pathX = as.numeric(V(ashape_graph)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  return(pathX)
}

getRangeArea <- function(myspecies){
  
  speciesData <- subset(modelSummaries_Limits,Species==myspecies)
  
  dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  years <- names(tapply(dat$mean,dat$Year,sum))
  
  nuGrids <- as.numeric(tapply(dat$mean,dat$Year,sum))
  
  data.frame(Species = myspecies, Year = years, nuGrids)
  
}

getRangeExtents <- function(myspecies){
  
  speciesData <- subset(modelSummaries_Limits,Species==myspecies)
  
  dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  years <- names(tapply(dat$x_MTB,dat$Year,max))
  max_X <- as.numeric(tapply(dat$x_MTB,dat$Year,max))
  min_X <- as.numeric(tapply(dat$x_MTB,dat$Year,min))
  max_Y <- as.numeric(tapply(dat$y_MTB,dat$Year,max))
  min_Y <- as.numeric(tapply(dat$y_MTB,dat$Year,min))
  
  data.frame(Species = myspecies, Year = years, max_X, min_X, max_Y, min_Y)
  
}

plotAlphaHull <- function(myspecies){
  
  require(sp)
  require(rgdal)
  require(alphahull)
  
  speciesData <- subset(modelSummaries_Limits,Species==myspecies)
  
  #put hull around data for each year
  dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  #present in at least 10 MTBs in each year
  nuYears <- length(unique(dat$Year))
  nuMTBs <- tapply(dat$MTB,dat$Year,function(x)length(unique(x)))
    
  if(nuYears==2 & all(nuMTBs>10)){
    
  myYears <- sort(unique(speciesData$Year))
  
  par(mfrow=c(1,2))

  #start year  
  ydat <- dat[dat$Year==min(myYears),c("x_MTB","y_MTB")]
  ah <- ashape(ydat, alpha = 20000)
  plot(ah, main = myspecies)
  
  #end year
  ydat <- dat[dat$Year==max(myYears),c("x_MTB","y_MTB")]
  ah <- ashape(ydat, alpha = 20000)
  plot(ah, main = myspecies)
  
  }
  
}


makeAlphaShapes <- function(myspecies){
  ashape2poly(ah)
}

getConvexArea <- function(species){
  
  require(sp)
  require(rgdal)
  
  speciesData <- subset(modelSummaries_Limits,Species==myspecies)
  
  #put hull around data for each year
  dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  #present in at least 10 MTBs in each year
  myYears <- sort(unique(dat$Year))
  nuYears <- length(myYears)
  nuMTBs <- tapply(dat$MTB,dat$Year,function(x)length(unique(x)))
  
  if(nuYears==2 & all(nuMTBs>2)){
    
    rangeHull <- sapply(myYears,function(x){
      ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
        ch <- chull(ydat)
        coords <- ydat[c(ch, ch[1]), ] 
        #plot(ydat, pch=19)
        #lines(coords, col="red")
        sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)),
                                   proj4string=CRS("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
        rangeSize <- raster::area(sp_poly)
        rangeSize <- rangeSize/31000#area in m2 of MTBQ 
        return(rangeSize)
    })
   data.frame(Species = myspecies, Year = myYears, rangeHull = rangeHull) 
  }
}

### applying functions #####

for(species in allspecies){
  plotAlphaHull(species)
}




