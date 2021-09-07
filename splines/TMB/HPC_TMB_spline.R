#based on 
#https://github.com/r-glennie/occuR
#https://r-glennie.github.io/occuR/
#remotes::install_github("r-glennie/occuR")

#correlation of site occupancy over space and time is induced by allowing occupancy probability to be a smooth function of space and time

#library(occuR)
library(occuR, lib.loc="/gpfs0/home/bowler")
library(tidyverse)

myfolder <- "/data/idiv_ess/Odonata"#HPC

### load in the full dataset ####

adultData <- readRDS(paste(myfolder,"adultData.rds",sep="/"))
load(paste(myfolder,"mtbqsDF.RData",sep="/"))
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

# subset  

df <- subset(adultData, Year>=1990  & Year<2017)
#table(df$Month)
df <- subset(df, Month %in% 4:10)

# get species

speciesTaskID <- read.delim(paste(myfolder,"speciesTaskID_adult.txt",sep="/"),as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)]
#myspecies <- "Sympetrum danae" #test case

#get MTB
df$MTB <- sapply(as.character(df$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,(len-1))})

### function to get data for a species ####

all(df$MTB %in% mtbqsDF$MTB)

df$State <- mtbsDF$Counties[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$State))

#check missing naturruam data
mtbsDF$Natur[is.na(mtbsDF$Natur)] <- mtbsDF$MTB_Natur[is.na(mtbsDF$Natur)]
df$Natur <- mtbsDF$Natur[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$Natur))

mtbsDF$CoarseNatur[is.na(mtbsDF$CoarseNatur)] <- mtbsDF$MTB_CoarseNatur[is.na(mtbsDF$CoarseNatur)]
df$CoarseNatur <- mtbsDF$CoarseNatur[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$CoarseNatur))

mtbsDF$MTB_MidNatur <- gdata::trim(mtbsDF$MTB_MidNatur)
df$MidNatur <- mtbsDF$MTB_MidNatur[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$MidNatur))

#define a visit
df$visit <- paste(df$MTB,df$Date,df$Beobachter,sep="_")

#get occurence matrix  - detection and non-detection
getOccurrenceMatrix<-function(df){
  out <- reshape2::acast(df,visit~Species,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
  out[out>0]<-1
  return(out)
}
occMatrix <- getOccurrenceMatrix(df)

#get list length
getListLength<-function(df){
  df %>% 
    group_by(visit,Date,MTB) %>%
    summarise(nuSpecies=length(unique(Species)),
              nuRecords=length(Species)) %>%
    arrange(.,visit)
}
listlengthDF <- getListLength(df)

#check rows of occuMatrix match visits
all(listlengthDF$visit==row.names(occMatrix))

#add on some indices
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- lubridate::year(listlengthDF$Date)
listlengthDF$yday <- lubridate::yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))

#add natur raum
listlengthDF$CoarseNaturraum <- mtbqsDF$CoarseNatur[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
subset(listlengthDF,is.na(CoarseNaturraum))

listlengthDF$Naturraum <- mtbqsDF$Natur[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$nnIndex <- as.numeric(factor(listlengthDF$Naturraum))
subset(listlengthDF,is.na(Naturraum))

listlengthDF$MidNaturraum <- mtbqsDF$MTB_MidNatur[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$mnIndex <- as.numeric(factor(listlengthDF$MidNaturraum))
subset(listlengthDF,is.na(MidNaturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)

#add on species data to the listlength object
listlengthDF$Species <- occMatrix[,myspecies]
table(listlengthDF$Species)

#rename
visit_data <- listlengthDF

#add coordinates
visit_data$x <- mtbsDF$x_MTB[match(visit_data$MTB, mtbsDF$MTB)]
visit_data$y <- mtbsDF$y_MTB[match(visit_data$MTB, mtbsDF$MTB)]

#make coordinates smaller
visit_data$x <- visit_data$x/10000
visit_data$y <- visit_data$y/1000000

#multi-occasion occupancy

#need to make visit data with columns: site, occassion and obs
visit_data <- visit_data[,c("MTB","Year","visit","Species",
                            "singleList","yday","CoarseNaturraum","x","y")]

names(visit_data)[1:4] <- c("site","occasion","visit","obs")

#need to make visit be indexed from i to n within each site and occasion
visit_data <- visit_data %>%
              group_by(site, occasion) %>%
              mutate(visit = as.numeric(as.factor(visit)))%>%
              ungroup()

visit_data$occasion <- as.numeric(as.factor(visit_data$occasion))

#need to make site data with "site" and "occasion"
site_data <- unique(visit_data[,c("site","occasion","CoarseNaturraum","x","y")])

#basic model
m0 <- fit_occu(list(psi ~ 1, p ~ 1), as.data.table(visit_data), as.data.table(site_data))
m0

# #two dimension
# m_spline2d <- fit_occu(list(psi ~ t2(x,y,bs = "ts", k=10), p ~ 1),
#                      as.data.table(visit_data), as.data.table(site_data))
# m_spline2d
# #k=5 worked pretty well
# #k=10 more wiggly.
# #k=15 strange lines
# 
# #predictions
# pred_xy <- predict(m_spline2d,
#                    as.data.table(visit_data),
#                    data.table(occasion = 1, x = mtbsDF$x_MTB/10000, y = mtbsDF$y_MTB/1000000),
#                    nboot = 1000)
# 
# summary(pred_xy$psi)
# mtbsDF$preds <- pred_xy$psi[,1]
# 
# ggplot(mtbsDF) +
#   geom_point(aes(x = x_MTB, y = y_MTB, colour = preds)) +
#   theme_bw() +
#   scale_colour_viridis_c("Occupancy")


# spatio-temporal effect

m_spline3d <- fit_occu(list(psi ~ t2(x, y, occasion, bs = c("ts", "cs"), k=c(5,2)), p ~ 1),
                       as.data.table(visit_data), as.data.table(site_data))

# Dataset for predictions
siteInfo_NAs <- mtbsDF
nuSites <- nrow(siteInfo_NAs)
siteInfo_NAs <- siteInfo_NAs %>% slice(rep(1:n(), each = length(unique(df$Year))))
siteInfo_NAs$Year <- rep(sort(unique(df$Year)),nuSites)
siteInfo_NAs <- siteInfo_NAs[,c("Year","MTB","MTB_CoarseNatur","x_MTB","y_MTB")]
siteInfo_NAs$occassion <- as.numeric(as.factor(siteInfo_NAs$Year))

pred_xyt <- predict(m_spline3d, 
                    as.data.table(visit_data), 
                    data.table(occasion = siteInfo_NAs$occassion, 
                               x = siteInfo_NAs$x_MTB/10000, 
                               y = siteInfo_NAs$y_MTB/1000000), 
                    nboot = 1000)

#save predicted values
siteInfo_NAs$preds <- pred_xyt$psi[,1]
saveRDS(siteInfo_NAs,file=paste0("TMB_spline_",myspecies,".rds"))