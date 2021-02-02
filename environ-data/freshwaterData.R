
#get population data
library(plyr)
############################################################################################

#read in dragonfly dataset

adultData <- readRDS("derived-data/adultData_allStates_Dec2020.rds")

##########################################################################################

#format date
library(lubridate)
adultData$Date<-as.Date(adultData$Date)
adultData$Year <- year(adultData$Date)
adultData$yday <- yday(adultData$Date)
adultData$week <- week(adultData$Date)
adultData$Month <- month(adultData$Date)
adultData$Day <- day(adultData$Date)
adultData <- subset(adultData,!is.na(Date))
adultData <- subset(adultData,yday!=1)

df <- subset(adultData, Year>=1975)

##########################################################################################

#get german county boundaries
germanAdmin <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/AdminBoundaries/gadm36_DEU_1_sp.rds")

#get MTBQs
library(rgdal)
mtbqs <- readOGR(dsn="C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/MTB_Q Informations/MTBQ_shapefile",
                 layer="MTBQ_25833")

##################################################################################

#read in data from ralf schafer
tdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/WaterQuality/Data_monitoring/Data_monitoring"
samples <- read.csv(paste(tdir,"samples.csv",sep="/"))
sites <- read.csv(paste(tdir,"sites.csv",sep="/"))
head(samples)
head(sites)

#plot the sites
library(sp)
coordinates(sites) <- c("easting","northing")
proj4string(sites) <- CRS("+init=epsg:31467")
#sites <- spTransform(sites, CRS(proj4string(germanAdmin)))
#plot(sites,add=T)

#overlay the coordinates with the mtqb map
sites <- spTransform(sites, CRS(proj4string(mtbqs)))
plot(mtbqs)
plot(sites,add=T,col="red")

#get the mtbq of each site
sitesMTB_Q <- over(sites,mtbqs)
sitesMTB_Q$Q <- NA
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="NW")]<-1
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="NO")]<-2
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="SW")]<-3
sitesMTB_Q$Q[which(sitesMTB_Q$Quadrant=="SO")]<-4
sitesMTB_Q$MTB_Q <- paste0(as.character(sitesMTB_Q$Value),
                        as.character(sitesMTB_Q$Q))
head(sitesMTB_Q)

#########################################################################################

#add the MTBQ to the site file
sites <- read.csv(paste(tdir,"sites.csv",sep="/"))
sites$MTB_Q <- sitesMTB_Q$MTB_Q
samples$MTB_Q <- sites$MTB_Q[match(samples$site_id,sites$site_id)]

#merge with samples
library(lubridate)
samples$Date <- as.Date(samples$date)
samples$Year <- year(samples$Date)
samples$yday <- yday(samples$Date)
samples$week <- week(samples$Date)
samples$Month <- month(samples$Date)
samples <- subset(samples,!is.na(date))
samples <- subset(samples,Year>=1975 & Year<2019)

summary(samples$Year)

########################################################################################

#get info on which MTBS have water qualities samples (least 5 years)
waterSummary <- ddply(samples,.(MTB_Q),summarise,
                      nuYears=length(unique(Year)),
                      TimeSpan=max(Year)-min(Year))
waterSummary <- subset(waterSummary, nuYears>=5 & TimeSpan>=10)

#get info on which MTBS
dfSummary <- ddply(df, .(MTB_Q),summarise,
                   nuYears=length(unique(Year)),
                   TimeSpan=max(Year)-min(Year))
dfSummary <- subset(dfSummary, nuYears>=5 & TimeSpan>=10)

#get samples overlapping
samples <- subset(samples, MTB_Q %in% waterSummary$MTB_Q)
samples <- subset(samples, MTB_Q %in% dfSummary$MTB_Q)
nrow(samples)

#nrow(samples)
#[1] 7658027
# length(unique(samples$MTB_Q))
#[1] 1014
write.table(samples,file="samplesforDragonflies.txt",sep="\t",row.names=FALSE)

########################################################################################

#examine water quality patterns
tdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Spatial_data/WaterQuality/Data_monitoring/Data_monitoring"
sites <- read.csv(paste(tdir,"psm_sites.csv",sep="/"),sep=";")
sites_info <- read.csv(paste(tdir,"psm_sites_info.csv",sep="/"),sep=";")
maxtu <- read.csv(paste(tdir,"psm_maxtu.csv",sep="/"),sep=";")

head(sites)
head(sites_info)
head(maxtu)

#samples extracted
load(paste(tdir,"Samples_extracted.RData",sep="/"))
str(samples_fin_sub)
head(samples_fin_sub)

#add year
library(lubridate)
samples_fin_sub$date <- as.Date(samples_fin_sub$date)
samples_fin_sub$Year <- year(samples_fin_sub$date)

#add site data


#look at variable names
unique(samples_fin_sub$name)

#plot histograms
library(ggplot2)
ggplot(data=subset(samples_fin_sub, name %in% unique(samples_fin_sub$name)[c(15:31)]))+
  geom_histogram(aes(value))+
  facet_wrap(~name,scales="free",ncol=5)
  
#which are the most commonly sampled variables
library(plyr)
varSummary <- ddply(samples_fin_sub,.(name),summarise,
                    nuDates = length(unique(date[!is.na(value)])),
                    nuSites = length(unique(site_id[!is.na(value)])),
                    nuYears = length(unique(Year[!is.na(value)])))
varSummary
  
  
#example correlations among variables
library(reshape2)
dataCast <- dcast(samples_fin_sub,sample_id+site_id+date~name,value.var="value")
cors <- cor(dataCast[,4:ncol(dataCast)],use="pairwise.complete.obs")
corsM <- melt(cors)
corsM[corsM==1] <- NA
library(corrplot)
corrplot(cors,method = "color",tl.col = "black",tl.cex=0.65)

subset(corsM,abs(value)>0.7)

ggplot(corsM)+
  geom_tile(aes(x=Var1,y=Var2,fill=value))+
  theme

#PCA
df <- dataCast[,4:ncol(dataCast)]
df = kNN(df)

names(df) <- gsub(" ","",names(df))
colnames(df) <- paste("var", 1:31, sep="")
princomp(~., data = df, 
         cor = TRUE, na.action=na.exclude)

summary(fit) 
loadings(fit) 
plot(fit,type="lines") 
fit$scores 
biplot(fit)

# impute missing values
library(missMDA)
# estimate number of components
nb <- estim_ncpPCA(df, ncp.min=0, ncp.max=5)
# actual impute
rr.impute <- imputePCA(df, ncp=2)

# Run pca
pca.fit <- prcomp(rr.impute$completeObs,scale=TRUE)

library(factoextra)
fviz_pca_var(pca.fit,repel = TRUE)

#time-series

#mean means per year
vars <- c("Wassertemperatur","Ammonium-Stickstoff","Nitrat","Sauerstoffgehalt")
summaryYear <- ddply(samples_fin_sub,.(name,Year,site_id),summarise,meanVal = median(value,na.rm=T))

ggplot(subset(summaryYear,name %in% vars),aes(x=Year,y=meanVal))+
        geom_line(aes(group=site_id))+
        facet_wrap(~name,scales="free")

ggplot(subset(summaryYear,name %in% vars),aes(x=Year,y=meanVal))+
  stat_smooth()+
  facet_wrap(~name,scales="free")


quantilesSummary <- ddply(samples_fin_sub,.(name),summarise,
                          min=min(value,na.rm=T),
                          lowerQ=quantile(value,0.25,na.rm=T),
                          median=quantile(value,0.5,na.rm=T),
                          mean=mean(value,na.rm=T),
                          upperQ=quantile(value,0.75,na.rm=T),
                          max=max(value,na.rm=T))


#get trend at each site_id
trends <- ddply(summaryYear,.(name,site_id),function(x){
  
  if(length(unique(x$Year))>10){
    lm1 <- lm(meanVal ~ Year, data=x)
    temp <- data.frame(t(summary(lm1)$coef[2,]))
    temp$minYear <- min(x$Year)
    temp$maxYear <- max(x$Year)
    return(temp)
  }
  
})

#merge with site data
trends <- merge(trends,sites,by="site_id")

ggplot(subset(trends,name==vars[4]),
       aes(x=easting, y=northing))+
  geom_point(aes(colour=Estimate))+
  scale_colour_gradient2(low="red",mid="white",high="blue",midpoint=0)+
  ggtitle(vars[4])

########################################################################################