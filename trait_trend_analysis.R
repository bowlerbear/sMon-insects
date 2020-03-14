#linking traits to trends

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/sparta_wrapper_functions.R')
library(plyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(wesanderson)

###merge data###################################################################

setwd("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects")

#get traits data
load("alltraits.RData")
#limit to those with complete cases?

#get trends data
trendEstimates <- readRDS("model-outputs/modelTrends_trends.rds")
trendEstimates$Species <- gsub("out_dynamic_nation_naturraum_adult_","",trendEstimates$file)
trendEstimates$Species <- gsub(".rds","",trendEstimates$Species)
names(trendEstimates)[which(names(trendEstimates)=="mean")]<-"trend"

#merge in one data frame
trendEstimates$Species[!trendEstimates$Species %in% alltraits$Species]
trendEstimates <- merge(trendEstimates,alltraits,by="Species")
nrow(trendEstimates)

###trend summary##############################################################

#add trend
trendEstimates$Trend <- "insignificant"
trendEstimates$Trend[trendEstimates$X2.5.>0 & trendEstimates$X97.5.>0]<-"significant increase"
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.<0]<-"significant decrease"
table(trendEstimates$Trend)

summary(trendEstimates$trend[trendEstimates$Trend=="significant increase"])
#1.0034478^35 

summary(trendEstimates$trend[trendEstimates$Trend=="significant decrease"])
#1-(1-0.003953)^35

###taxonomy###################################################################################

trendEstimates$Suborder[trendEstimates$Species=="Oxygastra curtisii"] <- "Anisoptera"

table(trendEstimates$Suborder)

table(trendEstimates$Suborder,trendEstimates$Trend)

chisq.test(table(trendEstimates$Suborder,trendEstimates$Trend))
#ns

###linear models######################################################################

#use dev time or voltinim
summary(lm(trend~VoltinismProp,data=trendEstimates,
           weights=1/sd))#use this one
summary(lm(trend~devTime,data=trendEstimates,
           weights=1/sd))

#german of europe range size
summary(lm(trend~germanRange,data=trendEstimates,
           weights=1/sd))
summary(lm(trend~nuEuroGrids,data=trendEstimates,
           weights=1/sd))#use thie one
  

#build multiple regression model
library(arm)
lm1 <- lm(trend ~ rescale(VoltinismProp) + rescale(nuEuroGrids) + 
                     rescale(Flight_start) + rescale(meanTemp) + 
                    rescale(medHw) + 
                     Habitat + HabitatBreadth,data=trendEstimates,
                                   weights=1/sd)
summary(lm1)

#extract and plot coefficients
coefDF <- data.frame(cbind(summary(lm1)$coefficients,confint(lm1)))
coefDF$Param <- c("int","voltinism","range size","flight start",
                  "temp pref","wing length","habitat","habitat breadth")

g1 <- ggplot(subset(coefDF,Param!="int"))+
  geom_crossbar(aes(x=Param,y=Estimate,ymin=X2.5..,ymax=X97.5..),width=0.3)+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  ylab("Effect on trend")+xlab("Trait")

lm1 <- lm(trend ~ rescale(nuEuroGrids) + 
            rescale(meanTemp) + 
            rescale(Flight_start) + 
            rescale(medHw),data=trendEstimates,
          weights=1/sd)

###plots#################################################################

#plot trait vs trends

q1<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=TMean,y=trend,colour=State))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("Temperature preference")+
  facet_grid(~State)+stat_smooth(method="lm")+
  theme(legend.position="none")
q2<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=Habitat.y,y=trend,fill=State))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylab("Population trend")+xlab("Habitat preference")+
  facet_grid(~State)+
  theme(legend.position="none")
q3<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=total,y=trend,colour=State))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("Range size (Germany)")+
  facet_grid(~State)+stat_smooth(method="lm")+
  theme(legend.position="none")
q4<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=nuEuroGrids,y=trend,colour=State))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("Range size (Europe)")+
  facet_grid(~State)+stat_smooth(method="lm")+
  theme(legend.position="none")
q5<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=development.time..mean.,y=trend,colour=State))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("Development time")+
  facet_grid(~State)+stat_smooth(method="lm")+
  theme(legend.position="none")
q6<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=wing_size,y=trend,colour=State))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("Wing length")+
  facet_grid(~State)+stat_smooth(method="lm")+
  theme(legend.position="none")

ggplot(trendEstimates,aes(x=SZP^2,y=trend))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("SZP")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

library(cowplot)
plot_grid(q1,q2,q3,q4,q5,q6,align="v",ncol=1)


#temp and habitat preference
q1<-ggplot(trendEstimates,aes(x=meanTemp,y=trend))+
  geom_point()+
  theme_bw()+
  geom_hline(yintercept=0,colour="red",linetype="dashed")+
  ylab("Population trend")+xlab("Temperature preference")+
  stat_smooth(method="lm",se=F)+
  theme(legend.position="none")

q2<-ggplot(trendEstimates,aes(x=medHw,y=trend))+
  geom_point()+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_hline(yintercept=0,colour="red",linetype="dashed")+
  ylab("Population trend")+xlab("Wing length")+
  theme(legend.position="none")

q3<-ggplot(trendEstimates,aes(x=Flight_start,y=trend))+
  geom_point()+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_hline(yintercept=0,colour="red",linetype="dashed")+
  ylab("Population trend")+xlab("Start of flight period")+
  theme(legend.position="none")

q4 <- ggplot(trendEstimates,aes(x=nuEuroGrids,y=trend))+
  geom_point()+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_hline(yintercept=0,colour="red",linetype="dashed")+
  ylab("Population trend")+xlab("European range size")+
  theme(legend.position="none")

library(cowplot)
plot_grid(q1,q2,q3,q4,align="h",nrow=1)

####CWM#################################################################

#total average

load("randomMatrix.RData")
library(reshape2)
randomMatrixM<-melt(randomMatrix,id=c("Year","Species"))
randomMatrixM<-arrange(randomMatrixM,Species,Year)
randomMatrixM$Species[!randomMatrixM$Species %in% alltraits$Species]
randomTrends <- merge(randomMatrixM,alltraits,by="Species")

#wing length
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(medHw,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g1<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Wing length")+
  theme(legend.position="none")

#voltinism
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(VoltinismProp,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g2<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Voltinism")+
  theme(legend.position="none")


#temp sd
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(sdTemp,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g3<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Temp range")+
  theme(legend.position="none")

#tmean
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(meanTemp,value))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g4<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Temp pref")+
  theme(legend.position="none")


#european distribution
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(nuEuroGrids,value))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g5<-ggplot(subset(tmeansMeans,Year>1980))+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Range size")+
  theme(legend.position="none")

#habitat preference
randomTrends$Habitat.y <- ifelse(randomTrends$Habitat=="Standing",0,1)
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(Habitat.y,value))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g6<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Running water")+
  theme(legend.position="none")


#flight period start
#check missing value
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(Flight_start,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g7<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Start of flight period")+
  theme(legend.position="none")

#SZP
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(SZP,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

g8<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("SZP")+
  theme(legend.position="none")

plot_grid(g4,g1,g7,g5,nrow=1)

####species richness#################################################

load("randomMatrix.RData")

#library(reshape2)
#randomMatrixM<-melt(randomMatrix,id=c("Year","Species"))
#randomMatrixM<-arrange(randomMatrixM,Species,Year)

out <- ddply(randomMatrix,.(Year),function(x){
  numcolwise(sum)(x)})

#get mean and 95% CI across communities
out$meanRichness<-apply(out[,2:1001],1,median)
out$lowerRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.025))
out$upperRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.975))

g3 <- ggplot(subset(out,Year>1980))+
  geom_line(aes(x=Year,y=meanRichness))+
  geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness),alpha=0.3)+
  theme_bw()+
  ylab("Average species richness")


library(cowplot)
plot_grid(g1,g2,g3,nrow=1)

###disimmilarity##############################################

load("randomMatrix.RData")

library(vegan)
randomMatrixM<-melt(randomMatrix,id=c("Year","Species"))
randomMatrixM<-arrange(randomMatrixM,Species,Year)

nuSites <- 500
randomMatrixM$predProp <- round(randomMatrixM$value*nuSites)

out <- ddply(randomMatrixM,.(variable),function(x){
  moo2 <- acast(x,Year~Species,value.var="predProp")
  diss <- as.numeric(as.matrix(vegdist(moo2,method="bray",binary=FALSE))[,1])
  year <- names(as.matrix(vegdist(moo2,binary=TRUE))[,1])
  data.frame(year,diss)
})

out$year <- as.numeric(as.character(out$year))
out <- subset(out,year>1980)

#take average across all
out <- ddply(out,.(year),summarise,
             meanD = mean(diss),
             lowerCI = quantile(diss,0.025),
             upperCI = quantile(diss,0.975))

out$Year <-as.numeric(as.character(out$year))
g2 <- ggplot(out)+
  geom_line(aes(x=Year,y=meanD))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Average dissimilarity")

###diversity#######################################################################

data(BCI)
head(BCI)
dim(BCI)
diversity(BCI)

H <- diversity(acast(BCI,Year~Species,value.var="predProp"))

out <- ddply(randomMatrixM,.(variable),function(x){
  moo2 <- acast(x,Year~Species,value.var="predProp")
  div <- diversity(moo2) 
  year <- sort(unique(x$Year))
  data.frame(year,div)
})

out$year <- as.numeric(as.character(out$year))
out <- subset(out,year>1980)

#take average across all
out <- ddply(out,.(year),summarise,
             meanD = mean(div),
             lowerCI = quantile(div,0.025),
             upperCI = quantile(div,0.975))

out$Year <-as.numeric(as.character(out$year))

g3 <- ggplot(out)+
  geom_line(aes(x=Year,y=meanD))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Average diversity")

plot_grid(g1,g3,g2,nrow=1)

####clustering############################################################################################

#convert time series into a list
annualDFS <- subset(annualDF,Year>1979)
myTS <- dlply(annualDFS,.(Species),
              function(x){x[,"mean"]})

#preprocessing - mean
myTS <- lapply(myTS,function(x){
  x/mean(x)
})

#preprocessing - year 1
myTS <- lapply(myTS,function(x){
  x/x[1]
})

#combine and plot again
temp <- ldply(myTS,function(x){
  data.frame(Index=x,Year=1:length(x))
})
temp$Year <- temp$Year + min(annualDFS$Year)-1
temp$Species <-temp$.id
temp2 <- temp
qplot(Year,Index,data=temp,colour=Species)+
  theme(legend.position = "none")

annualDFS <- subset(annualDF,Year>1979)

###dwtclust#######################################################

library("dtwclust")

#SBD
hc_sbd <- tsclust(myTS, type = "partitional", k=2:20L,
                  preproc = zscore, 
                  distance = "sbd",centroid = "shape")
#DTW
hc_sbd <- tsclust(myTS, ype="partitional", k= 2:20L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "dba")
#pam
hc_sbd <- tsclust(myTS, ype="partitional", k= 2:20L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "pam")

#comparing cluster numbers
names(hc_sbd) <- paste0("k_", 2L:20L)
temp <- sapply(hc_sbd, cvi, type = "internal")
#plot each one
temp <- data.frame(temp)
temp$Param <- row.names(temp)
library(reshape2)
temp <- melt(temp,id="Param")
temp$variable <- as.numeric(gsub("k_","",temp$variable))
qplot(variable,value,data=temp) + facet_wrap(~Param,scales="free")

#get inflexion point for each plot
#inflexion point is maximum absolute second derivative
deriv <- function(x, y) diff(y) / diff(x)
middle_pts <- function(x) x[-1] - diff(x) / 2
out <- ldply(unique(temp$Param),function(x){
  temp2 <- subset(temp,Param==x)
  firstderiv <- deriv(temp2$variable, temp2$value)
  plot(firstderiv ~ temp2$variable[-1])
  second_d <- deriv(middle_pts(temp2$variable), firstderiv)
  plot(second_d ~ temp2$variable[-c(1:2)])
  df <- data.frame(second_d = second_d, 
           midpts = middle_pts(middle_pts(temp2$variable)),
           Param=x)
  subset(df,second_d==max(second_d))
})
median(out$midpts)#5

###pick K#####################################################################

#SBD
hc_sbd <- tsclust(myTS, type = "partitional", k=4L,
                  preproc = zscore, 
                  distance = "sbd",centroid = "shape")

#DTW
hc_sbd <- tsclust(myTS, ype="partitional", k= 3L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "dba")

#pam
hc_sbd <- tsclust(myTS, ype="partitional", k= 9L,
                  preproc = zscore, 
                  distance = "dtw", centroid = "pam")

plot(hc_sbd)
plot(hc_sbd, type = "sc")
#plot(hc_sbd, type = "series", clus = 1L)
plot(hc_sbd, type = "series")
plot(hc_sbd, type = "centroids")

#check
clusterDF <- data.frame(Species=names(myTS),
                        cluster=hc_sbd@cluster)

annualDFS$cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]
plotCluster(annualDFS)

###plot each cluster#######################################################

#order cluster by numbers of species
clusts <- table(clusterDF$cluster)
clusterOrder <- rev(names(clusts)[order(clusts)])

#examine cluster by cluster
ggplot(subset(annualDFS,cluster==1))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")

ggplot(subset(annualDFS,cluster==2))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")

ggplot(subset(annualDFS,cluster==3))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")

ggplot(subset(annualDFS,cluster==4))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")

ggplot(subset(annualDFS,cluster==5))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")

###plot clusters##############################################################

myCentroids <- data.frame(Year=rep(sort(unique(annualDFS$Year)),length(hc_sbd@centroids)),
                          Cluster=rep(1:length(hc_sbd@centroids),each=length(unique(annualDFS$Year))),
                          ts=do.call(c,hc_sbd@centroids))
myCentroids$Cluster <- factor(myCentroids$Cluster,levels=1:length(clusterOrder))

#smooth predicted series
ggplot(data=myCentroids,aes(x=Year,y=ts))+
  geom_smooth(aes(colour=factor(Cluster),fill=factor(Cluster)))+
  facet_wrap(~Cluster,nrow=1)+
  theme_bw()+
  theme(legend.position = "none")

#smooth underlying dynamics
temp2$Cluster <- clusterDF$cluster[match(temp2$Species,clusterDF$Species)]
ggplot(data=temp2)+
  geom_smooth(aes(x=Year,y=Index),size=rel(2))+
  facet_wrap(~Cluster,ncol=1)+
  theme_bw()+
  theme(legend.position = "none")

####bootstrap#############################################################

#bootstrap original values within each cluster at each step

#get year 1 mean value

#randomly pick value for each species

#then scale subsequent values by this

###final plots##################################################################

myOrder <- c(2,1,4,3)
mylabels <- c("increasing","increasing-decreasing",
              "decreasing-increasing","decreasing")

mycols <- wes_palette("Darjeeling1", length(myOrder))
myCentroids$Cluster <- factor(myCentroids$Cluster,levels=rev(myOrder))
levels(myCentroids$Cluster) <- rev(mylabels)

#plot as a smooth
ggplot(data=myCentroids,aes(x=Year,y=ts))+
  geom_smooth(aes(colour=factor(Cluster),fill=factor(Cluster)))+
  facet_wrap(~Cluster,nrow=1)+
  theme_bw()+
  scale_fill_manual(values=rev(mycols))+
  scale_colour_manual(values=rev(mycols))+
  theme(legend.position = "none")+
  ylab("relative occupancy prop")+
  theme(axis.title = element_text(size=rel(1.2)),
        axis.text = element_text(size=rel(1.2)),
        strip.text = element_text(size=rel(1.2)))

table(clusterDF$cluster)
#8 35 16 18 

#smooth underlying dynamics
temp2$Cluster <- clusterDF$cluster[match(temp2$Species,clusterDF$Species)]
myOrder <- c(1,2,3,4)
mylabels <- c("increasing","increasing-decreasing",
              "decreasing-increasing","decreasing")

mycols <- wes_palette("Darjeeling1", length(myOrder))
temp2$Cluster <- factor(temp2$Cluster,levels=myOrder)
levels(temp2$Cluster) <- mylabels

ggplot(data=temp2,aes(x=Year,y=Index))+
  geom_smooth(aes(colour=Cluster,fill=Cluster))+
  facet_wrap(~Cluster,ncol=1)+
  theme_bw()+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  theme(legend.position = "none")+
  ylab("relative occupancy prop")

###traits clusters########################################################

clusterDF$Species[!clusterDF$Species %in% alltraits$Species]
clusterDF <- merge(clusterDF,alltraits,by="Species",all.x=T)
clusterDF$cluster <- factor(clusterDF$cluster,levels=myOrder)
levels(clusterDF$cluster) <- mylabels

#scale between 0 and 1
stand <- function(x){
  x = x[!is.na(x)]
  (x-min(x))/(max(x)-min(x))
}

#get mean trait values for each cluster
clusterTraits <- ddply(clusterDF,.(cluster),summarise,
                       propDragon = median(ifelse(Suborder=="Anisoptera",1,0),na.rm=T),
                       meanFlightStart = median(Flight_start,na.rm=T),
                       propHabitat = mean(ifelse(Habitat=="Standing",1,0),na.rm=T),
                       propHabitatBreadth = mean(ifelse(HabitatBreadth=="generalist",1,0),na.rm=T),
                       meanRange = median(nuEuroGrids,na.rm=T),
                       meanDevTime = median(devTime,na.rm=T),
                       meanVoltinism = median(VoltinismProp,na.rm=T),
                       meanTemp = median(meanTemp,na.rm=T),
                       meanTempRange = median(sdTemp,na.rm=T),
                       meanHW = median(medHw,na.rm=T))

#plotting 
clusterTraitsM <- melt(clusterTraits,id=c("cluster"))                       
names(clusterTraitsM)[2:3] <- c("traits","meanValue")
clusterTraitsM$cluster <- factor(clusterTraitsM$cluster)

ggplot(clusterTraitsM)+
  geom_point(aes(x=traits,y=meanValue,colour=cluster),shape="square")+
  coord_flip()+
  theme(legend.position="top")+
  facet_wrap(~cluster)

#check for which traits the clusters differ
clusterDF$Habitat_Bin <- ifelse(clusterDF$Habitat=="Standing",1,0) 
anova(lm(Flight_start ~ factor(cluster),data=clusterDF))
anova(lm(nuEuroGrids ~ factor(cluster),data=clusterDF))
anova(lm(devTime ~ factor(cluster),data=clusterDF))
anova(lm(VoltinismProp ~ factor(cluster),data=clusterDF))
anova(lm(meanTemp ~ factor(cluster),data=clusterDF))###
anova(lm(sdTemp ~ factor(cluster),data=clusterDF))
anova(lm(medHw ~ factor(cluster),data=clusterDF))
anova(glm(Habitat_Bin ~ factor(cluster),family=binomial,data=clusterDF),test="Chisq")###
#table(clusterDF$Habitat,clusterDF$cluster)     
 
trendEstimates$cluster <- clusterDF$cluster[match(trendEstimates$Species,clusterDF$Species)]
summary(lm(trend~factor(cluster),data=trendEstimates,
           weights=1/sd))

#plot trait groups for each cluster

#plot wing length for each cluster
g1 <- ggplot(clusterDF)+
  geom_boxplot(aes(x=cluster,y=medHw,fill=cluster))+
  theme_bw()+ylab("Wing length (mm)")+xlab("Cluster")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

#plot temperature mean
g2 <- ggplot(clusterDF)+
  geom_boxplot(aes(x=cluster,y=meanTemp,fill=cluster))+
  theme_bw()+ylab("Temp pref")+xlab("")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

#plot habitat use
clusterDF$Habitat_Bin <- ifelse(clusterDF$Habitat=="Standing",0,1)
clusterDF$HabitatBreadth_Bin <- ifelse(clusterDF$HabitatBreadth=="generalist",1,0)
habitatSummary <- ddply(clusterDF,.(cluster),summarise,
                        RunningWater=sum(Habitat_Bin),
                        HabitatGeneralist=sum(HabitatBreadth_Bin),
                        nuSpecies = length(Species))

g3 <- ggplot(habitatSummary)+
  geom_bar(aes(x=cluster,y=RunningWater/nuSpecies,fill=cluster),stat="identity")+
  theme_bw()+ylab("Running water use")+xlab("Cluster")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

g4 <- ggplot(habitatSummary)+
  geom_bar(aes(x=cluster,y=HabitatGeneralist/nuSpecies,fill=cluster),stat="identity")+
  theme_bw()+ylab("Habitat generalism")+xlab("")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()
  
plot_grid(g1,g2,g3,g4,ncol=2)

###change point analysis######################################################################

library(changepoint)

for(i in 1)
m1=c(annualDF$mean[annualDF$Species=="Aeshna affinis"])
m1.amoc=cpt.mean(m1)
cpts(m1.amoc)
plot(m1.amoc)

ddply(annualDF,.(Species),function(x){
  m1.amoc=cpt.mean(x$mean)
  print(cpts(m1.amoc))
})
#none...
