#linking traits to trends

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/sparta_wrapper_functions.R')
library(plyr)
library(ggplot2)
library(reshape2)
library(cowplot)
###merge data###################################################################

setwd("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects")

#get traits data
load("alltraits.RData")
#limit to those with complete cases?

#get trends data
#load("model-outputs/trendEstimatesNational.RData")
#load("model-outputs/trendEstimates.RData")
trendEstimates <- readRDS("model-outputs/modelTrends_nation_state_trends.rds")
trendEstimates$Species <- gsub("out_dynamic_nation_state_adult_","",trendEstimates$file)
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
1.0034478^35 

summary(trendEstimates$trend[trendEstimates$Trend=="significant decrease"])
1-(1-0.003953)^35

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

#cluster them into groups of species with similar dyanamics

# modelSummary2 <- nationTrends
# library(reshape2)

# #cluster on occurences#
# 
# mydata <- ddply(modelSummary2,.(Species),function(x){
#   require(zoo)
#   roll.mean <- rollmean(x$mean,5,align="center")
#   len <- length(roll.mean)
#   data.frame(Species=rep(unique(x$Species),len),roll.mean,Year=1:len)
# })
# 
# mydata <- acast(mydata,Species~Year,value.var="roll.mean")
# mydata[is.na(mydata)]<-0
# #apply cluster analysis
# d <- dist(mydata, method = "euclidean") 
# fit <- hclust(d, method="ward.D2")
# plot(fit) # display dendogram
# 
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares") 
# 
# #pam
# library(cluster)
# fit <- kmeans(mydata, 6) 
# 
# # get cluster means
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster) 
# mydata$Species <- row.names(mydata)
# 
# #add to dataframe 
# modelSummary2$cluster <- mydata$fit.cluster[match(modelSummary2$Species,mydata$Species)]
# 
# 
# #order
# modelSummary2$cluster <- as.factor(modelSummary2$cluster)
# 
# modelSummary2$cluster <- factor(modelSummary2$cluster, levels=c("2","4","5","1","6","3"))
# levels(modelSummary2$cluster) <- c("1","2","3","4","5","6")
# ddply(modelSummary2,.(cluster),summarise,nuS = length(unique(Species)))
# levels(modelSummary2$cluster) <- c("18 sp","15 sp","12 sp","6 sp","11 sp","13 sp")
# 
# unique(subset(modelSummary2,cluster==6)$Species)
# #"Aeshna juncea"          "Coenagrion hastulatum"  "Gomphus pulchellus"     "Ischnura pumilio"      
# #"Lestes dryas"           "Leucorrhinia rubicunda"
# 
# #plot
# ggplot(modelSummary2,aes(x=Year,y=mean))+
#   geom_line(aes(colour=Species))+
#   facet_wrap(~cluster,scales="free")+
#   #geom_smooth(aes(colour=Species),se=F)+
#   theme_bw()+ ylab("Occupancy")+
#   theme(legend.position="none")
# 
# 
# 
# #cluster on growth rates##
# 
# library(zoo)
# gee<- c(1:10)
# rollmean(gee,5,align="right")
# 
# mydata <- ddply(modelSummary2,.(Species),function(x){
#   len = length(x$Year)
#   #roll.mean <- rollmean(x$mean,5,align="right")
#   #len = length(roll.mean)
#   growth = x$mean[2:len]/x$mean[1]
#   data.frame(Species=rep(unique(x$Species),(len-1)),growth,Year=2:len)
# })
# mydata <- acast(mydata,Species~Year,value.var="growth")
# mydata[is.na(mydata)]<-0
# 
# #apply cluster analysis
# d <- dist(mydata, method = "euclidean") 
# fit <- hclust(d, method="ward.D2")
# plot(fit) # display dendogram
# 
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares") 
# 
# #choose cluster
# fit <- kmeans(mydata, 6) 
# table(fit$cluster)
# 
# # get cluster means
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster) 
# mydata$Species <- row.names(mydata)
# 
# #add to dataframe 
# modelSummary2$cluster <- mydata$fit.cluster[match(modelSummary2$Species,mydata$Species)]
# table(modelSummary2$cluster)
# 
# #plot
# ggplot(modelSummary2)+
#   geom_line(aes(x=Year,y=mean,colour=Species))+
#   facet_wrap(~cluster)+
#   theme_bw()+ ylab("Occupancy")+
#   theme(legend.position="none")
# 
# #plot change in occupancy
# mydataM <- melt(mydata,id=c("fit.cluster","Species"))
# 
# mydataM$Year <- as.numeric(gsub("X","",mydataM$variable))+1979
# mydataM$Year <- as.numeric(gsub("X","",mydataM$variable))
# 
# ggplot(mydataM)+
#   geom_line(aes(x=Year,y=value,colour=Species))+
#   facet_wrap(~fit.cluster,nrow=2)+
#   theme_bw()+ ylab("Occupancy")+
#   theme(legend.position="none")
# 
# #with order
# mydataM$fit.cluster <- as.factor(mydataM$fit.cluster)
# mydataM$fit.cluster <- factor(mydataM$fit.cluster, levels=c("5","1","6","2","4","3"))
# levels(mydataM$fit.cluster) <- c("1","2","3","4","5","6")
# ddply(mydataM,.(fit.cluster),summarise,nuS = length(unique(Species)))
# levels(mydataM$fit.cluster) <- c("18 sp","15 sp","12 sp","6 sp","11 sp","13 sp")
# 
# unique(subset(mydataM,cluster==6)$Species)
# 
# #two plots
# ggplot(subset(mydataM,fit.cluster %in% c("18 sp","15 sp","12 sp")))+
#   geom_line(aes(x=Year,y=value,colour=Species))+
#   facet_wrap(~fit.cluster,nrow=1)+
#   theme_bw()+ ylab("Occupancy")+
#   geom_smooth(aes(x=Year,y=value),se=F,colour="black")+
#   theme(legend.position="none")
# 
# ggplot(subset(mydataM, fit.cluster %in% c("6 sp","11 sp","13 sp")))+
#   geom_line(aes(x=Year,y=value,colour=Species))+
#   facet_wrap(~fit.cluster,nrow=1)+
#   geom_smooth(aes(x=Year,y=value),se=F,colour="black")+
#   theme_bw()+ ylab("Occupancy")+
#   theme(legend.position="none")

#using proper time series clustering

#convert time series into a list
annualDFS <- subset(annualDF,Year>1982)
myTS <- dlply(annualDFS,.(Species),
              function(x){x[,"mean"]})

#preprocessing
myTS <- lapply(myTS,function(x){
  len <- length(x)
  x[2:len]/mean(x)
})
annualDFS <- subset(annualDF,Year>1983)
#or smoothing first???

#different options:

###dwtclust#######################################################

library("dtwclust")

#hierarchical
hc_sbd <- tsclust(myTS, type = "h", k = 2:20L,
                  preproc = zscore, seed = 899,
                  distance = "sbd", centroid = shape_extraction,
                  control = hierarchical_control(method = "average"))

hc_sbd <- tsclust(myTS, type = "h", k = 7L,
                  preproc = zscore, seed = 899,
                  distance = "sbd", centroid = shape_extraction,
                  control = hierarchical_control(method = "average"))


hc_sbd <- tsclust(myTS, type = "h", k = 2:20L,
                  preproc = NULL, seed = 899,
                  distance = "sbd", centroid = shape_extraction,
                  control = hierarchical_control(method = "average"))


hc_sbd <- tsclust(myTS, type = "h", k = 5L,
                  preproc = NULL, seed = 899,
                  distance = "sbd", centroid = shape_extraction,
                  control = hierarchical_control(method = "average"))

#partitional
hc_sbd <- tsclust(myTS, type = "partitional", k=2:20L,
                  preproc = zscore, 
                  distance = "sbd",centroid = "shape")

hc_sbd <- tsclust(myTS, type = "partitional", k=5L,
                  preproc = zscore, 
                  distance = "sbd",centroid = "shape")

hc_sbd <- tsclust(myTS, type = "partitional", k=2:20L,
                  preproc = NULL, 
                  distance = "sbd",centroid = "shape")

hc_sbd <- tsclust(myTS, type = "partitional", k=4L,
                  preproc = NULL, 
                  distance = "sbd",centroid = "shape")

#other
hc_sbd <- tsclust(myTS, ype="partitional", k= 4L,
                  distance = "dtw_basic", centroid = "dba")

hc_sbd <- tsclust(myTS, ype="partitional", k= 4L,
                  distance = "dtw_basic", centroid = "pam")


#By default, the dendrogram is plotted in hierarchical clustering
plot(hc_sbd)
plot(hc_sbd, type = "sc")
#plot(hc_sbd, type = "series", clus = 1L)
plot(hc_sbd, type = "series")
plot(hc_sbd, type = "centroids")

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

#check
clusterDF <- data.frame(Species=names(myTS),
                        cluster=hc_sbd@cluster)

annualDFS$cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]
plotCluster(annualDFS)

###TSclust##############################################################

library("TSclust")

#scale each time series
mysTS <- lapply(myTS,function(x)as.numeric(scale(x)))

#dissimilarity possibilities
?diss

IP.dis <-diss(mysTS, "ACF", p = 0.05)
IP.dis <-diss(mysTS, "COR")
IP.dis <-diss(mysTS, "DTWARP")
IP.dis <-diss(mysTS, "FRECHET")
IP.dis <-diss(mysTS, "PDC")
IP.dis <- diss(mysTS, "INT.PER")


#hierarchical
IP <- hclust(IP.dis)
plot(IP)
IP <- cutree(IP, k = 6)

#pam
IP <- pam(IP.dis, k = 6)$clustering

#make data frame - from hclust
clusterDF <- data.frame(Species=names(myTS),
                        cluster=as.numeric(IP))

annualDF$cluster <- clusterDF$cluster[match(annualDF$Species,clusterDF$Species)]
plotCluster(annualDF)

###plot clusters#######################################################

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

ggplot(subset(annualDFS,cluster==7))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")

#plot each cluster series
myCentroids <- data.frame(Year=rep(sort(unique(annualDFS$Year)),length(hc_sbd@centroids)),
                          Cluster=rep(1:length(hc_sbd@centroids),each=length(unique(annualDFS$Year))),
                          ts=do.call(c,hc_sbd@centroids))
myCentroids$Cluster <- factor(myCentroids$Cluster,levels=1:length(clusterOrder))

ggplot(data=myCentroids,aes(x=Year,y=ts))+
  geom_smooth(aes(colour=factor(Cluster),fill=factor(Cluster)))+
  facet_wrap(~Cluster,ncol=1)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("z-score occupancy prop")

table(clusterDF$cluster)
#1  2  3  4 
#23 26 14 14



###traits clusters########################################################

clusterDF$Species[!clusterDF$Species %in% alltraits$Species]
clusterDF <- merge(clusterDF,alltraits,by="Species",all.x=T)

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
