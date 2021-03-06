#run first 2 boxes of trait_trend analysis
####get time series###########################################################################################

library(boot)

#convert time series into a list
annualDFS <- subset(annualDF,Year>1979)

#also remove Anax ephippiger, Boyeria irene and Oxygastra curtisii??
myTS <- dlply(annualDFS,.(Species),
              function(x){x[,"mean"]})

###preprocessing###########################################

myTSgr <- lapply(myTS,function(x){
  x/x[1]
})

myTSlogit <- lapply(myTS,function(x){
  logit(x)
})

#myTSlogit <- lapply(myTS,function(x){
#  logit(x)-logit(x)[1]
#})


myTSlgr <- lapply(myTS,function(x){
  log(x/x[1])
})

myTSl <- lapply(myTS,function(x){
  log(x)
})

myTSly <- lapply(myTS,function(x){
  x/x[37]
})

myTSm <- lapply(myTS,function(x){
  x/median(x)
})

#relative rate change
#e.g., diff(log(interest.rates), 1)
myTSrr <- lapply(myTS,function(x){
  diff(log(x),1)
})

#z-score
myTSz <- lapply(myTS,function(x){
  as.numeric(scale(x))
})

#combine and plot again
temp <- ldply(myTS,function(x){
  data.frame(Index=x,Year=1:length(x))
})
temp$Year <- temp$Year + min(annualDFS$Year)-1
#temp$Species <-temp$.id
temp2 <- temp
qplot(Year,Index,data=temp,colour=Species)+
  theme(legend.position = "none")

annualDFS <- subset(annualDF,Year>1979)

###TSclust################################################

library(TSclust)

#example:
#D1 <- diss(stocks, "COR")
#summary(D1)
#geometric P1 and P2 series are the closest ones
#structure P1 and P3 series become the closest ones
#we want want structural??

#direct values
#try dCORT and dCOR:1
IP.dis <- diss(myTSlogit, "CORT")
IP.dis <- diss(myTSm, "COR")
IP.dis <- diss(myTS, "COR")
IP.dis <- diss(myTSlogit, "COR")
IP.dis <- diss(myTSgr, "COR")
IP.dis <- diss(myTSlgr, "COR")
IP.dis <- diss(myTSz, "COR")
fit<- hclust(IP.dis)
plot(fit)
IP.clus <- pam(IP.dis, k = 5)$clustering

#relative rate change
IP.dis <- diss(myTSrr, "ACF")#performs badly
IP.dis <- diss(myTSrr, "PER")#performs badly
IP.dis <- diss(myTSrr, "AR.PIC")#performs badly
fit <- hclust(IP.dis)
plot(fit)
IP.clus <- pam(IP.dis, k = 6)$clustering

IP.clus <- pam(IP.dis, k = 4)$clustering
table(IP.clus)

####cluster evaluation##############################################

IP.clus <- list()
for(i in 2:20){
  IP.clus[[(i-1)]] <- pam(IP.dis, k = i)$clustering
}

IP.clusList <- ldply(IP.clus, function(x){
  getClusterStats(mydiss=IP.dis,
                  myclustering=x)})
IP.clusList$Cluster <- 2:20
plotTSclust(IP.clusList)

#find out when differences are closest to zero
findStationary <- function(x){
  min(which(diff(x)>0))
}

lapply(IP.clusList[,2:6],findStationary)

###best models#######################################################

getClusterStats(mydiss=IP.dis,myclustering=IP.clus)
cluster.stats(IP.dis,IP.clus)

#myTSgr 
#minCluster avBetween  avWithin  silWidth      dunn  sepIndex
#8  1.468354 0.8608622 0.1717292 0.2336572 0.4113406

#myTSlogit - this is the best!!!
#minCluster avBetween  avWithin  silWidth      dunn  sepIndex
#9  1.494418 0.8232091 0.2224314 0.2829076 0.45848822

#myTSlgr
#minCluster avBetween  avWithin silWidth      dunn  sepIndex
#6  1.469914 0.8133375 0.239221 0.2396137 0.3960618

#myTSm
#minCluster avBetween  avWithin  silWidth      dunn  sepIndex
#8  1.468354 0.8608622 0.1717292 0.2336572 0.4113406

#myTSz
#minCluster avBetween  avWithin  silWidth      dunn  sepIndex
#8  1.468354 0.8608622 0.1717292 0.2336572 0.4113406

#getClusterStats(mydiss=hc_sbd@distmat,myclustering=hc_sbd@cluster)

###cluster data frame#################################################

clusterDF <- data.frame(Species=names(IP.clus),
                        cluster=as.numeric(IP.clus))

annualDFS$cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]

table(clusterDF$cluster)
#1  2  3  4  5 
#35  9 12 10 11

plotCluster(annualDFS)

### mean correlation ####################################

#get pair-wise correlation of all species in each cluster

annualDFS_cluster1 <- subset(annualDFS,cluster==1)
temp <- acast(annualDFS_cluster1,Year~Species,value.var="mean")
temp <- cor(temp)
temp[lower.tri(temp,diag=FALSE)] <- NA
allCors <- reshape2::melt(temp,na.rm=T)
summary(allCors$value)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.05851  0.54236  0.72422  0.69127  0.85944  1.00000

annualDFS_cluster1 <- subset(annualDFS,cluster==2)
temp <- acast(annualDFS_cluster1,Year~Species,value.var="mean")
temp <- cor(temp)
temp[lower.tri(temp,diag=FALSE)] <- NA
allCors <- reshape2::melt(temp,na.rm=T)
summary(allCors$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.1449  0.3027  0.5157  0.5133  0.5909  1.0000 

annualDFS_cluster1 <- subset(annualDFS,cluster==3)
temp <- acast(annualDFS_cluster1,Year~Species,value.var="mean")
temp <- cor(temp)
temp[lower.tri(temp,diag=FALSE)] <- NA
allCors <- reshape2::melt(temp,na.rm=T)
summary(allCors$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1940  0.5553  0.6653  0.6780  0.8000  1.0000 

annualDFS_cluster1 <- subset(annualDFS,cluster==4)
temp <- acast(annualDFS_cluster1,Year~Species,value.var="mean")
temp <- cor(temp)
temp[lower.tri(temp,diag=FALSE)] <- NA
allCors <- reshape2::melt(temp,na.rm=T)
summary(allCors$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.1193  0.2535  0.3935  0.4757  0.7217  1.0000 

annualDFS_cluster1 <- subset(annualDFS,cluster==5)
temp <- acast(annualDFS_cluster1,Year~Species,value.var="mean")
temp <- cor(temp)
temp[lower.tri(temp,diag=FALSE)] <- NA
allCors <- reshape2::melt(temp,na.rm=T)
summary(allCors$value)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4180  0.7210  0.8570  0.8157  0.9198  1.0000 

#across all
temp <- acast(annualDFS,Year~Species,value.var="mean")
temp <- cor(temp)
temp[lower.tri(temp,diag=FALSE)] <- NA
allCors <- reshape2::melt(temp,na.rm=T)
summary(allCors$value)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.95250 -0.43898  0.07889  0.06906  0.57334  1.00000

###plot each cluster#######################################################

#order cluster by numbers of species
#clusts <- table(clusterDF$cluster)
#clusterOrder <- rev(names(clusts)[order(clusts)])

#examine cluster by cluster
ggplot(subset(annualDFS,cluster==1))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")
#general increase

ggplot(subset(annualDFS,cluster==2))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")
#mostly low high low

ggplot(subset(annualDFS,cluster==3))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")
#decrease

ggplot(subset(annualDFS,cluster==4))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")
#decrease then increase

ggplot(subset(annualDFS,cluster==5))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")
#decline at start 

ggplot(subset(annualDFS,cluster==6))+
  geom_line(aes(x=Year,y=mean))+
  facet_wrap(~Species,scales="free")


#rescale each to units of standard deviation
annualDFS_scaled <- ddply(annualDFS,.(Species),function(x){
  x$mean_scaled <- as.numeric(scale(x$mean))
  return(x)
})

#change order
annualDFS_scaled$cluster <- factor(annualDFS$cluster,levels=c(1,4,2,3,5))
levels(annualDFS_scaled$cluster) <- c(1,2,3,4,5)

ggplot(annualDFS_scaled)+
  geom_line(aes(x=Year,y=mean_scaled,colour=Species),
            alpha=0.7)+
  facet_wrap(~cluster)+
  scale_color_discrete()+
  theme_classic()+
  ylab("Scaled occupancy probability")+
  theme(legend.position = "none")

ggsave("plots/Cluster_species_ts.png",width=7,height=5)

###cluster panel####################################

#rescale indicies to 1
annualDFS <- ddply(annualDFS,.(Species,Code),function(x){
  rescaleFactor = 1/x$mean[x$Year==1980]
  x$sMean = x$mean * rescaleFactor
  x$sSD = x$sd * rescaleFactor
  x$rescale <- rescaleFactor
  x$mean1 <- x$mean[x$Year==1980]
  x$meanMean <- mean(x$mean)
  x$sdMean <- sd(x$mean)
  
  return(x)
})

### outliers ##########################################

annualDFS$Cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]

ggplot(annualDFS)+
  geom_line(aes(x=Year,y=mean/mean1,colour=Species))+
  facet_wrap(~cluster,scales="free")+
  theme(legend.position = "none")
  
#identify outliers
#need formal way to identify outliers
#(1)get standard deviation of occupancies of species in each group 
speciesSD <- ddply(annualDFS,.(Cluster,Species),summarise,
                   speciesSD = sd(mean),
                   speciesCOV=sd(mean)/mean(mean))
hist(speciesSD$speciesCOV[speciesSD$Cluster==2])
speciesSD <- arrange(speciesSD,speciesCOV)

library(StatMeasures)
outliers(speciesSD$speciesCOV[speciesSD$Cluster==1])
speciesSD$Species[speciesSD$Cluster==1][c(1,13)]

#get list of outlier species
probSpecies <- as.character(unlist(dlply(speciesSD,.(Cluster),function(x){
  ids <- outliers(x$speciesCOV)$idxOutliers
  return(x$Species[ids])
})))

probSpecies <- c(probSpecies,"Anax ephippiger","Boyeria irene","Oxygastra curtisii")
#"Sympetrum meridionale" "Lestes barbarus"       "Sympetrum flaveolum"   "Coenagrion scitulum"   "Anax ephippiger"      
#[6] "Boyeria irene"         "Oxygastra curtisii" 

annualDFS <- subset(annualDFS,!Species %in% probSpecies)

myorder <- c("1","4","2","3","5")
myrandomCI$Cluster <- factor(myrandomCI$Cluster,levels=myorder)
myrandomCI$ClusterF <- myrandomCI$Cluster
mylabels <- c("(1) increasing","(2) increasing late","(3) mixed","(4) decreasing late","(5) decreasing early")

levels(myrandomCI$ClusterF) <- mylabels

mycols <- wes_palette("Zissou1", 
                      n=10,
                      type="continuous")[c(1,4,5,8,10)]

#run code below to get graph

ggsave("plots/RandomCI_5_wo_outliers.png",width=7,height=4)
ggsave("plots/RandomCI_5_wo_outliers_reclusters.png",width=7,height=4)
####bootstrap#############################################################

# #add cluster to original data frame
# annualDFS$Cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]
# 
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# 
# #bootstrap original values within each cluster at each step
# library(boot)
# 
# myMean <- function(data, indices) {
#   
#   #randomly pick number
#   data$rMean <- apply(data,1,function(x)rnorm(1,
#                                               as.numeric(x["mean"]),
#                                               as.numeric(x["sd"])))
#   
#   #rescale number 
#   data$srMean <- data$rMean/data$mean1
#   
#   #pick bootstrp selection
#   d <- data[indices,] 
#   temp <-  gm_mean(d$srMean)
#   
#   return(temp)
# }
# 
# # bootstrapping with 1000 replications
# applyBoot <- function(mydata){
#   results <- boot(data=mydata, statistic=myMean,
#                   R=5000)
#   out <- boot.ci(results, type="bca")
#   data.frame(lowerQ = out$bca[1,4],
#              upperQ = out$bca[1,5])
# }
# 
# #applyBoot(annualDFS)
# 
# mybootCI <- ddply(annualDFS,.(Year,Cluster),function(x)applyBoot(x))
# 
# 
# #plot
# ggplot(mybootCI)+
#   geom_ribbon(aes(x=Year,ymin=lowerQ,ymax=upperQ))+
#   facet_wrap(~Cluster,scales="free")

###random scale####################################################################

#add cluster to original data frame
annualDFS$Cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]

gm_mean = function(x, na.rm=TRUE){
     exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
myMean <- function(mydata) {
  
  mydata$rMean <- apply(mydata,1,function(x)rnorm(1,
                                           as.numeric(x["mean"]),
                                           as.numeric(x["sd"])))
  
  mydata$rMean <- mydata$rMean/mydata$mean1
  
  return(gm_mean(mydata$rMean))
  
}


# draw random number 1000 times and calculate the geometric mean
applyRandomCI <- function(mydata,n.sim=1000){
  
  myMeans <- NA
  
  for(i in 1:n.sim){
    myMeans[i] <- myMean(mydata)
  }
  
  data.frame(lowerQ=quantile(myMeans,c(0.025)),
             medianQ=quantile(myMeans,c(0.5)),
             upperQ=quantile(myMeans,c(0.975)))
  
}

#
myrandomCI <- ddply(annualDFS,.(Year,Cluster),function(x)applyRandomCI(x))

require(wesanderson)

#5 clusters
myorder <- c("1","4","2","3","5")
myrandomCI$Cluster <- factor(myrandomCI$Cluster,levels=myorder)
myrandomCI$ClusterF <- myrandomCI$Cluster
mylabels <- c("(1) increasing","(2) increasing late","(3) mixed","(4) decreasing late","(5) decreasing early")

levels(myrandomCI$ClusterF) <- mylabels

mycols <- wes_palette("Zissou1", 
                      n=10,
                      type="continuous")[c(1,4,5,8,10)]

#4 clusters
myorder <- c("1","4","2","3")
myrandomCI$Cluster <- factor(myrandomCI$Cluster,levels=myorder)
myrandomCI$ClusterF <- myrandomCI$Cluster
mylabels <- c("(1) increasing","(2) increasing late","(4) mixed","(5) decreasing")
levels(myrandomCI$ClusterF) <- mylabels

mycols <- wes_palette("Zissou1", 
                      n=10,
                      type="continuous")[c(1,4,8,10)]
ggsave("plots/RandomCI_4clusters.png",width=7,height=4)

#6 clusters
myorder <- c("1","5","4","2","3","6")
myrandomCI$Cluster <- factor(myrandomCI$Cluster,levels=myorder)
myrandomCI$ClusterF <- myrandomCI$Cluster
mylabels <- c("(1) increasing","(2) increasing late","(3) increasing early","(4) mixed","(5) decreasing late","(6) decreasing early")

levels(myrandomCI$ClusterF) <- mylabels

mycols <- wes_palette("Zissou1", 
                      n=10,
                      type="continuous")[c(1,4,5,7,8,10)]

ggsave("plots/RandomCI_6clusters.png",width=7,height=4)

#plot
gA <- ggplot(myrandomCI)+
  geom_ribbon(aes(x=Year,ymin=lowerQ,ymax=upperQ,
                  fill=factor(Cluster)))+
  facet_wrap(~ClusterF,scales="free",ncol=3)+
  scale_fill_manual(values=rev(mycols))+
  theme_classic()+
  theme(legend.position = "none")
  
gA
# #plus GAM derivatives
# myrandomCI$Signif <- clusterDeriv$Signif[match(interaction(myrandomCI$Cluster,myrandomCI$Year),interaction(clusterDeriv$Cluster,clusterDeriv$data))]
# 
# gA <- ggplot(myrandomCI)+
#   geom_ribbon(aes(x=Year,ymin=lowerQ,ymax=upperQ,
#                   fill=factor(Cluster)))+
#   facet_wrap(~ClusterF,scales="free",ncol=3)+
#   scale_fill_manual(values=rev(mycols))+
#   theme_classic()+
#   ylab("Mean occupancy index")+
#   theme(legend.position = "none")+
#   geom_point(data=subset(myrandomCI,Signif!="insig"),aes(x=Year,y=upperQ),color="black",size=0.5,shape=18)+
#   geom_point(data=subset(myrandomCI,Signif!="insig"),aes(x=Year,y=lowerQ),color="black",size=0.5,shape=18)


#add on geometric mean of the assemblage-level results



# ###GAM on fits############################################
# 
# #https://stats.stackexchange.com/questions/84325/calculating-bootstrap-confidence-intervals-on-second-derivatives-of-a-gam-object
# library(gratia)
# library(mgcv)
# 
# mod <- gam(medianQ ~ s(Year), data = subset(myrandomCI,Cluster==1), method = "REML")
# 
# ## first derivative
# firstD <- derivatives(mod, order=1,type = "central",interval="confidence")
# ggplot(firstD)+
#   geom_ribbon(aes(x=data,ymin=lower,ymax=upper))+
#   geom_hline(yintercept=0,slope=0,linetype="dashed")
# 
# clusterDeriv <- ddply(myrandomCI,.(Cluster),function(x){
#   mod <- gam(medianQ ~ s(Year), data = x, method = "REML")
#   mydf <- expand.grid(Year=sort(unique(x$Year)))
#   derivatives(mod, order=1,type = "central",interval="confidence",newdata=mydf)
# })
# 
# #second derivative
# secondD <- derivatives(mod, order=2,type = "central",interval="confidence",
#                        newdata=data.frame(Year=unique(annualDFS$Year)))
# ggplot(secondD)+
#   geom_ribbon(aes(x=data,ymin=lower,ymax=upper))+
#   geom_hline(yintercept=0,linetype="dashed")
# 
# clusterDeriv <- ddply(myrandomCI,.(Cluster),function(x){
#   mod <- gam(medianQ ~ s(Year), data = x, method = "REML")
#   mydf <- expand.grid(Year=sort(unique(x$Year)))
#   derivatives(mod, order=2,type = "central",interval="confidence",newdata=mydf)
# })
# 
# clusterDeriv$Signif <- "insig"
# clusterDeriv$Signif[clusterDeriv$lower>0 & clusterDeriv$upper>0] <- "Positive"
# clusterDeriv$Signif[clusterDeriv$lower<0 & clusterDeriv$upper<0] <- "Negative"
# 
# ###GAM on data############################################
# 
# library(gratia)
# library(mgcv)
# library(gamm4)
# 
# mod <- gamm4(mean ~ s(Year),random=~(1|Species),data = subset(annualDFS,Cluster==3), REML=T)
# 
# ## first derivative
# firstD <- derivatives(mod$gam, order=1,type = "central",interval="confidence",newdata=data.frame(Year=unique(annualDFS$Year)))
# 
# ggplot(firstD)+
#   geom_ribbon(aes(x=data,ymin=lower,ymax=upper))+
#   geom_hline(yintercept=0,linetype="dashed")
# 
# #fit to each cluster
# clusterDeriv <- ddply(annualDFS,.(Cluster),function(x){
#   #mod <- gam(mean ~ s(Year)+Species, data = x, method = "REML")
#   mod <- gamm4(mean ~ s(Year),random=~(1|Species), data = x,REML = TRUE)
#   mydf <- expand.grid(Year=sort(unique(x$Year)))
#   derivatives(mod$gam, order=1,type = "central",interval="confidence",newdata=mydf)
# })
# 
# clusterDeriv$Signif <- "Signif"
# clusterDeriv$Signif[clusterDeriv$lower<0 & clusterDeriv$upper>0] <- "Insignif"
# 
# #second derivative
# secondD <- derivatives(mod$gam, order=2,type = "central",interval="confidence",
#                        newdata=data.frame(Year=unique(annualDFS$Year)))
# ggplot(secondD)+
#   geom_ribbon(aes(x=data,ymin=lower,ymax=upper))+
#   geom_hline(yintercept=0,linetype="dashed")
# 
# clusterDeriv <- ddply(annualDFS,.(Cluster),function(x){
#   #mod <- gam(mean ~ s(Year)+Species, data = x, method = "REML")
#   mod <- gamm4(mean ~ s(Year),random=~(1|Species), data = x,REML = TRUE)
#   mydf <- expand.grid(Year=sort(unique(x$Year)))
#   derivatives(mod$gam, order=2,type = "central",interval="confidence",newdata=mydf)
# })
# 
# clusterDeriv$Signif <- "Signif"
# clusterDeriv$Signif[clusterDeriv$lower<0 & clusterDeriv$upper>0] <- "Insignif"
# 
# 
# ###traits cluster test###################################
# 
#merge clusters with traits
clusterDF$Species[!clusterDF$Species %in% trendEstimates$Species]
clusterDF <- merge(clusterDF,trendEstimates,by="Species",all.x=T)
clusterDF$cluster <- factor(clusterDF$cluster)

#plot habitat use
habitatSummary <- ddply(clusterDF,.(cluster),summarise,
                        RunningWater=sum(river),
                        Bog=sum(bog),
                        HabitatGeneralist=sum(Generalism),
                        nuSpecies = length(Species))


#order clusters
clusterDF$cluster <- factor(clusterDF$cluster,levels=myorder)
habitatSummary$cluster <- factor(habitatSummary$cluster,levels=myorder)
levels(clusterDF$cluster) <- 1:5
levels(habitatSummary$cluster) <- 1:5



#how does cluster relate to population trends?
trendEstimates$cluster <- clusterDF$cluster[match(trendEstimates$Species,clusterDF$Species)]
summary(lm(trend~factor(cluster),data=trendEstimates))
#Multiple R-squared:  0.5863,	Adjusted R-squared:  0.5634 

###trait cluster plots##########################################################

#reverse the clusters
clusterDF$cluster <- factor(clusterDF$cluster,levels=5:1)
habitatSummary$cluster <- factor(habitatSummary$cluster,levels=5:1)

#plot wing length for each cluster
g1 <- ggplot(clusterDF)+
  geom_boxplot(aes(x=cluster,y=medHw,fill=cluster))+
  theme_bw()+ylab("Wing length")+xlab("Cluster")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

#temperature mean
g2 <- ggplot(clusterDF)+
  geom_boxplot(aes(x=cluster,y=meanTemp,fill=cluster))+
  theme_bw()+ylab("Temp pref")+xlab("")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

#voltinism
g3 <- ggplot(clusterDF)+
  geom_boxplot(aes(x=cluster,y=VoltinismProp,fill=cluster))+
  theme_bw()+ylab("Voltinism")+xlab("")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

#flight start
g4 <- ggplot(clusterDF)+
  geom_boxplot(aes(x=cluster,y=Flight_start,fill=cluster))+
  theme_bw()+ylab("Flight start")+xlab("")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

#habitat variables
g5 <- ggplot(habitatSummary)+
  geom_bar(aes(x=cluster,y=RunningWater/nuSpecies,fill=cluster),stat="identity",colour="black")+
  theme_bw()+ylab("Proportion of riverine species")+xlab("Cluster")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

g6 <- ggplot(habitatSummary)+
  geom_bar(aes(x=cluster,y=Bog/nuSpecies,fill=cluster)
           ,stat="identity",color="black")+
  theme_bw()+ylab("Proportion of bog species")+xlab("")+
  scale_fill_manual(values=mycols)+
  theme(legend.position = "none")+
  coord_flip()

gB <- plot_grid(g1,g2,g5,g6,ncol=2)
gB

plot_grid(gA,gB,
          axis="tr",
          align="v",
          labels=c("A","B"),
          ncol=1)

ggsave("plots/RandomCI_traits_5.png",width=7,height=8)

###single cluster tests#############################################

chisq.test(clusterDF$cluster,clusterDF$Habitat)
chisq.test(clusterDF$cluster,clusterDF$Overwintering)
chisq.test(clusterDF$cluster,clusterDF$Generalism)
chisq.test(clusterDF$cluster,clusterDF$river)
chisq.test(clusterDF$cluster,clusterDF$bog)

#cluster association one by one
clusterDF <- cbind(clusterDF,model.matrix(~cluster-1,data=clusterDF))

glm1 <- glm(cluster1 ~ medHw + meanTemp + river + bog, family=binomial,data=clusterDF)
summary(glm1)

glm1 <- glm(cluster2 ~  medHw + meanTemp + river + bog, family=binomial,data=clusterDF)
summary(glm1)

glm1 <- glm(cluster3 ~  medHw + meanTemp + river + bog, family=binomial,data=clusterDF)
summary(glm1)

glm1 <- glm(cluster4 ~  medHw + meanTemp + river + bog, family=binomial,data=clusterDF)
summary(glm1)

glm1 <- glm(cluster5 ~  medHw + meanTemp + river + bog, family=binomial,data=clusterDF)
summary(glm1)

glm1 <- glm(cluster6 ~  medHw + meanTemp + river + bog, family=binomial,data=clusterDF)
summary(glm1)

###lda#############################################################

#http://dwoll.de/rexrepos/posts/regressionMultinom.html#using-vglm-from-package-vgam
#http://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/

mysample <- sample(1:nrow(clusterDF),0.25*nrow(clusterDF))
trainData <- clusterDF[-mysample,]
testData <- clusterDF[mysample,]

library(MASS)
lda1 <- lda(cluster ~ nuEuroGrids + meanTemp + Flight_start +
              VoltinismProp + river + bog  Flow + Generalism,
            data=trainData)
lda1
plot(lda1) 

# Make predictions
predictions <- lda1 %>% predict(testData)

# Model accuracy
mean(predictions$class==testData$cluster)

#Model accuracy split by class
testData$acc <- predictions$class==testData$cluster
tapply(testData$acc,testData$cluster,mean)
#driven by group 1
table(trainData$cluster)

#Better than chance?
comp <- cbind(predictions$class,testData$cluster)
library(irr)
agree(comp, tolerance=0)
kappa2(comp)

#run over multiple possible random samples and get distribution of accuracy

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

###dwtclust####################################################

#https://damien-datasci-blog.netlify.com/post/time-series-clustering-with-dynamic-time-warp/

library("dtwclust")

#SBD
hc_sbd <- tsclust(myTS, type = "partitional", k=2:20L,
                  preproc = zscore, 
                  distance = "sbd",centroid = "shape")
#DTW
hc_sbd <- tsclust(myTS, ype="partitional", k= 2:20L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "dba")

hc_sbd <- tsclust(myTS, ype="partitional", k= 2:20L,
                  preproc = zscore, 
                  distance = "dtw", centroid = "dba")

#pam
hc_sbd <- tsclust(myTS, ype="partitional", k= 2:20L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "pam")

#GAK
hc_sbd <- tsclust(myTS, ype="partitional", k= 2:20L,
                  preproc = zscore, 
                  distance = "gak", centroid = "pam")

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
out <- ldply(unique(temp$Param),function(x){getBestCluster(temp,x)})
median(out$midpts)
#https://en.wikipedia.org/wiki/Stationary_point

###DTW pick K#####################################################################

library("dtwclust")

#SBD
hc_sbd <- tsclust(myTS, type = "partitional", k=6L,
                  preproc = zscore, 
                  distance = "sbd",centroid = "shape")

#DTW
hc_sbd <- tsclust(myTS, type="partitional", k= 6L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "dba")

hc_sbd <- tsclust(myTS, ype="partitional", k= 6L,
                  preproc = zscore, 
                  distance = "dtw", centroid = "dba")

#pam
hc_sbd <- tsclust(myTS, type="partitional", k= 5L,
                  preproc = zscore, 
                  distance = "dtw_basic", centroid = "pam")

#gak
hc_sbd <- tsclust(myTS, type="partitional", k= 6L,
                  preproc = zscore, 
                  distance = "gak", centroid = "pam")

plot(hc_sbd)
plot(hc_sbd, type = "sc")
#plot(hc_sbd, type = "series", clus = 1L)
plot(hc_sbd, type = "series")
plot(hc_sbd, type = "centroids")

#check
clusterDF <- data.frame(Species=names(myTS),cluster=hc_sbd@cluster)

annualDFS$cluster <- clusterDF$cluster[match(annualDFS$Species,clusterDF$Species)]

plotCluster(annualDFS)

###other plots############################################

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

library(boot)

myMean <- function(data, indices) {
  
  #pick bootstrap selection
  d <- data[indices,]
  
  #randomly pick occurence number
  d$rMean <- apply(d,1,function(x)rnorm(1,
                                        as.numeric(x["mean"]),
                                        as.numeric(x["sd"])))
  
  #fit a gam to each species
  temp <- ddply(d,.(Species),function(x){
    gam1 <- loess(rMean ~ Year, data=x)
    x$fits <- gam1$fitted
    return(x)
  })
  
  #get means of the fits per year
  as.numeric(tapply(temp$fits,temp$Year,gm_mean))
  
}

# bootstrapping with 1000 replications
applyBoot <- function(mydata){
  
  results <- boot(data=mydata, statistic=myMean,
                  R=5000,strata=factor(mydata$Species))
  
  ldply(1:37,function(i){
    out <- boot.ci(results, type="perc",index=i)
    data.frame(Year=i,
               lowerQ = out$perc[1,4],
               upperQ = out$perc[1,5])})
}

applyBoot(annualDFS)

mybootCI <- dlply(annualDFS,.(Cluster),function(x)applyBoot(x))


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

