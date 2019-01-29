#linking traits to trends

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/R/sparta_wrapper_functions.R')

######################################################################

#get traits data
load("alltraits.RData")
#limit to those with complete cases?

load("trendEstimates.RData")

#########################################################################

#plot of traits vs trends
trendEstimates <- merge(trendEstimates,alltraits,by="Species")

#temperature preference plot
ggplot(trendEstimates,aes(x=TMean,y=trend))+
  geom_point()+
  facet_wrap(State~Stage)

summary(lm(trend~TMean,data=subset(trendEstimates,Stage=="adult"),
                                   weights=1/trend_sd))#sig, temp mean higher tends
#positive effect

#habitat preference
ggplot(trendEstimates,aes(x=Habitat.y,y=trend))+
  geom_boxplot()+
  facet_wrap(State~Stage)

summary(lm(trend~Habitat.y,data=subset(trendEstimates,Stage=="adult"),
           weights=1/trend_sd))#sig - running, higher trends

#range size
ggplot(trendEstimates,aes(x=nuEuroGrids,y=trend))+
  geom_point()+
  facet_wrap(State~Stage)

summary(lm(trend~nuEuroGrids,data=subset(trendEstimates,Stage=="adult"),
           weights=1/trend_sd))#ns

#generation time
ggplot(trendEstimates,aes(x=development.time..mean.,y=trend))+
  geom_point()+
  facet_wrap(State~Stage)

summary(lm(trend~development.time..mean.,data=subset(trendEstimates,Stage=="adult"),
           weights=1/trend_sd))#marginal


####################################################################

#plot trait vs trends
q1<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=TMean,y=trend,colour=State))+
  geom_point()+
  theme_bw()+
  ylab("Population trend")+xlab("Temperature preference")+
  facet_grid(~State)+stat_smooth(method="lm")+
  theme(legend.position="none")
q2<-ggplot(subset(trendEstimates,Stage=="adult"),aes(x=Habitat.y,y=trend,fill=State))+
  geom_boxplot()+
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

library(cowplot)
plot_grid(q1,q2,q3,q4,q5,q6,align="v",ncol=1)

library(gridExtra)
grid.arrange(q5,q6,ncol=1)

#####################################################################

#load and merge with random communities
load("randomMatrix.RData")
randomMatrixM<-melt(randomMatrix,id=c("State","Year","Species"))
randomMatrixM<-arrange(randomMatrixM,Species,Year,State)
randomMatrixM$Species[!randomMatrixM$Species %in% alltraits$Species]
randomTrends <- merge(randomMatrixM,alltraits,by="Species")

randomTrends$State <- as.factor(randomTrends$State)
levels(randomTrends$State)<-c("Bavaria","North Rhine-Westphalia","Saarland","Saxony","Schleswig Holstein")


#get average traits value per state and year
#tmean
tmeansMeans <- ddply(randomTrends,.(State,Year,variable),summarise,
                     tmean = weighted.mean(TMean,value))
tmeansMeans <- ddply(tmeansMeans,.(State,Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))
q1<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Temperature preference")+
  facet_grid(~State)+
  theme(legend.position="none")

#german distribution
tmeansMeans <- ddply(randomTrends,.(State,Year,variable),summarise,
                     tmean = weighted.mean(total,value))
tmeansMeans <- ddply(tmeansMeans,.(State,Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

q2<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Geographic range (Germany)")+
  facet_grid(~State)+
  theme(legend.position="none")

#european distribution
tmeansMeans <- ddply(randomTrends,.(State,Year,variable),summarise,
                     tmean = weighted.mean(nuEuroGrids,value))
tmeansMeans <- ddply(tmeansMeans,.(State,Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

q3<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Geographic range (Europe)")+
  facet_grid(~State)+
  theme(legend.position="none")

#habitat preference
randomTrends$Habitat.y <- ifelse(randomTrends$Habitat.y=="Standing",0,1)
tmeansMeans <- ddply(randomTrends,.(State,Year,variable),summarise,
                     tmean = weighted.mean(Habitat.y,value))
tmeansMeans <- ddply(tmeansMeans,.(State,Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

q4<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Habitat preference (lotic v lentic)")+
  facet_grid(~State)+
  theme(legend.position="none")

#development time
tmeansMeans <- ddply(randomTrends,.(State,Year,variable),summarise,
                     tmean = weighted.mean(development.time..mean.,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(State,Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

q5<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Development time")+
  facet_grid(~State)+
  theme(legend.position="none")


#wing length
tmeansMeans <- ddply(randomTrends,.(State,Year,variable),summarise,
                     tmean = weighted.mean(wing_size,value,na.rm=T))
tmeansMeans <- ddply(tmeansMeans,.(State,Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

q6<-ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean,colour=State))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI,fill=State),alpha=0.3)+
  theme_bw()+
  ylab("Wing length")+
  facet_grid(~State)+
  theme(legend.position="none")

library(cowplot)
plot_grid(q1,q2,q3,q4,q5,q6,align="v",ncol=1)

library(gridExtra)
grid.arrange(q5,q6,ncol=1)

############################################################################
