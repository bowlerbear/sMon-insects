library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

###get model summaries#################################################################################
### MTBQ ##########################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_Hessen/6826255"
#fixed subsetting code, with eta, ecoregion 1, simple initial values - works for all!!

#get updated models with more recent data
mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_He/MTBQ"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_He_naturraum_adult_")

#remove species with large uncertainty
#modelDF <- subset(modelDF,!Species %in% c("Leucorrhinia albifrons"))

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF <- arrange(annualDF,Species,ParamNu)

#for MTBQ
annualDF$Year <- annualDF$ParamNu + 1979
#it depends on the number of years
annualDF <- ddply(annualDF,.(Species),function(x){
  len <- length(x$Year)
  if(len==39){
    x$Year <- 1981:2019
  }else if(len==38){
    x$Year <- c(1981:1995,1997:2019)
  }
  return(x)
})


#how many years do we have for each species
ddply(annualDF,.(Species),summarise,nuYears = length(Year))
#lowest is 38 for Anax ephippiger

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6)+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  xlab("Jahr") + ylab("MTBQ-Präsenz-Anteil")+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTBQ_Hessen_trends.png",width=10,height=9)

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6,scales = "free_y")+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  xlab("Jahr") + ylab("MTBQ-Präsenz-Anteil")+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTBQ_Hessen_freey_trends.png",width=10,height=9)


#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
#4, but only slightly

#add values of first and last year
trendsDF$annual1980 <- annualDF$mean[annualDF$Year==1981]
trendsDF$annual2019 <- annualDF$mean[annualDF$Year==2019]

write.csv(trendsDF,file="state-analysis/Hessen/trends_MTBQ.csv",row.names=FALSE)

### MTB64 #################################################################

#get models at the scale of MTB64 - 2016 onwards
mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_He/MTB64"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_He_MTB64_naturraum_adult_")

#annual time series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF <- arrange(annualDF,Species,ParamNu)

#for MTB64
annualDF$Year <- annualDF$ParamNu + 2005
plotTS(annualDF)
table(annualDF$Rhat<1.1)
annualDF$Rhat[annualDF$Rhat>1.1]
#FALSE  TRUE 
#14   840

#how many years do we have for each species
ddply(annualDF,.(Species),summarise,nuYears = length(Year))
#lowest is 14

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6)+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  xlab("Jahr") + ylab("MTB64-Präsenz-Anteil")+
  scale_x_continuous(breaks=c(2006,2011,2016),labels=c(2006,2011,2016))+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTB64_Hessen_trends.png",width=9,height=9)

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6,scales = "free_y")+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  xlab("Jahr") + ylab("MTB64-Präsenz-Anteil")+
  scale_x_continuous(breaks=c(2006,2011,2016),labels=c(2006,2011,2016))+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTB64_Hessen_freey_trends.png",width=9,height=9)

#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
#4, but only slightly

#add values of first and last year
trendsDF$annual2006 <- annualDF$mean[annualDF$Year==2006]
trendsDF$annual2019 <- annualDF$mean[annualDF$Year==2019]

#add values of first and last year
write.csv(trendsDF,file="state-analysis/Hessen/trends_MTB64.csv",row.names = FALSE)

### trends comparison ########################################################

allTrendsDF <- merge(trendsDF_MTB64,trendsDF_MTBQ,by=c("Species","Code"))
ggplot(allTrendsDF)+
  geom_text(aes(x = mean.x,
                 y= mean.y,
                label= Code),size=2.5)+
  xlab("MTB64 trend (2006-)")+ylab("MTBQ trend(1980-)")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme_bw()
        
ggsave("state-analysis/trend_scale_comparison.png",width=6,height=5)

###merge trait data############################################################

setwd("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects")

#get traits data
load("alltraits.RData")

trendEstimates <- trendsDF

names(trendEstimates)[which(names(trendEstimates)=="mean")]<-"trend"
trendEstimates$Species[!trendEstimates$Species %in% alltraits$Species]

trendEstimates <- merge(trendEstimates,alltraits,by="Species",all.x=T)
nrow(trendEstimates)#63

###trend classify #############################################################

summary(trendEstimates$trend)
hist(trendEstimates$trend)#pretty normal

#null
trendEstimates$Trend <- "insignificant"

#trend is annual proportional change in number of sites
#mean 1% change over the 35 year - less than total change 5%
nu <- 0.05/ length(unique(annualDF$Year))

#standard deviation to large to detect a 5% change

#stable
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.>0
                     & abs(trendEstimates$trend) < nu]<-"stable"

#increasing
trendEstimates$Trend[trendEstimates$X2.5.>0 & trendEstimates$X97.5.>0]<-"significant increase"

#decreasing
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.<0]<-"significant decrease"

#summary
table(trendEstimates$Trend)

### winner/loser #####################################################

trendEstimates$Species <- factor(trendEstimates$Species,levels=c(trendEstimates$Species[order(trendEstimates$trend)]))
trendEstimates$Trend[trendEstimates$Trend=="stable"] <- "insignificant"
trendEstimates$Trend <- factor(trendEstimates$Trend)
levels(trendEstimates$Trend)

trendEstimates$Trend <- factor(trendEstimates$Trend,
                               levels=c("significant decrease",
                                        #"stable",
                                        "insignificant",
                                        "significant increase"))

g2 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=trend,fill=Trend),
           stat="identity",width=rel(1))+
  theme_classic()+  
  theme(axis.text.x = element_text(angle=90, size=7,hjust=0.95,vjust=0.2))+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  ylab("mittlere jährliche Veränderung")+xlab("Arten")

### total change ####################################################

subset(annualDF,Year==1980)
nrow(subset(annualDF,Year==2019))
subset(annualDF,Year==2018)
unique(annualDF$Year)

totalChange <- ddply(annualDF,.(Species),summarise,
                     change=mean[Year==2019]/mean[Year==1981],
                     perChange = (((mean[Year==2019]-mean[Year==1981])/mean[Year==1981])*100),
                     growth=((mean[Year==2019]/mean[Year==1981])-1)*100)

totalChange <- ddply(annualDF,.(Species),summarise,
                     change=mean[Year==2019]/mean[Year==2006],
                     perChange = (((mean[Year==2019]-mean[Year==2006])/mean[Year==2006])*100),
                     growth=((mean[Year==2019]/mean[Year==2006])-1)*100)

totalChange$change
summary(totalChange$change)
summary(totalChange$perChange)
summary(totalChange$growth)

#order
totalChange <- arrange(totalChange,change)
totalChange$perChange[totalChange$perChange>1000] <- 500

#check all there
totalChange$Species[!totalChange$Species %in% trendEstimates$Species]
trendEstimates$Species[!trendEstimates$Species %in% totalChange$Species]

#reorder
trendEstimates$change <- totalChange$change[match(trendEstimates$Species,totalChange$Species)]
trendEstimates$Species <- factor(trendEstimates$Species,
                                 levels=c(totalChange$Species[order(totalChange$change)]))

g3 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=change),
           stat="identity",width=rel(1))+
  theme_classic()+  
  scale_y_log10(breaks = c(0.01,0.1,1,10,100), labels=c(0.01,0.1,1,10,100))+
  theme(axis.text.x = element_text(angle=90, size=7,hjust=0.95,vjust=0.2))+
  ylab("Verhältnis (2019/1980)")+xlab("Arten")
  #ylab("Verhältnis (2019/2006)")+xlab("Arten")

### save together #################################

cowplot::plot_grid(g2,g3,ncol=1,labels=c("A","B"))

ggsave("state-analysis/Hessen/Verhältnis_MTBQ.png",width=7,height=8)

ggsave("state-analysis/Hessen/Verhältnis_MTB64.png",width=7,height=8)

### regression line predictions for each year #####

#so for each year, you want the predicted y along the regression line? Are you sure you want this? As seen in the case of Gomplus flavipes, the annual values can lead to strange when species changes are not linear (overestimates in the 1980s and 2010s). I regard the regression line as just a summary of the change over the study period.

#if you do want it: the column labelled "Trend" already contains the slope of the regression line but I can also give you the intercept of the line, the "a" in y = a + bx to calculate the annual points ("b" is the Trend).

#yes, exactly - the regression takes into account the uncertainty of the estimate each year so years with high uncertainty - usually lower amount of data - have lower weight on the regression

# The paper is printed in DIN-A5 format with 300 dpi. With a page width of 120mm, this results in an image width of 1440 pixels. This means that the current size of the individual diagrams is already OK, 

#But they are comparatively "coarse" and a somewhat smaller font size would be advantageous. 

#I would like to set up the diagrams in three columns and compare your results with the calculated grid frequencies according to the attached pattern. It would be good if we could agree on font type and font size. 
#I have used Calibri 20pt as a sample, but I am open. However, the font should not be much larger. 
# Your plots are very clear as a compilation, but in this form they are not suitable for regrouping. 
#Individual diagrams, each with complete axis labelling, would be more advantageous! 
#Probably easiest if I add on the final graph column - "RFreq in %" and create the 3 graphs in the same format at the same time.
#RFreq in %" is simply the % of MTBQs in which the species is recorded each year?

### MTBQ #############################################################

#get trends based on annualDF and simple lm

annualDF_trends <- ddply(annualDF,.(Species,Code),function(x){
  lm1 <- lm(mean ~ Year, weights = 1/sd, data = x)
  summary(lm1)$coef[2,]
})

#compare with 
allTrends <- merge(trendsDF,annualDF_trends,by=c("Species","Code"))
qplot(mean,Estimate,data=allTrends)

#add to predicted values to the dataset
annualDF_regressionFits <- ddply(annualDF,.(Species,Code),function(x){
  lm1 <- lm(mean ~ Year, weights = 1/sd, data = x)
  myYear <- 1980:2019
  mySpecies <- unique(x$Species)
  predValues <- predict(lm1,newdata=data.frame(Year = myYear))
  data.frame(Species=mySpecies,
             Year=myYear,
             Regression_prediction=predValues)
})

write.csv(annualDF_regressionFits,file="state-analysis/Hessen/Regression_fits_MTBQ.csv",row.names=FALSE)

### MTB64 ###########################################################

#as above

annualDF_regressionFits <- ddply(annualDF,.(Species,Code),function(x){
  lm1 <- lm(mean ~ Year, weights = 1/sd, data = x)
  myYear <- 2006:2019
  mySpecies <- unique(x$Species)
  predValues <- predict(lm1,newdata=data.frame(Year = myYear))
  data.frame(Species=mySpecies,
             Year=myYear,
             Regression_prediction=predValues)
})

write.csv(annualDF_regressionFits,file="state-analysis/Hessen/Regression_fits_MTB64.csv",row.names=FALSE)

### RD data ##########################################################

#MTBQ
annualDF_MTBQ <- annualDF
#MT64
annualDF_MTB64 <- annualDF
#RF
annualDF_RF <- read.csv("state-analysis/Hessen/Rasterfrequenz_A_gesamt_DB.csv",sep=";")
annualDF_RF <- annualDF_RF[,c(1,8:20)] 
annualDF_RF <- melt(annualDF_RF,id="Rasterfrequenz")
names(annualDF_RF) <- c("Species","Year","RF")
annualDF_RF$Year <- as.numeric(gsub("X","",annualDF_RF$Year))

#draw plot for each species
library(ggthemes)
species_MTBQ <- sort(unique(annualDF_MTBQ$Species))
species_MTB64 <- sort(unique(annualDF_MTB64$Species))
species_MTBQ[!species_MTBQ %in% species_MTB64]
#"Sympetrum depressiusculum" "Sympetrum pedemontanum"
species_MTB64[!species_MTB64 %in% species_MTBQ]
#0
species_RF <- sort(unique(annualDF_RF$Species))
species_RF[!species_RF %in% species_MTB64]
#0

speciesList <- sort(unique(species_MTB64,species_RF))

for(i in 1:length(speciesList)){
  

mySpecies <- speciesList[i]

#same scale for the occupancy models 

#find the min 
temp <- rbind(subset(annualDF_MTBQ,Species==mySpecies),
              subset(annualDF_MTB64,Species==mySpecies))
minY <- min(temp$X2.5.)*100
maxY <- max(temp$X97.5.)*100

g1 <- ggplot(subset(annualDF_MTBQ,Species==mySpecies),aes(x=Year,y=mean*100))+
  geom_line()+
  geom_ribbon(aes(x=Year,ymin=X2.5.*100,ymax=X97.5.*100),alpha=0.15)+
  stat_smooth(method="lm",se=F,linetype="dashed")+
  theme_few()+
  ylim(minY,maxY)+
  #ylim(0,100)+
  ylab("Besetzungsgrad (MTB-Q) %")+xlab("Jahr")+
  ggtitle(mySpecies)

g2 <- ggplot(subset(annualDF_MTB64,Species==mySpecies),aes(x=Year,y=mean*100))+
  geom_line()+
  geom_ribbon(aes(x=Year,ymin=X2.5.*100,ymax=X97.5.*100),alpha=0.15)+
  stat_smooth(method="lm",se=F,linetype="dashed")+
  theme_few()+
  ylim(minY,maxY)+
  #ylim(0,100)+
  ylab("Besetzungsgrad (MTB/64) %")+xlab("Jahr")+
  scale_x_continuous(breaks=c(2005,2010,2015,2020),labels=c(2005,2010,2015,2020),limits=c(2005,2020))+
  ggtitle("")

g3 <- ggplot(subset(annualDF_RF,Species==mySpecies),aes(x=Year,y=RF))+
  geom_line()+
  theme_few()+
  #ylim(0,100)+
  ylab("Rasterfrequenz (MTB/64) %")+xlab("Jahr")+
  stat_smooth(method="lm",se=F,linetype="dashed")+
  scale_x_continuous(breaks=c(2005,2010,2015,2020),labels=c(2005,2010,2015,2020),limits=c(2005,2020))+
  ggtitle("")

cowplot::plot_grid(g1,g2,g3,ncol=3)

ggsave(paste0("state-analysis/Hessen/speciesPlots/",mySpecies,".png"),width=9,height=3)

}



### end ##############################################################
