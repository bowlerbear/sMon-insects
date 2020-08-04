library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')


###get model summaries#################################################################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_Hessen/6826255"
#fixed subsetting code, with eta, ecoregion 1, simple initial values - works for all!!

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_Hessen_adult_")


#remove species with large uncertainty
modelDF <- subset(modelDF,!Species %in% c("Leucorrhinia albifrons"))

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#44  2189

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Species,ncol=6)+
  xlab("Jahr") + ylab("MTBQ-Präsenz-Anteil")+
  theme_bw()

#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
#4, but only slightly

#some species disappear after 2014?????

###merge trait data############################################################

setwd("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects")

#get traits data
load("alltraits.RData")

trendEstimates <- trendsDF

names(trendEstimates)[which(names(trendEstimates)=="mean")]<-"trend"
trendEstimates$Species[!trendEstimates$Species %in% alltraits$Species]


trendEstimates <- merge(trendEstimates,alltraits,by="Species",all.x=T)
nrow(trendEstimates)#61


###trend classify #############################################################

summary(trendEstimates$trend)
hist(trendEstimates$trend)#pretty normal

#null
trendEstimates$Trend <- "insignificant"

#trend is annual proportional change in number of sites
#mean 1% change over the 35 year - less than total change 5%
#0.05/36
#[1] 0.001388889

#standard deviation to large to detect a 5% change

#stable
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.>0
                     & abs(trendEstimates$trend) < 0.001388889]<-"stable"

#increasing
trendEstimates$Trend[trendEstimates$X2.5.>0 & trendEstimates$X97.5.>0]<-"significant increase"

#decreasing
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.<0]<-"significant decrease"

#summary
#table(trendEstimates$Trend)
#significant decrease significant increase               stable 
#19                   16                   20

### Red list plot ####################################################

table(trendEstimates$GermanRedList)

#or add numbers of species to axes label??
trendEstimates$GermanRedListL <- factor(trendEstimates$GermanRedList)
levels(trendEstimates$GermanRedListL)
trendEstimates$GermanRedListL <- factor(trendEstimates$GermanRedListL,
                                        levels=c("*","V","R","3","2","1"))

levels(trendEstimates$GermanRedListL) <- 
  c("* (42)", "V (5)",
    "R (1)","3 (8)",
    "2 (4)","1 (1)")

g1 <- ggplot(trendEstimates)+
  geom_boxplot(aes(x=GermanRedListL,y=trend))+
  geom_hline(yintercept=0,linetype="dashed")+
  coord_flip()+
  theme_classic()+
  ylab("Langzeitveränderung")+xlab("Rote-Liste Kategorie Deutschlands")


### winner/loser #####################################################

trendEstimates$Species <- factor(trendEstimates$Species,levels=c(trendEstimates$Species[order(trendEstimates$trend)]))
trendEstimates$Trend[trendEstimates$Trend=="stable"] <- "insignificant"
trendEstimates$Trend <- factor(trendEstimates$Trend,levels=c("significant decrease","insignificant","significant increase")) 

g2 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=trend,fill=Trend),
           stat="identity",width=rel(1))+
  theme_classic()+  
  theme(axis.text.x = element_text(angle=90, size=7,hjust=0.95,vjust=0.2))+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  ylab("Langzeitveränderung")+xlab("Arten")

### total change ####################################################

totalChange <- ddply(annualDF,.(Species),summarise,
                     change=mean[Year==2014]/mean[Year==1980],
                     growth=((mean[Year==2014]/mean[Year==1980])-1)*100)

summary(totalChange$change)

totalChange <- arrange(totalChange,change)

trendEstimates$change <- totalChange$change[match(trendEstimates$Species,totalChange$Species)]

trendEstimates$Species <- factor(trendEstimates$Species,
                                 levels=c(totalChange$Species[order(totalChange$change)]))


g3 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=change),
           stat="identity",width=rel(1))+
  theme_classic()+  
  scale_y_log10(breaks=c(0,0.01,0.1,1,10,100,100),
                labels=c(0,0.01,0.1,1,10,100,100))+
  theme(axis.text.x = element_text(angle=90, size=7,hjust=0.95,vjust=0.2))+
  ylab("Verhältnis der Präsenz zwischen 2016 und 1980")+xlab("Arten")

### end ##############################################################