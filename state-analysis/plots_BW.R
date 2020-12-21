library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

###get model summaries#################################################################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_BW/7458244"
#fixed subsetting code, with eta, ecoregion 1

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_BW_naturraum_adult_")


#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")

#for MTBQ
annualDF$Year <- annualDF$ParamNu + 1979

#remove Anax epi
annualDF <- subset(annualDF,Species!="Anax ephippiger")

#how many years do we have for each species
ddply(annualDF,.(Species),summarise,nuYears = length(Year))
#39 for all

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~Code,ncol=6)+
  facet_wrap(~Code,ncol=6,scales = "free_y")+
  xlab("Jahr") + ylab("MTBQ-Präsenz-Anteil")+
  theme_bw()

ggsave("state-analysis/tsplot_MTBQ_BW.png",width=10,height=11)
ggsave("state-analysis/tsplot_MTBQ_BW_freey.png",width=10,height=11)

#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
#4, but only slightly

write.csv(trendsDF,file="state-analysis/trends_MTBQ.rds")
write.csv(trendsDF,file="state-analysis/trends_MTB64.rds")

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
nrow(trendEstimates)#61

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

### Red list plot ####################################################

# table(trendEstimates$GermanRedList)
# 
# #or add numbers of species to axes label??
# trendEstimates$GermanRedListL <- factor(trendEstimates$GermanRedList)
# levels(trendEstimates$GermanRedListL)
# trendEstimates$GermanRedListL <- factor(trendEstimates$GermanRedListL,
#                                         levels=c("*","V","R","3","2","1"))
# 
# levels(trendEstimates$GermanRedListL) <- 
#   c("* (42)", "V (5)",
#     "R (1)","3 (8)",
#     "2 (4)","1 (1)")
# 
# g1 <- ggplot(trendEstimates)+
#   geom_boxplot(aes(x=GermanRedListL,y=trend))+
#   geom_hline(yintercept=0,linetype="dashed")+
#   coord_flip()+
#   theme_classic()+
#   ylab("Langzeitveränderung")+xlab("Rote-Liste Kategorie Deutschlands")


### winner/loser #####################################################

trendEstimates$Species <- factor(trendEstimates$Species,levels=c(trendEstimates$Species[order(trendEstimates$trend)]))

trendEstimates$Trend <- factor(trendEstimates$Trend,
                               levels=c("significant decrease",
                                        "stable",
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

ggsave("state-analysis/mittlere_jährliche_Veränderung_MTB64.png",width=8,height=6)

### total change ####################################################


subset(annualDF,Year==1980)
nrow(subset(annualDF,Year==2019))
subset(annualDF,Year==2018)
unique(annualDF$Year)

totalChange <- ddply(annualDF,.(Species),summarise,
                     change=mean[Year==2019]/mean[Year==1981],
                     growth=((mean[Year==2019]/mean[Year==1981])-1)*100)
totalChange$change
summary(totalChange$change)

totalChange <- arrange(totalChange,change)

totalChange$Species[!totalChange$Species %in% trendEstimates$Species]
trendEstimates$Species[!trendEstimates$Species %in% totalChange$Species]

trendEstimates$change <- totalChange$change[match(trendEstimates$Species,totalChange$Species)]

trendEstimates$Species <- factor(trendEstimates$Species,
                                 levels=c(totalChange$Species[order(totalChange$change)]))

summary(totalChange$change)

g3 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=change),
           stat="identity",width=rel(1))+
  theme_classic()+  
  scale_y_log10(breaks=c(0,0.01,0.1,1,10,100),
                labels=c(0,0.01,0.1,1,10,100))+
  theme(axis.text.x = element_text(angle=90, size=7,hjust=0.95,vjust=0.2))+
  ylab("Verhältnis der Präsenz zwischen 2019 und 1980")+xlab("Arten")

ggsave("state-analysis/Verhältnis_MTBQ.png",width=8,height=6)

#combine together
### end ##############################################################