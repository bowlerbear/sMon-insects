library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

#models 2018 onwards and including mid naturraum as a random effect

###get model summaries#################################################################################
### MTBQ ##########################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_He/MTBQ_MidNatur"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_He_naturraum_adult_")

#remove species with large uncertainty
#modelDF <- subset(modelDF,!Species %in% c("Leucorrhinia albifrons"))

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF <- arrange(annualDF,Species,ParamNu)
table(annualDF$Rhat<1.1)
annualDF$Rhat[annualDF$Rhat>1.1]
#FALSE  TRUE 
#54   678

#for MTBQ
annualDF$Year <- annualDF$ParamNu + 2007

#how many years do we have for each species
ddply(annualDF,.(Species),summarise,nuYears = length(Year))
#12 for all

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6)+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  scale_x_continuous(breaks=c(2008,2012,2016),labels=c(2008,2012,2016))+
  xlab("Jahr") + ylab("MTBQ-Präsenz-Anteil")+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTBQ_Hessen_trends_updated.png",width=10,height=9)

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6,scales = "free_y")+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  scale_x_continuous(breaks=c(2008,2012,2016),labels=c(2008,2012,2016))+
  xlab("Jahr") + ylab("MTBQ-Präsenz-Anteil")+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTBQ_Hessen_freey_trends_updated.png",width=10,height=9)


#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
#4, but only slightly

#add values of first and last year
trendsDF$annual1980 <- annualDF$mean[annualDF$Year==2008]
trendsDF$annual2019 <- annualDF$mean[annualDF$Year==2019]

write.csv(trendsDF,file="state-analysis/Hessen/trends_MTBQ_updated.csv",row.names=FALSE)

### MTB64 #################################################################

#get models at the scale of MTB64 - 2016 onwards
mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_He/MTB64_MidNatur"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_He_MTB64_naturraum_adult_")

#annual time series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF <- arrange(annualDF,Species,ParamNu)

#for MTB64
annualDF$Year <- annualDF$ParamNu + 2007
table(annualDF$Rhat<1.1)
annualDF$Rhat[annualDF$Rhat>1.1]
#FALSE  TRUE 
#132   840
annualDF$Species[annualDF$Rhat>1.5]

#remove problematic species
probSpecies <- c("Aeshna juncea","Coenagrion hastulatum","Somatochlora flavomaculata")
annualDF <- subset(annualDF, !Species %in% probSpecies)

#how many years do we have for each species
ddply(annualDF,.(Species),summarise,nuYears = length(Year))
#all are 12

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6)+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  xlab("Jahr") + ylab("MTB64-Präsenz-Anteil")+
  scale_x_continuous(breaks=c(2006,2011,2016),labels=c(2008,2012,2016))+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTB64_Hessen_trends_updated.png",width=9,height=9)

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.3)+
  facet_wrap(~Code,ncol=6,scales = "free_y")+
  stat_smooth(aes(x=Year,y=mean),method="lm",se=FALSE,linetype="dashed")+
  xlab("Jahr") + ylab("MTB64-Präsenz-Anteil")+
  scale_x_continuous(breaks=c(2006,2011,2016),labels=c(2008,2012,2016))+
  theme_bw()

ggsave("state-analysis/Hessen/tsplot_MTB64_Hessen_freey_trends_updated.png",width=9,height=9)

#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
#4, but only slightly
trendsDF <- subset(trendsDF, !Species %in% probSpecies)

#add values of first and last year
trendsDF$annual2006 <- annualDF$mean[annualDF$Year==2008]
trendsDF$annual2019 <- annualDF$mean[annualDF$Year==2019]

#add values of first and last year
write.csv(trendsDF,file="state-analysis/Hessen/trends_MTB64_updated.csv",row.names = FALSE)

### end ##############################################################
