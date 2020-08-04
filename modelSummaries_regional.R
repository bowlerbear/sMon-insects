#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

plotTS_regional <- function(x){
  require(ggplot2)
  g1 <- ggplot(data=x,aes(x=Year,y=mean,group=factor(Site)))+
    geom_line(aes(x=Year,y=mean,colour=factor(Site)))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.,fill=factor(Site)),alpha=0.5)+
    facet_wrap(~Species)
  print(g1)
}

###Species list#################################################################################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

allSpecies <- read.delim("specieslist_odonata.txt",as.is=T)$Species

###Model summaries############################################################################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_sparta_regional"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_sparta_regional_nation_naturraum_adult_")

#annual tims series
annualDF <- getBUGSfitsII(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS_regional(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#46  2803



