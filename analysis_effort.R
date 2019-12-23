#analysis_effort
###########################################################################################
#get effort files
myPath <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs"
getFiles <- list.files(myPath)
effortFiles <- getFiles[grepl("effort",getFiles)]

#read in each file
library(plyr)
effortDF <- ldply(effortFiles,function(x){
  temp <- readRDS(paste(myPath,x,sep="/"))
  return(temp)
})

#get metadata information
effortDF$File <- sapply(effortDF$File,function(x)strsplit(x,"\\.rds")[[1]][1])
effortDF$State <- sapply(effortDF$File,function(x)strsplit(x,"_")[[1]][4])
effortDF$Species <- sapply(effortDF$File,function(x)strsplit(x,"_")[[1]][5])
effortDF$effortType <- sapply(effortDF$File,function(x)strsplit(x,"_")[[1]][2])

#look at parameters
unique(effortDF$Param)

#file names explained below
###########################################################################################
#effort (1): linear length and single list
bugs.data$Effort <- bugs.data[["nuSpecies"]]
saveRDS(data.frame(out$summary),file=paste0("outSummary_nuSpecies_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################
#effort (2): short and single list
bugs.data$Effort <- bugs.data[["shortList"]]
saveRDS(data.frame(out$summary),file=paste0("outSummary_speciesList_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################
#effort (3): linear length and single list, and expertise
bugs.data$Effort <- bugs.data[["nuSpecies"]]
bugs.data$Effort2 <- bugs.data[["expertise"]]
saveRDS(data.frame(out$summary),file=paste0("outSummary_expertise_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################
#effort (4): linear length and single list, and number of records
bugs.data$Effort <- bugs.data[["shortList"]]
bugs.data$Effort2 <- bugs.data[["nuRecs"]]
saveRDS(data.frame(out$summary),file=paste0("outSummary_nuRecs_",stage,"_",state,"_", myspecies,".rds"))
############################################################################################

#histograms

library(ggplot2)

#effort.p
effortDF_effort.p <- subset(effortDF,Param =="effort.p")
ggplot(effortDF_effort.p)+
  geom_histogram((aes(mean)))+
  facet_wrap(~effortType)


#single list
effortDF_single.p <- subset(effortDF,Param =="single.p")
ggplot(effortDF_single.p)+
  geom_histogram((aes(mean)))+
  facet_wrap(~effortType)


#effort2.p
effortDF_effort2.p <- subset(effortDF,Param =="effort2.p")
ggplot(effortDF_effort2.p)+
  geom_histogram((aes(mean)))+
  facet_wrap(~effortType)

#########################################################################################

#do other params depend on effort model
effortDF_phenol.p <- subset(effortDF,Param =="phenol.p")
ggplot(effortDF_phenol.p)+
  geom_histogram((aes(mean)))+
  facet_wrap(~effortType)+
  xlim(-0.2,0.2)

effortDF_phenol2.p <- subset(effortDF,Param =="phenol2.p")
ggplot(effortDF_phenol2.p)+
  geom_histogram((aes(mean)))+
  facet_wrap(~effortType)+
  xlim(-0.01,0.01)

#########################################################################################

#yearly estimates
effortDF_psi <- effortDF[grepl("psi",effortDF$Param),]
effortDF_psi$Year <- sapply(effortDF_psi$Param,function(x){
  sub(".*\\[([^][]+)].*", "\\1", x)})
effortDF_psi$Year <- as.numeric(effortDF_psi$Year)+1979

ggplot(subset(effortDF_psi,State=="NRW"))+
  geom_point(aes(x=Year,y=mean,colour=effortType),alpha=0.5)+
  facet_wrap(~Species)

ggplot(subset(effortDF_psi,State=="Sa"))+
  geom_line(aes(x=Year,y=mean,colour=effortType),alpha=0.5)+
  facet_wrap(~species)+
  scale_x_continuous(labels=c(1980,2010),breaks=c(1980,2010))+
  ylab("Occupancy probability")+
  theme_bw()+
  theme(legend.position="top")

ggplot(subset(effortDF_psi,State=="SAnhalt"))+
  geom_point(aes(x=Year,y=mean,colour=effortType),alpha=0.5)+
  facet_wrap(~Species)

#########################################################################################

#number of significant effort terms
effortDF$sig <- "none"
effortDF$sig[effortDF$X2.5.>0 & effortDF$X97.5.>0] <- "increase"
effortDF$sig[effortDF$X2.5.<0 & effortDF$X97.5.<0] <- "decrease"

out <- ddply(effortDF,.(Param,effortType),summarise,nuSigs = length(sig[sig!="none"]))

ggplot(out)+
  geom_bar(aes(x=Param,nuSigs,fill=effortType),stat="identity",position="dodge")+
  coord_flip()

#########################################################################################

#data just for one state
#effortDF <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/modelSummary_effort_tests_Odonata_adult_MV.rds")
effortDF <- readRDS("C:/Users/db40fysa/Nextcloud/sMon-Analyses/Git/sMon-insects/model-outputs/modelSummary_effort_tests_Odonata_adult_NS.rds")
effortDF$File <- sapply(effortDF$File,function(x)strsplit(x,"\\.rds")[[1]][1])
effortDF$State <- sapply(effortDF$File,function(x)strsplit(x,"_")[[1]][4])
effortDF$Species <- sapply(effortDF$File,function(x)strsplit(x,"_")[[1]][5])
effortDF$effortType <- sapply(effortDF$File,function(x)strsplit(x,"_")[[1]][2])
effortDF <- subset(effortDF, effortType=="nuSpecies")

#year estimates
effortDF_psi <- effortDF[grepl("psi",effortDF$Param),]
effortDF_psi$Year <- sapply(effortDF_psi$Param,function(x){
  sub(".*\\[([^][]+)].*", "\\1", x)})
effortDF_psi$Year <- as.numeric(as.character(effortDF_psi$Year)) + 1979

#abbreviate genus and species to 3 letters
effortDF_psi$Genus <- sapply(effortDF_psi$Species,function(x)strsplit(x," ")[[1]][1])
effortDF_psi$Genus <- sapply(effortDF_psi$Genus,function(x)substr(x,1,4))
effortDF_psi$Spec <- sapply(effortDF_psi$Species,function(x)strsplit(x," ")[[1]][2])
effortDF_psi$Spec <- sapply(effortDF_psi$Spec,function(x)substr(x,1,4))
effortDF_psi$species <- paste(effortDF_psi$Genus,effortDF_psi$Spec," ")

#plot the time series
library(ggplot2)
ggplot(effortDF_psi)+
  geom_line(aes(x=Year,y=mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  facet_wrap(~species)+
  theme_bw()+
  theme(strip.text = element_text(size = rel(0.75)))+
  scale_x_continuous(labels=c(1980,2010),breaks=c(1980,2010))+
  ylab("Occupancy probability")

#########################################################################################


