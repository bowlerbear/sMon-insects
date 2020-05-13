#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

###Species list#################################################################################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

allSpecies <- read.delim("specieslist_odonata.txt",as.is=T)$Species

###Model summaries############################################################################################

#model summaries:
#modelfile="/data/idiv_ess/Odonata/BUGS_dynamic_nation_naturraum_raumFEyear1_rw1.txt"
#single random walk

#modelfile="/data/idiv_ess/Odonata/BUGS_dynamic_nation_naturraum_raumFEyear1_rw.txt"
#double random walk

###Nation state model#########################################################################################

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_state"

#read in model summaries
modelDF <- getModelSummaries(mdir)

#species code
modelDF$Species <- gsub("out_dynamic_nation_state_adult_","",modelDF$File)
modelDF$Species <- gsub(".rds","",modelDF$Species)
modelDF$Genus <- sapply(modelDF$Species, function(x)strsplit(x," ")[[1]][1])
modelDF$Genus <- sapply(modelDF$Genus, function(x) substr(x,1,3))
modelDF$Spec <- sapply(modelDF$Species, function(x)strsplit(x," ")[[1]][2])
modelDF$Spec <- sapply(modelDF$Spec, function(x) substr(x,1,3))
modelDF$Code <- paste(modelDF$Genus,modelDF$Spec,sep="_")

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979

#plot
plotTS(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#1350  1351

#where is not converging!!
table(annualDF$Rhat<1.1,annualDF$Year)
table(annualDF$Rhat<1.1,annualDF$Species)#variable...fits for some not for others


#change between first and last year
changeDF <- ddply(annualDF,.(Species),summarise,
                  change=(mean[Year==2016]-mean[Year==1980])/mean[Year==1980],
                  change2=(mean[Year==2016]-mean[Year==1985])/mean[Year==1985])

summary(changeDF$change[changeDF$Species %in% trendEstimates$Species[trendEstimates$Trend=="significant increase"]])
summary(changeDF$change2[changeDF$Species %in% trendEstimates$Species[trendEstimates$Trend=="significant increase"]])
summary(changeDF$change[changeDF$Species %in% trendEstimates$Species[trendEstimates$Trend=="significant decrease"]])
summary(changeDF$change2[changeDF$Species %in% trendEstimates$Species[trendEstimates$Trend=="significant decrease"]])

#cut first few years
annualDF <- subset(annualDF, Year >=1982)

#traceplots
setwd(mdir)
out <- readRDS("out_dynamic_nation_state_adult_Aeshna affinis.rds")
tracePlot(out)
out <- readRDS("out_dynamic_nation_state_adult_Sympetrum danae.rds")
tracePlot(out)
out <- readRDS("out_dynamic_nation_state_adult_Anax imperator.rds")
tracePlot(out)

#maybe include random walk priors to help converge...?
ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean*100))+
  geom_ribbon(aes(x=Year,ymin=X2.5.*100,ymax=X97.5.*100),alpha=0.5)+
  facet_wrap(~Code)+
  theme_bw()+
  ylab("% occupancy")

#annual DF for each state
annualDF <- getBUGSfitsII(modelDF,param="psi.state")
annualDF$Year <- annualDF$Year + 1979

load("stateIndices.RData")
annualDF$State <- stateIndices$State[match(annualDF$State,stateIndices$stateIndex)]

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean*100,colour=State))+
  #geom_ribbon(aes(x=Year,ymin=X2.5.*100,ymax=X97.5.*100,fill=factor(State)),alpha=0.5)+
  facet_wrap(~Code)+
  ylim(0,100)+
  theme_bw()+
  theme(legend.position="top")+
  ylab("% Occupancy")

####Nation trends############################################################### 

#see analysis_HPC_trends

trends <- readRDS("model-outputs/modelTrends_nation_state_trends.rds")

trends$Species <- gsub("out_dynamic_nation_state_adult_","",
                       trends$file)
trends$Species <- gsub(".rds","",
                       trends$Species)
trends$Trend <- "insignificant"
trends$Trend[trends$X2.5.>0 & trends$X97.5.>0]<-"significant increase"
trends$Trend[trends$X2.5.<0 & trends$X97.5.<0]<-"significant decrease"
table(trends$Trend)
trends$Trend <- factor(trends$Trend,levels=c("significant decrease",
                                             "insignificant",
                                             "significant increase"))

table(trends$Trend)
#significant decrease        insignificant significant increase 
#14                   32                   27 

library(ggplot2)
g1 <- ggplot(trends)+
  geom_histogram(aes(x=mean,fill=Trend))+
  theme_bw()+
  xlab("trend estimate")+
  ylab("number of species")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_fill_viridis_d()


#match Davids
library(RColorBrewer)
blues_fun <- colorRampPalette(brewer.pal(11, "RdYlGn"))
mycols <- blues_fun(nrow(trends))
#plot as freq poly
trends$order <- rank(trends$mean)
g1 <- ggplot(aes(x=order,y=mean),data=trends)+
  geom_bar(aes(fill=mean),stat="identity",width=rel(1))+
  geom_hline(yintercept=0,colour='black')+
  scale_x_reverse()+
  scale_fill_continuous(name="Magnitude",low="#A50026",high="#006837")+
  xlab("Species")+
  ylab("Population trend")+
  theme_classic()+
  theme(axis.text.x = element_blank())
#theme(legend.position = "top")

#which species are declining
subset(trends,Trend=="significant decrease")
#Sympetrum, Coenagrion

#Zissou
library(wesanderson)
mycols <- wes_palette("Zissou1", nrow(trends), type = "continuous")
trends$order <- rank(trends$mean)

g1 <- ggplot(aes(x=order,y=mean),data=trends)+
  geom_col(aes(fill=as.factor(mean)),width=rel(1))+
  geom_hline(yintercept=0,colour='black')+
  scale_x_reverse()+
  scale_fill_manual(values=mycols)+
  xlab("Species")+
  ylab("Population trend")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        legend.position = "none")

####Naturraum models######################################################################

#naturraum analysis

#original models
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum/5914536"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF$Species <- gsub("out_dynamic_nation_naturraum_adult_","",modelDF$File)
modelDF$Species <- gsub(".rds","",modelDF$Species)

#national time series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
table(annualDF$Rhat<1.1)

#annual time series
annualDF <- getBUGSfitsII(modelDF,param="psi.raum")
annualDF$Year <- annualDF$Year + 1979

#convergence
table(annualDF$Rhat<1.1)#better!!
ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean,colour=factor(Species)))+
  facet_wrap(~State)+
  theme(legend.position="none")#each has a range

ggplot(annualDF)+
  geom_line(aes(x=Year,y=mean,colour=factor(State)))+
  facet_wrap(~Species)+
  theme(legend.position="none")

#plot as time-series
annualnaturraumDF <- ddply(annualDF,.(Year,State),summarise,nuSpecies=sum(mean))
regions <- c("Alps",
             "Alps foreland",
             "NE lowlands",
             "E uplands",
             "SW uplands",
             "W uplands",
             "NW lowlands")

annualnaturraumDF$naturraum <- as.factor(annualnaturraumDF$State)
levels(annualnaturraumDF$naturraum) <- regions
ggplot(subset(annualnaturraumDF,Year>=1985))+
  geom_line(aes(x=Year,y=nuSpecies,colour=naturraum))+
  facet_wrap(~naturraum)

ggplot(subset(annualnaturraumDF,Year>=1985))+
  geom_line(aes(x=Year,y=nuSpecies,colour=naturraum))

#do same for small-scale naturraum

###Naturraum trends######################################################

trends <- readRDS("model-outputs/modelTrends_naturraum_trends.rds")

trends$Trend <- "insignificant"
trends$Trend[trends$X2.5.>0 & trends$X97.5.>0]<-"significant increase"
trends$Trend[trends$X2.5.<0 & trends$X97.5.<0]<-"significant decrease"
table(trends$Trend)

#add species data
trends$Species <- gsub("out_dynamic_nation_naturraum_adult_","",trends$file)
trends$Species <- gsub(".rds","",trends$Species)

#no increases in 1.
#largest change in
regions <- c("Alps","Alps foreland","NE lowlands","E uplands","SW uplands",
             "W uplands","NW lowlands")

trends$naturraum <- as.factor(trends$naturraum)
levels(trends$naturraum) <- regions

library(ggplot2)
library(wesanderson)
ggplot(trends)+
  geom_histogram(aes(x=mean,fill=Trend))+
  theme_bw()+
  xlab("trend")+
  facet_wrap(~naturraum,scales="free_y")+
  ylab("number of species")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)[c(3,5,4)])

ggplot(trends)+
  geom_histogram(aes(x=mean,fill=Trend))+
  theme_bw()+
  xlab("trend")+
  facet_grid(Trend~naturraum)+
  ylab("number of species")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)[c(3,5,4)])

#plot as dotplot
ggplot(trends)+
  geom_point(aes(x=Species,y=mean))+
  facet_wrap(~naturraum,nrow=1)+
  coord_flip()+
  geom_hline(yintercept=0,colour="red")
#order species by mean trend, graph ugly anyhow

ggplot(trends)+
  geom_histogram(aes(x=mean,fill=Trend))+
  theme_bw()+
  xlab("trend")+
  ylab("number of species")+
  geom_vline(xintercept=0,linetype="dashed")+
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 5)[c(3,5,4)])


#get number of species increasing/decreasing in each region
summaryTrends <- ddply(trends,.(Trend,naturraum),summarise,nuSpecies=length(unique(Species)))
ggplot(subset(summaryTrends,Trend!="insignificant"))+
  geom_bar(aes(x=naturraum,y=nuSpecies,fill=Trend),position="dodge",stat="identity")+
  coord_flip()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 5)[c(1,3)])+
  theme_bw()

###Naturraum percol####################################################################

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_5000iter/6057269"

#read in all files
modelDF <- getModelSummaries(mdir)
modelDF$Species <- gsub("out_dynamic_nation_naturraum_adult_","",modelDF$File)
modelDF$Species <- gsub(".rds","",modelDF$Species)

#compare mean colonization and persistence rates
modelDF <- subset(modelDF, Param %in% c("meanPersist","meanColonize"))
colper <- dcast(modelDF,Species~Param,value.var="mean")
colper$trend <- trends$mean[match(colper$Species,trends$Species)]
colper$TrendSig <- trends$Trend[match(colper$Species,trends$Species)]

g2 <- ggplot(colper,aes(x=meanColonize,y=meanPersist,color=TrendSig))+
  geom_point()+
  scale_x_log10()+
  scale_colour_viridis_d()+
  xlab("colonization probability")+
  ylab("persistence probability")+
  theme_bw()+
  theme(legend.position = "none")

###Nation spline#######################################################################

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_spline/5952992"

#read in all files
modelDF <- getModels(mdir)
modelDF$Species <- gsub("outSummary_dynamicspline_adult_","",modelDF$File)
modelDF$Species <- gsub(".rds","",modelDF$Species)

#compare mean colonization and persistence rates
modelDF$Param <- sapply(modelDF$Param,function(x)ifelse(grepl("persist",x),"persist","colonise"))
modelDF <- ddply(modelDF,.(Param,Species),summarise,meanVal=mean(mean,na.rm=T))

ggplot(modelDF,aes(x=Species,y=meanVal))+
  geom_point()+
  coord_flip()+
  facet_wrap(~Param,scales="free")

#compare colonization and persistence
library(reshape2)
colper <- dcast(modelDF,Species~Param,value.var="meanVal")
colper$trend <- trends$mean[match(colper$Species,trends$Species)]
colper$TrendSig <- trends$Trend[match(colper$Species,trends$Species)]

g2 <- ggplot(colper,aes(x=colonise/7,y=persist,color=TrendSig))+
  geom_point()+
  scale_x_log10()+
  scale_colour_viridis_d()+
  xlab("colonization probability")+
  ylab("persistence probability")+
  theme_bw()+
  theme(legend.position = "none")

library(cowplot)
plot_grid(g1,g2)

#species vary more in colonisation dyanmics
cor.test(colper$trend,colper$colonise)
cor.test(colper$trend,colper$persist)

###Range-expanding species###################################################################

modelSummary <- data.frame(readRDS(paste(mdir,"outSummary_dynamicspline_adult_Crocothemis erythraea.rds",sep="/")))
modelSummary$Param <- row.names(modelSummary)

#colonize
modelSummaryP <- subset(modelSummary,grepl("colonize",modelSummary$Param))
modelSummaryP$Index <- sapply(modelSummaryP$Param,function(x){sub(".*\\[([^][]+)].*", "\\1", x)})
modelSummaryP$Index <- as.numeric(modelSummaryP$Index)

#add x and y data
load("siteInfo.RData")
length(unique(siteInfo$siteIndex))
nrow(modelSummaryP)
modelSummaryP$MTB_Q <- siteInfo$MTB_Q[match(modelSummaryP$Index,siteInfo$siteIndex)]
summary(modelSummaryP)
sum(is.na(modelSummaryP$MTB_Q))

#add to mtbqd data frame
load("mtbqsDF.RData")
head(mtbqsDF)
mtbqsDF$colonize <- modelSummaryP$mean[match(mtbqsDF$MTB_Q,modelSummaryP$MTB_Q)]

#plotting
ggplot(mtbqsDF)+
  geom_point(aes(x=x,y=y,colour=colonize))+
  scale_color_gradient2(low="blue",high='red',midpoint=median(mtbqsDF$colonize,na.rm=T))

###Problematic species check#####################################################

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_problemSpecies/Odonata_adult_nation_naturraum/6280560"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_dynamic_nation_naturraum_adult_")
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS(annualDF)
table(annualDF$Rhat<1.1)

#trend estimates
trendDF <- getBUGSfits(modelDF,param="regres.psi")
trendDF[,c("Species","mean","X2.5.","X97.5.")]
table(trendDF$Rhat<1.1)

#colonize/persist
persistDF <- getBUGSfits(modelDF,param="meanPersist")
colonizeDF <- getBUGSfits(modelDF,param="meanColonize")
table(persistDF$Rhat<1.1)
table(colonizeDF$Rhat<1.1)

#still problematic

###Sparta models##########################################################

a = 10


#run on R server -  original sparta package function
# mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/sparta-outputs"
# 
# modelDF <- getSpartaModels(mdir)
# modelDF$Species <- gsub(".rdata","",modelDF$File)
# annualDF <- getBUGSfits(modelDF,param="psi.fs")
# annualDF$Year <- annualDF$ParamNu + 1979
# plotTS(annualDF)
# #look very nice!!!
# table(annualDF$Rhat<1.1)
# #seems to be fewer problems

#sparta models - run on HPC using own sparta jags file with naturraum as fixedeffect

source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6288453"
#mistake in data processing

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6323394"
#didnt run for all species

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_wo_eta/6329258"
#without eta

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_wo_rw/6354978"
#without eta and rw (missing Calopteryx splendens)
#trends between this and the last one are after similar

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6417285"
#fixed subsetting code, with eta, ecoregion 1 and ecoregion 2 - but only works for a subset

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6441706"
#fixed subsetting code, with eta, ecoregion 1 - but only works for a subset

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6466710"
#fixed subsetting code, with eta, ecoregion 1, simple initial values - works for all!!

mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6474342"
#fixed subsetting code, with eta, ecoregion 1, ecoregion 2,
#simple initial values - works for all!!

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_adult_")

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#46  2803

#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]

###Nation naturraum#############################

#included rw on persist and colonization
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_rw1/6336641"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]
#missing quite a few
#error in nodes for eta
#Failure to calculate log density

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_dynamic_nation_naturraum_adult_")

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#519  1368

#trends
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
#FALSE  TRUE 
#16    35

trendsDF$Rhat[trendsDF$Rhat>1.1]

###Double random walk####################################################

#include rw also on observation year effect
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_rw/6288917"
#mistake in data processing

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]
#none missing

allSpecies[!sapply(allSpecies,function(x)any(grepl(x,speciesFiles)))]
#none missing

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,myfile="out_dynamic_nation_naturraum_adult_")

#annual tims series
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS(annualDF)
table(annualDF$Rhat<1.1)
#FALSE  TRUE 
#832  2017

table(annualDF$Rhat<1.1,annualDF$Species)

#example a few
list.files(mdir)
mod <- readRDS(paste(mdir,speciesFiles[1],sep="/"))

###Posterior draws############################################

library(rjags)
mymodel <- readRDS(paste(mdir,"out_sparta_nation_naturraum_adult_Sympetrum vulgatum.rds",sep="/")) 
#mysamples <- jags.samples(mymodel, "psi.fs", n.iter=1000)

#https://mjskay.github.io/tidybayes/
#library(tidybayes)
#mymodel$sims.list
#mymodel$samples - mcmc.list

#each year is in a different column
#mymodel$sims.lis[,1]
hist(mymodel$sims.list$psi.fs[,1])
summary(mymodel$sims.list$psi.fs[,1])

#growth
growth <- mymodel$sims.list$psi.fs[,37]/mymodel$sims.list$psi.fs[,1]
hist(growth)
summary(growth)
quantile(growth,c(0.025,0.5,0.975))

mymodel <- readRDS(paste(mdir,myfile,sep="/")) 
getGrowth <- function(mymodel){
  mymodel$sims.list$psi.fs[,37]/mymodel$sims.list$psi.fs[,1]
  quantile(growth,c(0.025,0.5,0.975))
}

#read in each model and apply function
#year, species and 1000 draw in subsequent columns

#sample psi.fs

myfiles <- list.files(mdir)

randomMatrix <- ldply(myfiles,function(myfile){
  
mymodel <- readRDS(paste(mdir,myfile,sep="/"))
  
out <- ldply(1:37,function(i){
  sims <- mymodel$sims.list$psi.fs[sample(1:56250,1000),i]
  cbind(Year=i,sims,Run=1:1000)
})  

out$Species <- myfile

return(out)

})  

randomMatrix$Species <- gsub("out_sparta_nation_naturraum_adult_",
                             "",randomMatrix$Species) 
randomMatrix$Species <- gsub(".rds",
                             "",randomMatrix$Species)
randomMatrix <- dcast(randomMatrix,Species+Year~Run,value.var="sims")
save(randomMatrix,file="randomMatrix.RData")

###get z########################################################

#HPC_update file
#Calculate parameters for each site??? take too long
#see what richness looks like

plotTS <- function(x){
  require(ggplot2)
  g1 <- ggplot(x)+
    geom_line(aes(x=Year,y=mean))+
    geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
    facet_wrap(~Code,ncol=6,scales="free_y")+
    ylab("Predicted occupancy")
  print(g1)
}
ggsave("plots/ts_scaled.png",height=10,width=7)

