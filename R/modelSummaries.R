#summaryPlots
library(rgdal)
library(ggplot2)
library(plyr)
library(reshape2)
library(ggthemes)

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

###Species list#################################################################################################

mySpecies <- read.delim("model-auxfiles/speciesTaskID_adult.txt",as.is=T)$Species

allSpecies <- read.delim("specieslist_odonata.txt",as.is=T)$Species

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

### Dynamic nation naturraum#############################

#included rw on persist and colonization
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_rw1/6336641"

#updates model fixing the site subset error
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum/6513010"
#missing 8 species - errors in node eta (site-level random effect)
#still have the problems at the start of the time series - upward bias
#overly smoothed compared to simple model

#as above but without random priors and eta
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum/6524659"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

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

###FINAL: Sparta models##########################################################

#sparta models - run on HPC using own sparta jags file with naturraum as fixedeffect

source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')

# mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs#/Odonata_adult_nation_naturraum_sparta/6466710"
# #fixed subsetting code, with eta, ecoregion 1, simple initial values - works for all!!
# #fixed effects are dnorm(0, 0.001)
# #intercepts are dnorm(0, 0.001)
# 
# mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/6474342" 
#fixed subsetting code, with eta, ecoregion 1, ecoregion 2,
#simple initial values - works for all!!

#updated data with revised files for Hessen and BW
#mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/7469396"
#"Coenagrion armatum"      "Coenagrion ornatum"      "Cordulegaster bidentata" "Sympetrum flaveolum"     "Sympetrum fonscolombii"
#problemSpecies <- c("Coenagrion armatum","Coenagrion ornatum","Cordulegaster bidentata","Sympetrum flaveolum","Sympetrum fonscolombii")
#error message is:
#Non-finite boundary in truncated normal

#updated data - with unbounded priors
#mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/7479651"
#10 missing

#updated with more bounded priors, dnorm(0,0.25) and midraum #instead of raum
#mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/7495806"
#problemSpecies <- c("Erythromma lindenii","Ischnura elegans","Lestes sponsa","Orthetrum albistylum","Sympetrum striolatum")

#updated with uniform priors on sd
#mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/7502229"

#with wide normal priors dnorm(0,0.001), 30,000
#mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs#/Odonata_adult_nation_naturraum_sparta/7503029"
#preds of Boyeria irene are quite different to last time...

#updated another 20,000 iterations - used for revisions
mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta/7503029/7578391"

#again with wide normal priors dnorm(0,0.001) and 50,000 iter (on slurm)
#mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs#/Odonata_adult_nation_naturraum_sparta/slurm/=Odonata_adult_nation_naturraum/55297"
#"Erythromma lindenii" (33) "Ischnura elegans"  (39)   "Lestes sponsa"(43)        
#"Orthetrum albistylum" (56) "Sympetrum striolatum" (76)

#Error in checkForRemoteErrors(val) : 
#  3 nodes produced errors; first error: LOGIC ERROR:
#  Non-finite boundary in truncated normal
#Please send a bug report to martyn_plummer@users.sourceforge.net

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_adult_")

#annual time series in occupancy
annualDF <- getBUGSfits(modelDF,param="psi.fs")
annualDF$Year <- annualDF$ParamNu + 1979
plotTS(annualDF)
table(annualDF$Rhat<1.1)
summary(annualDF$Rhat[annualDF$Rhat>1.1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.104   1.127   1.163   1.164   1.196   1.261

annualDF$Species[annualDF$Rhat>1.1]
#FALSE  TRUE 
#46  2803 

#trends in occupancy
trendsDF <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDF$Rhat<1.1)
trendsDF$Rhat[trendsDF$Rhat>1.1]
trendsDF$Species[trendsDF$Rhat>1.1]

#annual detection probabilities
annualDF <- getBUGSfits(modelDF,param="annual.p")
annualDF$Year <- annualDF$ParamNu + 1979
table(annualDF$Rhat<1.1)
annualDF$Species[annualDF$Rhat>1.1]

#mean detection probability
detprobDF <- getBUGSfits(modelDF,param="mean.p")
table(detprobDF$Rhat<1.1)
detprobDF$Rhat[detprobDF$Rhat>1.1]
detprobDF$Species[detprobDF$Rhat>1.1]
summary(detprobDF$mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01379 0.19957 0.24247 0.25958 0.31237 0.48551

#bpv
bpvDF <- getBUGSfits(modelDF,param="bpv")
summary(bpvDF$mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1267  0.4317  0.4788  0.4643  0.5018  0.9322
sort(bpvDF$mean)
#which ones are at the extreme
subset(bpvDF,mean<0.15)#Orthetrum albistylum
subset(bpvDF,mean>0.85)#Oxygastra curtisii 

### plotting ##################################

#Increasing the size of the plots, even if it takes up more space, as they are difficult to read
#Add the  % change to the title for each species
#Colour and/or group/order the plots according to the time series cluster that they fall in

#occupancy model
write.csv(annualDF,file="derived-data/annualOccupancies.csv")
plotTS_scales(annualDF)

ggsave("plots/ts_scaled.png",height=10,width=7)

#add occupancy data as rug plots
speciesAnnualObs <- readRDS("speciesAnnualObs.rds")
speciesAnnualObs$Code <- annualDF$Code[match(speciesAnnualObs$Species,annualDF$Species)]
speciesAnnualObs$dummy <- 0

#first 40
species1 <- sort(unique(annualDF$Species))[1:40]
plotTSwithRugs(subset(annualDF,Species %in% species1),subset(speciesAnnualObs,Species %in% species1)) 
ggsave("plots/ts_rugged_group1_scaled.png",height=10,width=7)

#next 40
species1 <- sort(unique(annualDF$Species))[41:77]
plotTSwithRugs(subset(annualDF,Species %in% species1),subset(speciesAnnualObs,Species %in% species1)) 
ggsave("plots/ts_rugged_group2_scaled.png",height=10,width=7)

#annual detection model
plotDetections(annualDF)
ggsave("plots/ts_detection.png",height=9,width=7)

#mean detection probability plot
detprobDF <- arrange(detprobDF,desc(mean))
detprobDF$Species <- factor(detprobDF$Species,levels=detprobDF$Species)
summary(detprobDF$mean)
ggplot(detprobDF)+
  geom_bar(aes(x=Species,y=mean),stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,hjust=1))+
  ylab("Mean detection probability")+
  xlab("")+
  coord_flip()
ggsave("plots/mean_detection.png",height=10,width=5)
  

### plot time series separately #############################

#used for shiny app by Edwin

head(annualDF)

for(i in mySpecies){
ggplot(subset(annualDF, Species == i))+
  #geom_point(aes(x = Year, y = mean))+
  geom_line(aes(x = Year, y = mean))+
  geom_ribbon(aes(x = Year, ymin = X2.5., ymax = X97.5.), alpha = 0.7)+
  theme_few()+ylab("Occupancy proportion")+ylim(0,1)
  
  ggsave(file=paste0("plots/species/",i,"_ts.png"),width = 5,height = 4)
}


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
  mymodel$sims.list$psi.fs[1:4500,37]/mymodel$sims.list$psi.fs[,1]
  quantile(growth,c(0.025,0.5,0.975))
}

#read in each model and apply function
#year, species and 1000 draw in subsequent columns

#sample psi.fs

myfiles <- list.files(mdir)

randomMatrix <- ldply(myfiles,function(myfile){
  
mymodel <- readRDS(paste(mdir,myfile,sep="/"))
  
out <- ldply(1:37,function(i){
  sims <- mymodel$sims.list$psi.fs[sample(1:4500,1000),i]
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
#Get z for each site
mdir <- "C:/Users/db40fysa/Nextcloud/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_sparta_upated/6522666"

#check files
outA <- readRDS(paste(mdir,"out_sparta_updated_z_adult_Aeshna affinis.rds",sep="/"))
outB <- readRDS(paste(mdir,"out_sparta_updated_z_adult_Aeshna caerulea.rds",sep="/"))

#site dimension for each species
mydims <- as.numeric()
for(i in 1:length(list.files(mdir))){
  outB <- readRDS(paste(mdir,list.files(mdir)[i],sep="/"))
  mydims[i] <- dim(outB)[2]
  rm(outB)
}

#check with number of sites calculated direct from data
load("allSpeciesSites.RData")
temp <- ddply(allSpeciesSites,.(Species),summarise,
              nuSites = length(unique(MTB_Q)))
all(temp$nuSites==mydims)#TRUE!!!!

#get list of common MTBQs for all species
mtbqSummary <- ddply(allSpeciesSites,.(MTB_Q),
                     summarise,nuSpecies=length(unique(Species)))
mtbqSummary <- subset(mtbqSummary,nuSpecies==77)
nrow(mtbqSummary)
#3243
allSpeciesSites <- subset(allSpeciesSites,MTB_Q %in% mtbqSummary$MTB_Q)
allSpeciesSites <- arrange(allSpeciesSites,Species,siteIndex)


#read in file for first species
out1 <- readRDS(paste(mdir,"out_sparta_updated_z_adult_Aeshna affinis.rds",sep="/"))

#subset to common files
out1 <- out1[,allSpeciesSites$siteIndex[allSpeciesSites$Species=="Aeshna affinis"],]
dim(out1)

myspecies <- unique(allSpeciesSites$Species)

#loop through each subsequent file, add them together,
for(i in 2:50){
  
    #read in file  
    outB <- readRDS(paste(mdir,list.files(mdir)[i],sep="/"))
    
    #restrict to common sites
    outB <- outB[,allSpeciesSites$siteIndex[allSpeciesSites$Species==myspecies[i]],]
    
    #add files together
    out1 <- out1 + outB
   
    #clean
    rm(outB)
}
#done in the RStudio Server
#files too big for here

#read in output done on the server
outAll <- readRDS("z_ALL.rds")
dim(outAll)#sim,site,year
max(outAll)
mean(outAll)#21

#we want mean number of species at a site in each year

#get number of predicted species at each site in each year
richness <- apply(outAll,c(1,3),mean) #per sim/year
dim(richness)

#get mean and 95% across sims
meanRichness <- apply(richness,2,mean)
lowerQRichness <- apply(richness,2,function(x)quantile(x,0.025))
upperQRichness <- apply(richness,2,function(x)quantile(x,0.975))
Year <- 1980:2016
richnessDF <- data.frame(Year,meanRichness,lowerQRichness,upperQRichness)

ggplot(richnessDF)+
  geom_line(aes(x=Year,y=meanRichness))+
  geom_ribbon(aes(x=Year,ymin=lowerQRichness,ymax=upperQRichness))
#almost the same as before....
#save as the random matrix...or keep with original?


### for iDiv defense ####
library(ggthemes)

myspecies <- "Sympetrum danae"
myspecies <- "Crocothemis erythraea"
year <- 1990:2016

#year 1
ggplot(subset(annualDF,Species==myspecies & Year == 1990))+
  ylim(0,0.65) +
  scale_x_continuous(breaks = c(1990,2000,2010),labels= c(1990,2000,2010),limits=c(1990,2016))+
  geom_linerange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+
  theme_few()+
  xlab("Year")+ylab("Mean occupancy")
ggsave(paste0("gifs/time-series/CE/timeseries_year",1990,".png"),width=4,height=3)

#next years
for(i in 2:length(year)){
  
ggplot(subset(annualDF,Species==myspecies & Year >= 1990 & Year <= year[i]))+
  ylim(0,0.65) +
  scale_x_continuous(breaks = c(1990,2000,2010),labels= c(1990,2000,2010),limits=c(1990,2016))+
  geom_line(aes(x = Year, y = mean))+
  geom_ribbon(aes(x=Year,ymin=X2.5.,ymax=X97.5.),alpha=0.5)+
  theme_few()+
  xlab("Year")+ylab("Mean occupancy")

ggsave(paste0("gifs/time-series/CE/timeseries_year",year[i],".png"),width=4,height=3)
}
  
  
  

