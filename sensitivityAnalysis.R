#sensitivityAnalysis

#compare models with and without the phenology by year detection

###  standard modelDF ##################################

#run script on modelSummaries file 
#to get annualDF and trendsDF

### phenology change ###################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_phenologyChange/7464688" 
#two species are missing: "Coenagrion ornatum"      "Cordulegaster bidentata"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_phenologyChange_adult_")

#annual tims series
annualDFpc <- getBUGSfits(modelDF,param="psi.fs")
annualDFpc$Year <- annualDFpc$ParamNu + 1979
plotTS(annualDFpc)
table(annualDFpc$Rhat<1.1)
#FALSE  TRUE 
#78  2697

#trends
trendsDFpc <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDFpc$Rhat<1.1)
trendsDFpc$Rhat[trendsDFpc$Rhat>1.1]
#FALSE  TRUE 
#6    69 

### comparison #########################################

trendsDFcomparison <- merge(trendsDF,trendsDFpc,by="Species")
qplot(mean.x,mean.y,data=trendsDFcomparison)
cor.test(trendsDFcomparison$mean.x,trendsDFcomparison$mean.y)
cor.test(trendsDFcomparison$sd.x,trendsDFcomparison$sd.y)

annualDFcomparison <- merge(annualDF,annualDFpc,by=c("Species","Year"))
qplot(mean.x,mean.y,data=annualDFcomparison)+
  facet_wrap(~Species)
                                                          
temp <- ddply(annualDFcomparison,.(Species),function(x){
  cor(x$mean.x,x$mean.y)
})
median(temp$V1)


### mtbq subset ########################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_sparta_MTBQ_eachDecade/7503038"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]
#"Crocothemis erythraea" "Oxygastra curtisii"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_MTBQ_eachDecadeadult_")

#annual time series in occupancy
annualDFs <- getBUGSfits(modelDF,param="psi.fs")
annualDFs$Year <- annualDFs$ParamNu + 1979
plotTS(annualDFs)
table(annualDFs$Rhat<1.1)
summary(annualDFs$Rhat[annualDF$Rhat>1.1])

#trends in occupancy
trendsDFs <- getBUGSfits(modelDF,param="regres.psi")

### comparison #####################################

trendsDFcomparison <- merge(trendsDF,trendsDFs,by="Species")
qplot(mean.x,mean.y,data=trendsDFcomparison)
cor.test(trendsDFcomparison$mean.x,trendsDFcomparison$mean.y)
cor.test(trendsDFcomparison$sd.x,trendsDFcomparison$sd.y)

qplot(sd.x,sd.y,data=trendsDFcomparison)
qplot(mean.x,mean.y,data=trendsDFcomparison)


### dynamic model ##################################

mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_dynamic/7531072"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_dynamic_nation_naturraum_adult_")

#annual time series in occupancy
annualDFs <- getBUGSfits(modelDF,param="psi.fs")
annualDFs$Year <- annualDFs$ParamNu + 1979
plotTS(annualDFs)
table(annualDFs$Rhat<1.1)
summary(annualDFs$Rhat[annualDFs$Rhat>1.1])

#trends in occupancy
trendsDFs <- getBUGSfits(modelDF,param="regres.psi")

### comparison ####################################

trendsDFcomparison <- merge(trendsDF,trendsDFs,by="Species")
qplot(mean.x,mean.y,data=trendsDFcomparison)
cor.test(trendsDFcomparison$mean.x,trendsDFcomparison$mean.y)
cor.test(trendsDFcomparison$sd.x,trendsDFcomparison$sd.y)

qplot(sd.x,sd.y,data=trendsDFcomparison)
qplot(mean.x,mean.y,data=trendsDFcomparison)


### detection model #################################

#including coarse and mid naturruam on detection probability
mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_detModel/7564941"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]
#"Erythromma lindenii"  "Ischnura elegans"     "Lestes sponsa"        "Orthetrum albistylum" "Sympecma paedisca" [6] "Sympetrum striolatum"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_detModel_adult_")

#annual time series in occupancy
annualDFs <- getBUGSfits(modelDF,param="psi.fs")
annualDFs$Year <- annualDFs$ParamNu + 1979
plotTS(annualDFs)
table(annualDFs$Rhat<1.1)
summary(annualDFs$Rhat[annualDFs$Rhat>1.1])

#trends in occupancy
trendsDFs <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDFs$Rhat<1.1)
#FALSE  TRUE 
#9    62

#annual detection probabilities
annualDF <- getBUGSfits(modelDF,param="annual.p")
annualDF$Year <- annualDF$ParamNu + 1979
plotDetections(annualDF)
table(annualDF$Rhat<1.1)
annualDF$Species[annualDF$Rhat>1.1]
#worse convergence

#bpv
bpvDF <- getBUGSfits(modelDF,param="bpv")
summary(bpvDF$mean)
#better...

### comparison ####################################

trendsDFcomparison <- merge(trendsDF,trendsDFs,by="Species")
qplot(mean.x,mean.y,data=trendsDFcomparison)
cor.test(trendsDFcomparison$mean.x,trendsDFcomparison$mean.y)
cor.test(trendsDFcomparison$sd.x,trendsDFcomparison$sd.y)

qplot(sd.x,sd.y,data=trendsDFcomparison)
qplot(mean.x,mean.y,data=trendsDFcomparison)

### excluding naturalist data #######################

#including coarse and mid naturruam on detection probability
mdir <- "C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/model-outputs/Odonata_adult_nation_naturraum_statesOnly/7563237"

#do we have the models for all species?
speciesFiles <- list.files(mdir)
mySpecies[!sapply(mySpecies,function(x)any(grepl(x,speciesFiles)))]
#"Nehalennia speciosa"    "Somatochlora alpestris" "Sympetrum meridionale"  "Sympetrum sanguineum"   "Sympetrum striolatum"

#read in model summaries
modelDF <- getModelSummaries(mdir)
modelDF <- getCodeFromFile(modelDF,
                           myfile="out_sparta_nation_naturraum_statesOnly_adult_")

#annual time series in occupancy
annualDFs <- getBUGSfits(modelDF,param="psi.fs")
annualDFs$Year <- annualDFs$ParamNu + 1979
plotTS(annualDFs)
table(annualDFs$Rhat<1.1)
summary(annualDFs$Rhat[annualDFs$Rhat>1.1])

#trends in occupancy
trendsDFs <- getBUGSfits(modelDF,param="regres.psi")
table(trendsDFs$Rhat<1.1)
#FALSE  TRUE 
#3    69

#annual detection probabilities
annualDF <- getBUGSfits(modelDF,param="annual.p")
annualDF$Year <- annualDF$ParamNu + 1979
plotDetections(annualDF)
table(annualDF$Rhat<1.1)
annualDF$Species[annualDF$Rhat>1.1]
#worse convergence

#bpv
bpvDF <- getBUGSfits(modelDF,param="bpv")
summary(bpvDF$mean)

### comparison ####################################

trendsDFcomparison <- merge(trendsDF,trendsDFs,by="Species")
cor.test(trendsDFcomparison$mean.x,trendsDFcomparison$mean.y)
cor.test(trendsDFcomparison$sd.x,trendsDFcomparison$sd.y)

qplot(sd.x,sd.y,data=trendsDFcomparison)
qplot(mean.x,mean.y,data=trendsDFcomparison)

#almost the same

### end ###########################################