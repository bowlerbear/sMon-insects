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
