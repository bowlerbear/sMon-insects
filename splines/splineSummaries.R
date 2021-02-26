library(zoo)
library(plyr)

specieslist = read.delim("model-auxfiles/speciesTaskID_adult.txt")

#### simple GAMS #####

# see analysis HPC nation spline MTB
output <- readRDS("model-outputs/simpleGAMS_1990_2017.rds")
output$Year <- output$yearIndex + 1989

#average over 3 year periods
outputS <- ddply(output,.(Species,x,y),function(x){
  x$rollPred = rollmedian(x$preds, 3, align="left",na.pad = TRUE) 
  return(x)
})

#subset to selection of years
outputSS <- subset(outputS,Year %in% c(1990,2002,2015))

for(i in 1:length(specieslist$Species)){
  
  #plot for each species
  qplot(x,y,data=subset(outputSS,Species==specieslist$Species[i]),colour=rollPred)+
    facet_wrap(~Year,ncol=4)+
    scale_colour_viridis_c("Occupancy", option = "magma", direction = -1)+
    theme_void()
  
  ggsave(paste0("splines/simpleGAMS/GAM_",specieslist$Species[i],".png"),width=8,height=2.5)
  
}


#### Trends ######

#get from modelSummaries file
trendEstimates <- trendsDF
trendEstimates$Trend <- "insignificant"
trendEstimates$Trend[trendEstimates$X2.5.>0 & trendEstimates$X97.5.>0]<-"significant increase"
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.<0]<-"significant decrease"
table(trendEstimates$Trend)


#### GAMS by trend classification ####

output$TrendClass <- trendEstimates$Trend[match(output$Species,trendEstimates$Species)]
output$Trend <- trendEstimates$mean[match(output$Species,trendEstimates$Species)]

#mean per yearly occupancy per trend class
outputClass <- ddply(output,.(TrendClass,Year,x,y),summarise,preds = mean(preds))

#rolling medians
outputClass_S <- ddply(outputClass,.(TrendClass,x,y),function(x){
  x$rollPred = rollmedian(x$preds, 3, align="left",na.pad = TRUE) 
  return(x)
})

#plot
qplot(x,y,
      data=subset(outputClass_S,Year %in% c(1990,2002,2015)),
      colour=rollPred)+
  facet_grid(TrendClass~Year)+
  scale_colour_viridis_c("Occupancy", option = "magma", direction = -1)+
  theme_void()

#### regions of median change ####

outputChange <- ddply(outputS,.(Species,x,y),summarise,Change=preds[Year==2010]/preds[Year==1990])
outputChange$TrendClass <- trendEstimates$Trend[match(outputChange$Species,trendEstimates$Species)]

outputChange_Median <- ddply(outputChange,.(x,y),summarise,Change=median(Change))
qplot(x,y, data=outputChange_Median,colour=Change)+
  scale_colour_viridis_c("Occupancy", option = "magma", direction = -1)+
  theme_void()

#increasing and decreasing species
outputChange_Class <- ddply(outputChange,.(TrendClass,x,y),summarise,Change=median(Change))

qplot(x,y, data=subset(outputChange_Class,TrendClass=="significant increase"),colour=Change)+
  scale_colour_viridis_c("Occupancy", option = "magma", direction = -1)+
  theme_void()
#increasing species are increasing the most in the north

qplot(x,y, data=subset(outputChange_Class,TrendClass=="significant decrease"),colour=Change)+
  scale_colour_viridis_c("Occupancy", option = "magma", direction = -1)+
  theme_void()
#decreasing species are decreasing most in the south

#### regions of compositional change ####

#proportion of shared species

#### novel communities ##############

### combination of species not seen previously


#### end ####

