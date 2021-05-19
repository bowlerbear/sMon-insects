library(tidyverse)
library(cowplot)

#### overarching questions ####

#patterns of change - multi-faceting nature
#which communities most likely to change?
#property of the community or the environment

### distribution change summary ####

#plot map mean vs variance
#median vs IQR
#do metrics covary
#plot map with 2d legend

### output #####
list.files("model-outputs") %>% str_subset("simpleGAMS")

output <- readRDS("model-outputs/simpleGAMS_1990_2017.rds")
#output <- readRDS("model-outputs/simpleGAMS_1990_2017_k30.rds")
output$Year <- output$yearIndex + 1989

speciesMTB <- output %>%
              group_by(Species,MTB) %>%
              summarise(meanPred = mean(preds)) %>%
              filter(meanPred>0.01)#remove sites where species is probably absent

#### mtb trends #####

#get trend per MTBQ
mtbqTrends <- output %>%
  group_by(MTB,lon,lat,Species) %>%
  do(model = broom::tidy(lm(preds ~ Year, data = .))) %>% 
  unnest(model) %>%
  filter(term=="Year") %>%
  rename(trend="estimate",trend_se="std.error")

saveRDS(mtbqTrends,file="splines/mtbqTrends.rds")

#### species summary ####

mtbqTrends <- readRDS("splines/mtbqTrends.rds")

#remove MTB where species is not present
mtbqTrends <- inner_join(mtbqTrends,speciesMTB)

summary(mtbqTrends$trend)
quantile(mtbqTrends$trend,c(0.25,0.5,0.75))

#total 5% for stable
threshold <- 0.05/length(unique(output$Year))

speciesSummary <- mtbqTrends %>%
  group_by(MTB,lon,lat) %>%
  summarize(meanTrends = mean(trend),
            sdTrends = sd(trend),
            nuIncreasing = length(trend[trend>threshold]),
            nuStable = length(trend[trend>(-1*threshold)&trend<threshold]),
            nuDecreasing = length(trend[trend<(-1*threshold)]),
            propDeclining = mean(trend<0),
            propIncreasing = mean(trend>0),
            medianTrends = quantile(trend,0.5),
            lowerQ = quantile(trend,0.25),
            upperQ = quantile(trend,0.75)) %>%
  mutate(IQR = upperQ - lowerQ) %>%
  ungroup()

pairs(speciesSummary[,4:13])


#plot as maps
summary(speciesSummary$IQR)
hist(speciesSummary$IQR)
g1 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=IQR))+
  scale_color_viridis_c()+theme_void()

summary(speciesSummary$sdTrends)
hist(speciesSummary$sdTrends)
summary(speciesSummary$sdTrends)
speciesSummary$sdTrends[speciesSummary$sdTrends>0.006] <- 0.006
g2 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=sdTrends))+
  scale_color_viridis_c()+theme_void()

summary(speciesSummary$meanTrends)
hist(speciesSummary$meanTrends)
g3 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=meanTrends))+
  scale_color_viridis_c()+theme_void()

summary(speciesSummary$propDeclining)
hist(speciesSummary$propDeclining)
g4 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=propDeclining))+
  scale_color_viridis_c()+theme_void()

summary(speciesSummary$medianTrends)
hist(speciesSummary$medianTrends)
g5 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=medianTrends))+
  scale_color_viridis_c()+theme_void()

cowplot::plot_grid(g1,g2,g3,g4,g5,nrow=3)
cowplot::plot_grid(g3,g2,nrow=2)
cowplot::plot_grid(g5,g1,nrow=2)

#relationship between centre and variability

qplot(sdTrends,meanTrends, data = speciesSummary)+
  theme_bw()+
  geom_hline(yintercept=0, linetype="dashed",color="red")


qplot(IQR,medianTrends, data = speciesSummary)+
  theme_bw()+
  geom_hline(yintercept=0, linetype="dashed",color="red")

#### 2D map ####

#plot as 2D colour map

#put data into 9 values
speciesSummary$medianTrendsF <- as.numeric(cut_number(speciesSummary$propIncreasing,3))
speciesSummary$IQRF <- as.numeric(cut_number(speciesSummary$IQR,3))
speciesSummary$value <- speciesSummary$medianTrendsF
speciesSummary$value[speciesSummary$IQRF==2] <- speciesSummary$value[speciesSummary$IQRF==2] + 3
speciesSummary$value[speciesSummary$IQRF==3] <- speciesSummary$value[speciesSummary$IQRF==3] + 6
table(speciesSummary$value)

speciesSummary$meanTrendsF <- as.numeric(cut_number(speciesSummary$meanTrends,3))
speciesSummary$sdTrendsF <- as.numeric(cut_number(speciesSummary$sdTrends,3))
speciesSummary$value <- speciesSummary$meanTrendsF
speciesSummary$value[speciesSummary$sdTrendsF==2] <- speciesSummary$value[speciesSummary$sdTrendsF==2] + 3
speciesSummary$value[speciesSummary$sdTrendsF==3] <- speciesSummary$value[speciesSummary$sdTrendsF==3] + 6
table(speciesSummary$value)

bvColors=c("#be64ac","#8c62aa","#3b4994",
           "#dfb0d6","#a5add3","#5698b9",
           "#e8e8e8","#ace4e4","#5ac8c8")
legendGoal=reshape2::melt(matrix(1:9,nrow=3))

test <- ggplot(speciesSummary, aes(x=lon,y=lat,colour = as.factor(value)))+ geom_point(shape=15)
test <- test + scale_colour_manual(name="",values=bvColors)
test <- test + guides(colour = guide_legend(nrow = 3,
                                            override.aes = list(size=10)))
test <- test + theme_void() + theme(legend.text=element_blank())
test <- ggdraw(test) + draw_text(text = "IQR of trends -->",x=0.84,y=0.67, size=8)
test <- ggdraw(test) + draw_text(text = "Median trends -->",x=0.67,y=0.47,angle=270, size=8)
test 

#remove space in the legend??
lg <- test #would not be true when making a map
lg <- lg + theme(axis.title.x=element_text(size=rel(1),color=bvColors[3])) + xlab("More Var2 -->")
lg <- lg + theme(axis.title.y=element_text(size=rel(1),color=bvColors[3])) + ylab("More Var1 -->")
lg <- lg+theme(axis.text=element_blank())
lg <- lg+theme(line=element_blank())
#put both plots on a grid
ggdraw()+ draw_plot(lg,0.1,0.7,width=0.2,height=0.2) +draw_plot(test,0.3,0,width=.7,height=.7)

#### gains/loses ##################################################

g1 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=nuIncreasing))+
  scale_color_viridis_c()+theme_void()+coord_equal()
g2 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=nuDecreasing))+
  scale_color_viridis_c()+theme_void()+coord_equal()
g3 <- ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=nuStable))+
  scale_color_viridis_c()+theme_void()+coord_equal()

cowplot::plot_grid(g1,g2,g3,nrow=1)

#as 2D plot
speciesSummary$nuIncreasingF <- as.numeric(cut_number(speciesSummary$nuIncreasing,3))
speciesSummary$nuDecreasingF <- as.numeric(cut_number(speciesSummary$nuDecreasing,3))
speciesSummary$value <- speciesSummary$nuIncreasingF
speciesSummary$value[speciesSummary$nuDecreasingF==2] <- speciesSummary$value[speciesSummary$nuDecreasingF==2] + 3
speciesSummary$value[speciesSummary$nuDecreasingF==3] <- speciesSummary$value[speciesSummary$nuDecreasingF==3] + 6
table(speciesSummary$value)

library(cowplot)
bvColors=c("#be64ac","#8c62aa","#3b4994",
           "#dfb0d6","#a5add3","#5698b9",
           "#e8e8e8","#ace4e4","#5ac8c8")
legendGoal=reshape2::melt(matrix(1:9,nrow=3))

test <- ggplot(speciesSummary, aes(x=lon,y=lat,colour = as.factor(value)))+ geom_point(shape=15)
test <- test + scale_colour_manual(name="",values=bvColors)
test <- test + guides(colour = guide_legend(nrow = 3,
                                            override.aes = list(size=10)))
test <- test + theme(legend.text=element_blank())
test <- ggdraw(test) + draw_text(text = "Nu decreasing -->",x=0.84,y=0.71, size=8)
test <- ggdraw(test) + draw_text(text = "Nu increasing -->",x=0.67,y=0.53,angle=270, size=8)
test 

### richness change #####

#sum occu preds for each year and MTB (sum across species)
outputCommunity <- output %>%
                    group_by(MTB,lon,lat,Year) %>%
                    summarise(commOcc=sum(preds))
summary(outputCommunity$commOcc)

occChange <- outputCommunity %>%
  group_by(MTB,lon,lat) %>%
  do(model = broom::tidy(lm(commOcc ~ Year, data = .))) %>% 
  unnest(model) %>%
  filter(term=="Year") %>%
  rename(trend="estimate",trend_se="std.error")

g1 <- ggplot(occChange)+ geom_point(aes(x=lon,y=lat,colour=trend))+
  scale_color_viridis_c("Sp richness trend")+theme_void()

### composition change ####

library(vegan)
#years in rows and species in columns

#test with one MTB
output1 <- subset(output,MTB==916)
#organise data into a community matrix
commMatrix <- reshape2::acast(output1,Year~Species,value.var="preds")
#get the distance values
distanceMatrix <- vegdist(commMatrix,method="bray")
#extract comparison to the first year
indices <- as.matrix(distanceMatrix)[,1]
#put into a data frame
distanceDF <- data.frame(Year = as.numeric(as.factor(names(indices))), 
                         values = as.numeric(indices))
#do a linear regression
summary(lm(values ~ Year, data = distanceDF))$coefficients[2,1]

#package into a function

getDissimilarityTrend <- function(output1){
  commMatrix <- reshape2::acast(output1,Year~Species,value.var="preds")
  distanceMatrix <- vegan::vegdist(commMatrix,method="bray")
  indices <- as.matrix(distanceMatrix)[,1]
  distanceDF <- data.frame(Year = as.numeric(as.factor(names(indices))), 
                           values = as.numeric(indices))
  trend = summary(lm(values ~ Year, data = distanceDF))$coefficients[2,1]
}

outputDissimilarity <- output %>%
  group_by(MTB,lon,lat) %>%
  do(broom::tidy(getDissimilarityTrend(.)))

#some outliers
threshold <- quantile(outputDissimilarity$x,0.99)
outputDissimilarity$x[outputDissimilarity$x>threshold] <- threshold

g2 <- ggplot(outputDissimilarity)+ geom_point(aes(x=lon,y=lat,colour=x))+
  scale_color_viridis_c("Dissimilarity trend")+theme_void()

cowplot::plot_grid(g1,g2,nrow=2)


### 2D map ####

occChange$Dissimilarity <- outputDissimilarity$x
occChange$meanTrendsF <- as.numeric(cut_number(occChange$trend,3))
occChange$sdTrendsF <- as.numeric(cut_number(occChange$Dissimilarity,3))
occChange$value <- occChange$meanTrendsF
occChange$value[occChange$sdTrendsF==2] <- occChange$value[occChange$sdTrendsF==2] + 3
occChange$value[occChange$sdTrendsF==3] <- occChange$value[occChange$sdTrendsF==3] + 6
table(occChange$value)

library(cowplot)
bvColors=c("#be64ac","#8c62aa","#3b4994",
           "#dfb0d6","#a5add3","#5698b9",
           "#e8e8e8","#ace4e4","#5ac8c8")
legendGoal=reshape2::melt(matrix(1:9,nrow=3))

test <- ggplot(occChange, aes(x=lon,y=lat,colour = as.factor(value)))+ geom_point(shape=15)
test <- test + scale_colour_manual(name="",values=bvColors)
test <- test + guides(colour = guide_legend(nrow = 3,
                                            override.aes = list(size=10)))
test <- test + theme_void() + theme(legend.text=element_blank())
test <- ggdraw(test) + draw_text(text = "Incr dissimilarity -->",x=0.84,y=0.67, size=8)
test <- ggdraw(test) + draw_text(text = "Incr richness -->",x=0.67,y=0.47,angle=270, size=8)
test 


qplot(Dissimilarity, trend, data = occChange)+
  theme_bw()+
  geom_hline(yintercept=0, linetype="dashed",color="red")+
  ylab("species richness trend")+xlab("disimilarity trend")

#### range shift type ####

##range filler
#or
#range expander
#change in number of occupied sites
#versus change in areas of convex hull

rangeFun <- function(myspecies){

speciesData <- subset(output,Species==myspecies)

#sum all predictions
speciesData$predsArea <- speciesData$preds * 31000
totalArea <- as.numeric(tapply(speciesData$predsArea,speciesData$Year,sum))
#change to m2

#get area 
library(sp)
library(rgdal)

#put hull around data for each year
dat <- subset(speciesData, preds>0.01)[,c("x","y","Year","preds")]

myYears <- sort(unique(speciesData$Year))
rangeHull <- sapply(myYears,function(x){
ydat <- dat[dat$Year==x,c("x","y")]
if(nrow(ydat)>10){
ch <- chull(ydat)
coords <- ydat[c(ch, ch[1]), ] 
#plot(ydat, pch=19)
#lines(coords, col="red")
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)),
                           proj4string=CRS("+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
rangeSize <- raster::area(sp_poly)
rangeSize <- rangeSize/31000#area in m2 of MTBQ 
return(rangeSize)
}else{
  NA
}
})

return(data.frame(Species=myspecies,Year=myYears, totalArea, rangeHull))

}

temp <- rangeFun("Aeshna affinis")
qplot(Year, totalArea, data=temp)
qplot(Year, rangeHull, data=temp)

#apply the function to all species
speciesRanges <- plyr::ldply(sort(unique(output$Species)),rangeFun)

#how many missing values
speciesRanges %>%
  group_by(Year) %>%
  summarise(naHull = sum(is.na(rangeHull))) %>%
  print(n=40)

#exclude all those rare???
rareSpecies <- unique(subset(speciesRanges,is.na(rangeHull))$Species)

#now get the trends in each
areaChanges <- speciesRanges %>%
  filter(!Species %in% rareSpecies) %>%
  group_by(Species) %>%
  do(model = broom::tidy(lm(totalArea ~ Year, data = .))) %>% 
  unnest(model) %>%
  filter(term=="Year")%>%
  rename(area="estimate")

hullChanges <- speciesRanges %>%
  filter(!Species %in% rareSpecies) %>%
  filter(!is.na(rangeHull)) %>%
  group_by(Species) %>%
  do(model = broom::tidy(lm(rangeHull ~ Year, data = .))) %>% 
  unnest(model) %>%
  filter(term=="Year") %>%
  rename(hull="estimate")

#plotting together
all(areaChanges$Species==hullChanges$Species)
areaChanges$hull <- hullChanges$hull
ggplot(areaChanges,aes(x=area,y=hull))+
  geom_point()+
  #geom_text(aes(label=Species))+
  geom_hline(yintercept=0) + geom_vline(xintercept = 0)+
  xlab("trend in range area")+ylab("trend in range extent")+
  theme_bw()+
  geom_smooth(data=subset(areaChanges,area>0),method="lm")+
  geom_smooth(data=subset(areaChanges,area<0),method="lm")

#increasing species - are expanding their hull
#decreasing species are not changing in their hull

#### strongholds ##############################################

#change in species strongholds versus regions of rarity

#get region where species is common and where it is rare
#and get trend in each of those regions
output$MTB_CoarseNatur <- gsub("NaturrÃ¤ume Deutschland/","",
                               as.character(output$MTB_CoarseNatur))

#look at change within these regions, 1996 onwards
speciesNaturChange <- output %>%
  filter(Year>1995 & MTB_CoarseNatur!="Alpen") %>%
  group_by(Species,MTB_CoarseNatur,Year) %>%
  summarise(preds = sum(preds)) %>%
  group_by(Species,MTB_CoarseNatur) %>%
  do(model = broom::tidy(lm(preds ~ Year, data = .))) %>% 
  unnest(model) %>%
  filter(term=="Year")%>%
  rename(perChange="estimate")
hist(speciesNaturChange$perChange)

#add on overall trends
speciesNaturChange$overallTrend <- areaChanges$area[match(speciesNaturChange$Species,areaChanges$Species)]
speciesNaturChange$overallTrend <- speciesNaturChange$overallTrend/10000

#initial stronghold
speciesStrong <- output %>%
                  filter(Year<=1995 & MTB_CoarseNatur!="Alpen") %>%
                  group_by(Species,MTB_CoarseNatur) %>%
                  summarise(preds = mean(preds))%>%
                  group_by(Species)%>%
                  filter(preds == max(preds))
table(speciesStrong$MTB_CoarseNatur)
                  
#merge
speciesNaturChange_Strong <- inner_join(speciesStrong,speciesNaturChange)
speciesNaturChange_Strong$strongholdTrend <- speciesNaturChange_Strong$perChange
g1 <- qplot(overallTrend,strongholdTrend,data=speciesNaturChange_Strong) +
      geom_smooth(data=subset(speciesNaturChange_Strong,overallTrend<0),method="lm")+
      geom_smooth(data=subset(speciesNaturChange_Strong,overallTrend>0),method="lm")+
      theme_bw()

#initial weakhold
speciesWeak <- output %>%
  filter(Year<=1995 & MTB_CoarseNatur!="Alpen") %>%
  group_by(Species,MTB_CoarseNatur) %>%
  summarise(preds = mean(preds))%>%
  group_by(Species)%>%
  filter(preds == min(preds))
table(speciesWeak$MTB_CoarseNatur)

#merge
speciesNaturChange_Weak <- inner_join(speciesWeak,speciesNaturChange)
speciesNaturChange_Weak$weakholdTrend <- speciesNaturChange_Weak$perChange
g2 <- qplot(overallTrend,weakholdTrend,data=speciesNaturChange_Weak)+
  geom_smooth(data=subset(speciesNaturChange_Weak,overallTrend<0),method="lm")+
  geom_smooth(data=subset(speciesNaturChange_Weak,overallTrend>0),method="lm")+
  theme_bw()

cowplot::plot_grid(g1,g2,nrow=1)

#### 1000 cuts ###############################################

##distribution of trends across MTBs
mtbqTrends <- inner_join(mtbqTrends,areaChanges[,c("Species","area")])
mtbqTrends$Direction <- ifelse(mtbqTrends$area>0,"Winning","Losing")

mtbqTrendsSummary <- mtbqTrends %>%
                      group_by(MTB,Direction) %>%
                      summarise(lowerQ = quantile(trend,0.25),
                                upperQ = quantile(trend,0.75)) %>%
                      mutate(IQR = upperQ - lowerQ)

ggplot(mtbqTrendsSummary)+
  geom_boxplot(aes(x=Direction, y=IQR))+
  theme_bw()+
  ylab("Trend variability")

##distribution of trends across species

speciesTrendsSummary <- mtbqTrends %>%
  group_by(Species,Direction) %>%
  summarise(area = mean(area),
            median = quantile(trend,0.5),
            lowerQ = quantile(trend,0.25),
            upperQ = quantile(trend,0.75)) %>%
  mutate(IQR = upperQ - lowerQ)


speciesTrendsSummary$area <- abs(speciesTrendsSummary$area)/10000
ggplot(speciesTrendsSummary,aes(x=area, y=IQR))+
  geom_point(aes(colour=Direction))+
  theme_bw()+
  ylab("Spatial trend variability")+
  geom_smooth(aes(colour=Direction),method="lm")+
  xlab("|nationwide area change|")

#### novel assemblages ####

#Where are novel assemblages occurring?
#use analogous of functional trait uniqueness??
#how different a community is to any other community

#round pred occupancies per decade
floor_decade    = function(value){ return(value - value %% 10) }
output$Decade <- floor_decade(output$Year)

outputDecade <- output %>%
                group_by(MTB,lon,lat,Decade,Species) %>%
                summarise(preds=mean(preds))

library(vegan)
#organise data into a community matrix
commMatrix <- reshape2::acast(outputDecade,MTB+Decade~Species,value.var="preds")
#get the distance values
distanceMatrix <- as.matrix(vegdist(commMatrix,method="bray"))
distanceMatrix[1:5,1:5]
#select columns with 1990
distanceMatrix <- distanceMatrix[,grepl("_1990",colnames(distanceMatrix))]
#select rows with 2010
distanceMatrix <- distanceMatrix[grepl("_2010",rownames(distanceMatrix)),]

#get maximum dissimilarity for each MTB
distanceMatrix <- apply(distanceMatrix,1,max)

#put into a data frame
distanceDF <- data.frame(Year = names(distanceMatrix), 
                         values = as.numeric(distanceMatrix))
distanceDF$MTB <- sapply(distanceDF$Year,function(x)strsplit(x,"_")[[1]][1])

speciesSummary$Uniqueness <- distanceDF$values[match(speciesSummary$MTB,distanceDF$MTB)]
ggplot(speciesSummary)+ geom_point(aes(x=lon,y=lat,colour=Uniqueness))+
  scale_color_viridis_c("Dissimilarity")+theme_void()

#### changepoints ####


#get derivatives of temporal gams..


#### betadiversity ####

library(betapart)
#matrix (x) codifying the presence (1) or absence (0) of m species (columns) in n sites (rows)
#betapart.core(x)
#?beta.multi - This function computes the total dissimilarity across all n sites, #along with the turnover and nestedness components of that dissimilarity. The input #x may be a presence–absence matrix or a betapart object.
#beta.pair - pairwise between‐site values
#beta.temp(x, y, index.family)
# multiple‐site and pairwise partitions of beta diversity.

#where a is the number of shared species between two cells, b the number of species unique to the poorest site and c the number of species unique to the richest site. 

# Beta diversity, that is, the variation in species composition among sites, can be the result of species replacement between sites (turnover) and species loss from site to site (nestedness).


data(bbsData)
data(bbsData)
bbs.t <- beta.temp(bbs1980, bbs2000, index.family ='so')

# plotting root transformed components
with(bbs.t, plot(sqrt(beta.sim) ~ sqrt(beta.sne), type='n', ylab=expression(sqrt(beta[sim])), xlab=expression(sqrt(beta[sne]))))
with(bbs.t, text(y= sqrt(beta.sim), x=sqrt(beta.sne), labels=rownames(bbs1980)))
#sor = dissimilarity
#sim = turnover
#sne = nestedness
str(bbs.t)
plot(hclust(pair.s$beta.sne, method='average'), hang=-1, main='', sub='', xlab='')

#### GDM ####

#Get mean temperature for each grid - Got!!!
#https://rfunctions.blogspot.com/2015/08/calculating-beta-diversity-on-grid.html
#GDM - beta diversity - as a function of lat and lon
#https://cran.r-project.org/web/packages/gdm/vignettes/gdmVignette.pdf


#### LDA modelling ####


#topic modelling (LDA) to identify groups of species that co-occur in each decade?
use 
####rare species ####

#what are rare species doing versus common species?
#classify species into rare, medium and common

#based on data in first 5 years
speciesMeans <- output %>%
                filter(Year<1996) %>%
                group_by(Species,Year) %>%
                summarise(sumPred = sum(preds))%>%
                group_by(Species)%>%
                summarise(meanPred = median(sumPred))
speciesMeans$Commoness <- cut_number(speciesMeans$meanPred,4)
table(speciesMeans$Commoness)
levels(speciesMeans$Commoness) <- c("rare","intermediate rare",
                                    "intermediate common","common")
output$Commoness <- speciesMeans$Commoness[match(output$Species,
                                                     speciesMeans$Species)]

#for each species get occupancy Change between first and last year
occChange <- output %>%
  group_by(MTB,x,y,lat,lon,Species,Commoness) %>%
  summarise(occChange = log(preds[Year==2016]/preds[Year==1990]))%>%
  group_by(MTB,x,y,lat,lon,Commoness) %>%
  summarise(occChange = median(occChange))

threshold <- quantile(occChange$occChange,0.95)
occChange$occChange[occChange$occChange>threshold] <- threshold 
  
q1 <- ggplot(subset(occChange,Commoness=="rare"))+ geom_point(aes(x=lon,y=lat,colour=occChange))+scale_color_viridis_c()+theme_void()+coord_equal()
q2 <- ggplot(subset(occChange,Commoness=="intermediate rare"))+ geom_point(aes(x=lon,y=lat,colour=occChange))+scale_color_viridis_c()+theme_void()+coord_equal()
q3 <- ggplot(subset(occChange,Commoness=="intermediate common"))+ geom_point(aes(x=lon,y=lat,colour=occChange))+scale_color_viridis_c()+theme_void()+coord_equal()
q4 <- ggplot(subset(occChange,Commoness=="common"))+ geom_point(aes(x=lon,y=lat,colour=occChange))+scale_color_viridis_c()+theme_void()+coord_equal()

cowplot::plot_grid(q1,q2,q3,q4)

#### climate data ####

#Question - what is the best predictor mean temperature or change in temperature?

climdata <- readRDS("C:/Users/db40fysa/Nextcloud/sMon/Gis-DataData/Klima Deutschland/DWD/Climate_by_MTBQ/summary/allclimvarsSummary_MTBQ.rds")
#collapse to MTB leve;s
climdata_MTB <- climdata %>%
                group_by(Value) %>%
                summarise(springtempTrend = mean(springtempTrend),
                          springprecipTrend = mean(springprecipTrend)) %>%
                ungroup()%>%
                rename(MTB = "Value")
          
        

#PCA
mydata <- climdata[,grepl("Trend",names(climdata))]
mydata <- mydata[,!grepl("Trend_se",names(mydata))]
fit <- princomp(mydata, cor=TRUE)
summary(fit) 
loadings(fit) 
biplot(fit) 

#mean trend
library(mgcv)
speciesClims <- inner_join(speciesSummary,climdata_MTB)

gam1 <- gam(meanTrends ~ springprecipTrend + springtempTrend +
              lat + 
              s(lon,lat),
            data=speciesClims)
summary(gam1)

qplot(springprecipTrend,meanTrends,data=speciesClims)
qplot(springtempTrend,meanTrends,data=speciesClims)

#richness trends
speciesClims <- inner_join(occChange,climdata_MTB)

gam1 <- gam(trend ~ springprecipTrend + springtempTrend +
              lat + 
              s(lon,lat),
            data=speciesClims)
summary(gam1)

qplot(springprecipTrend,trend,data=speciesClims)
qplot(springtempTrend,trend,data=speciesClims)

#dissimilarity trend
names(outputDissimilarity)[4] <- "trend"
speciesClims <- inner_join(outputDissimilarity,climdata_MTB)

gam1 <- gam(trend ~ springprecipTrend + springtempTrend +
              lat + 
              s(lon,lat),
            data=speciesClims)
summary(gam1)

qplot(springprecipTrend,trend,data=speciesClims)
qplot(springtempTrend,trend,data=speciesClims)

#https://paul-buerkner.github.io/brms/reference/car.html
#https://stats.stackexchange.com/questions/277/spatial-statistics-models-car-vs-sar

#### end ####