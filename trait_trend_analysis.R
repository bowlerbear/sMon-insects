source('C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects/R/sparta_wrapper_functions.R')
library(plyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(wesanderson)
library(MuMIn)

###get data###################################################################

setwd("C:/Users/db40fysa/Nextcloud/sMon/sMon-Analyses/Odonata_Git/sMon-insects")

#get traits data
load("alltraits.RData")

#get trends data
trendEstimates <- trendsDF

###missing data##############################################

#few missing flight start and end dates
alltraits$Species[is.na(alltraits$Flight_start)]
#"Aeshna viridis", "Boyeria irene", "Epitheca bimaculata"

#available for other countries
#Aeshna viridis  - 7.3/8.3
#Boyeria irene - 7.1/8.1
#Epitheca bimaculata - 5.1/6.2
alltraits$Flight_start[alltraits$Species=="Aeshna viridis"] <- 7.3
alltraits$Flight_end[alltraits$Species=="Aeshna viridis"] <- 8.3
alltraits$Flight_start[alltraits$Species=="Boyeria irene"] <- 7.1
alltraits$Flight_end[alltraits$Species=="Boyeria irene"] <- 8.1
alltraits$Flight_start[alltraits$Species=="Epitheca bimaculata"] <- 5.1
alltraits$Flight_end[alltraits$Species=="Epitheca bimaculata"] <- 6.2

#missing data for develop time - not included anyhow
alltraits$Suborder[alltraits$Species=="Oxygastra curtisii"] <- "Anisoptera"
###merge###########################################################

names(trendEstimates)[which(names(trendEstimates)=="mean")]<-"trend"
trendEstimates$Species[!trendEstimates$Species %in% alltraits$Species]
trendEstimates <- merge(trendEstimates,alltraits,by="Species")
nrow(trendEstimates)#77
trendEstimates$Genus <- sapply(trendEstimates$Species,function(x){
  strsplit(x," ")[[1]][1]})

#trends calculated as per:
#https://github.com/BiologicalRecordsCentre/RangeChangeSims
#https://github.com/BiologicalRecordsCentre/RangeChangeSims/blob/master/Occupancy/Occupancy_model.r

###trend summary##############################################################

summary(trendEstimates$trend)
hist(trendEstimates$trend)#pretty normal

#null
trendEstimates$Trend <- "insignificant"

#trend is annual proportional change in number of sites
#mean 1% change over the 35 year - less than total change 5%
#0.05/36
#[1] 0.001388889

#standard deviation to large to detect a 5% change

#stable
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.>0
                     & abs(trendEstimates$trend) < 0.001388889]<-"stable"

#increasing
trendEstimates$Trend[trendEstimates$X2.5.>0 & trendEstimates$X97.5.>0]<-"significant increase"

#decreasing
trendEstimates$Trend[trendEstimates$X2.5.<0 & trendEstimates$X97.5.<0]<-"significant decrease"

#summary
table(trendEstimates$Trend)

#habitats of decreasing species - cols 50 -65
subset(trendEstimates,Trend=="significant decrease")[,50:65]
#all standing or slow-flowing species - pond species

#median increase
summary(trendEstimates$trend[trendEstimates$Trend=="significant increase"])

summary(trendEstimates$trend[trendEstimates$Trend=="significant decrease"])


###red list################################################################

table(trendEstimates$RedList)

#order the factor

trendEstimates$RedList <- factor(trendEstimates$RedList,
                                 levels=c("not listed","least concern","near threatened",
                                          "vulnerable","endangered","critically endangered"))


table(trendEstimates$RedList)

g1 <- ggplot(trendEstimates)+
  geom_boxplot(aes(x=RedList,y=trend))+
  geom_hline(yintercept=0,linetype="dashed")+
  coord_flip()+
  geom_text(x=1,y=0.01,label="1")+
  geom_text(x=2,y=0.01,label="43")+
  geom_text(x=3,y=0.01,label="6")+
  geom_text(x=4,y=0.01,label="12")+
  geom_text(x=5,y=0.01,label="6")+
  geom_text(x=6,y=0.01,label="9")+
  theme_classic()+
  ylab("Long-term trend")+xlab("Red list category")

#or add numbers of species to axes label??
trendEstimates$RedListL <- trendEstimates$RedList
levels(trendEstimates$RedListL) <- 
  c("not listed (1)", "least concern (43)",
    "near threatened (6)","vulnerable (12)",
    "endangered (6)","critically endangered (9)")

g1 <- ggplot(trendEstimates)+
  geom_boxplot(aes(x=RedListL,y=trend))+
  geom_hline(yintercept=0,linetype="dashed")+
  coord_flip()+
  theme_classic()+
  ylab("Long-term trend")+xlab("Red list category")

summary(lm(trend~RedList,data=subset(trendEstimates,RedList!="not listed")))
summary(lm(trend~RedList-1,data=subset(trendEstimates,RedList!="not listed")))

#no difference between near threatened and least concern
#but all threatened groups have significantly lower trends

###winner, losers####################################################################
trendEstimates$Species <- factor(trendEstimates$Species,levels=c(trendEstimates$Species[order(trendEstimates$trend)]))
trendEstimates$Trend[trendEstimates$Trend=="stable"] <- "insignificant"
trendEstimates$Trend <- factor(trendEstimates$Trend,levels=c("significant decrease","insignificant","significant increase")) 

g2 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=trend,fill=Trend),
           stat="identity",width=rel(1))+
  theme_classic()+  
  theme(axis.text.x = element_blank())+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  ylab("Long-term trend")+xlab("Species")


library(cowplot)
plot_grid(g2,g1,nrow=1)
ggsave("plots/Trend_summary.png",width=8,height=3)


#as interactive plot
#using ggiraph
library(ggiraph)
g2 <- ggplot(trendEstimates)+
  geom_bar_interactive(aes(x=Species,y=trend,fill=Trend,tooltip=Species,data_id=Species),
           stat="identity",width=rel(1))+
  theme_classic()+  
  theme(axis.text.x = element_blank())+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  ylab("Long-term trend")+xlab("Species")

x <- girafe(ggobj = g2,
       options = list(
         opts_hover(css = "fill:black;")))
if( interactive() ) print(x)
library(htmlwidgets)
library(plotly)
htmlwidgets::saveWidget(x, "Longterm_trends.html")

###prop change############################################################

totalChange <- ddply(annualDF,.(Species),summarise,
                     change=mean[Year==2016]/mean[Year==1980],
                     growth=((mean[Year==2016]/mean[Year==1980])-1)*100)

summary(totalChange$change)
#Crocothemis erythraea has the largest amount of change

totalChange <- arrange(totalChange,change)

trendEstimates$change <- totalChange$change[match(trendEstimates$Species,
                                                  totalChange$Species)]

trendEstimates$Species <- factor(trendEstimates$Species,levels=c(totalChange$Species[order(totalChange$change)]))


g3 <- ggplot(trendEstimates)+
  geom_bar(aes(x=Species,y=change),
           stat="identity",width=rel(1))+
  theme_classic()+  
  scale_y_log10(breaks=c(0,0.01,0.1,1,10,100),
                labels=c(0,0.01,0.1,1,10,100))+
  theme(axis.text.x = element_blank())+
  ylab("Total occupancy change")+xlab("Species")


plot_grid(g2,g3,g1,ncol=1,labels=c("A","B","C"))
ggsave("plots/Trend_summary.png",width=6,height=8)

#relationship between total change and trend
ggplot(trendEstimates)+
  geom_point(aes(x=trend,y=change))+
  theme_bw()+
  scale_y_log10()+
  ylab("Total occupancy change")+
  xlab("Long-term trend")

cor.test(trendEstimates$trend,log10(trendEstimates$change))
#0.6742095

#change of increasing species
trendEstimates <- arrange(trendEstimates,trend)
incr <- subset(trendEstimates,Trend=="significant increase")
summary(incr$change)
tail(incr[,c("Species","change","trend")])
subset(annualDF,Species=="Crocothemis erythraea" & Year %in% c(1980,2016))

#change of decreasing species
decr <- subset(trendEstimates,Trend=="significant decrease")
summary(decr$change)
head(decr[,c("Species","change","trend")])
subset(annualDF,Species=="Sympetrum danae" & Year %in% c(1980,2016))

###choose traits####################################################################

#Flight_start
hist(trendEstimates$Flight_start)
#VoltinismProp
hist(trendEstimates$VoltinismProp)
#nuEuroGrids
hist(trendEstimates$nuEuroGrids)
#meanTemp
hist(trendEstimates$meanTemp)
#medHw/medAb
hist(trendEstimates$medHw/trendEstimates$medAb)
trendEstimates$DispPot <- trendEstimates$medHw/trendEstimates$medAb
#Overwintering
table(trendEstimates$Overwintering)
#Generalism
table(trendEstimates$Generalism)
#Flow
table(trendEstimates$Flow)
trendEstimates$Flow[which(trendEstimates$Flow=="standing (slow-flowing)")] <- "standing"
trendEstimates$Flow[which(trendEstimates$Flow=="standing, slow-flowing")] <- "standing, running"
table(trendEstimates$Flow)

#reorganise flight start
trendEstimates$cFlight_start <- as.character(trendEstimates$Flight_start)
trendEstimates$cFlight_start <- gsub("1","0",trendEstimates$cFlight_start)
trendEstimates$cFlight_start <- gsub("3","66",trendEstimates$cFlight_start)
trendEstimates$cFlight_start <- gsub("2","33",trendEstimates$cFlight_start)
trendEstimates$Flight_start <- as.numeric(trendEstimates$cFlight_start) 

table(trendEstimates$Generalism)
trendEstimates$Generalism <- ifelse(trendEstimates$Generalism=="generalist",1,0)

#Egg or not
trendEstimates$winterEgg <- ifelse(trendEstimates$Overwintering=="Egg",1,0)
table(trendEstimates$winterEgg)
table(trendEstimates$Overwintering)

#turn habitat variables into binary variables
myhabitats <- names(trendEstimates)[52:65]
trendEstimates[,myhabitats] <- apply(trendEstimates[,myhabitats],
                                     2,function(x)ifelse(x>1,1,0))

#add on running column
trendEstimates$running <- rowSums(trendEstimates[,c("river","stream")])
trendEstimates$running <- ifelse(trendEstimates$running>=1,1,0)

#define traits of interest
mytraits <- c("Flight_start","VoltinismProp","meanTemp",
              "DispPot","medHw")

allTraits <- c(mytraits,myhabitats,"running")

###save file#############################################

write.csv(trendEstimates,file="traits_trends_Odonata_Jan2020.csv",
          row.names=FALSE)

write.csv(annualDF,file="annualDF_Jan2020.csv",
          row.names=FALSE)

###scale traits##############################################

library(arm)
scaledVars <- trendEstimates[,allTraits]
scaledVars <- data.frame(lapply(scaledVars,arm::rescale))
names(scaledVars) <- sapply(names(scaledVars),function(x)
  paste0("s",x))
trendEstimates <- cbind(trendEstimates,scaledVars)

mytraits <- c("sFlight_start","sVoltinismProp","smeanTemp","smedHw")

###trait corr################################################

cor.test(trendEstimates$meanTemp,trendEstimates$bog)
#significant

library(corrplot)

traitNums <- Filter(is.numeric, trendEstimates[,allTraits])
overwinterDummy <- model.matrix(~Overwintering-1,data=trendEstimates)
traitsCor <- cor(cbind(traitNums,overwinterDummy))
corrplot(traitsCor)
traitsCorMelt <- melt(traitsCor)
temp <- subset(traitsCorMelt,abs(value)>0.4 & value!=1)
subset(temp,!duplicated(value))

###Species-level############################################
###habitat effects###########################################

#relationship between habitat and river/stream categories
unique(trendEstimates$Flow)
subset(trendEstimates,Flow=="running")#not all use rivers
myhabitats <- c(myhabitats,"running")

#single trend models
temp <- ldply(myhabitats,function(x)
  summary(lm(as.formula(paste0("trend~",x)),data=trendEstimates))$coef[2,])
temp$habitat <- myhabitats
temp
#river, bog, then running

summary(lm(trend ~ river + bog,data=trendEstimates))
summary(lm(trend ~ running + bog,data=trendEstimates))
#river is more important than running water

#single change model
temp <- ldply(myhabitats,function(x)
  summary(lm(as.formula(paste0("change~",x)),data=trendEstimates))$coef[2,])
temp$habitat <- myhabitats
temp

summary(lm(change ~ river + bog,data=trendEstimates))
#more marginal        

### range size effects ##############################################################
summary(lm(trend ~ log(nuEuroGrids),data=trendEstimates))
summary(lm(trend ~ nuEuroGrids,data=trendEstimates))


###trend models######################################################################

#use weights or not?
summary(trendEstimates$sd)
qplot(nuEuroGrids,trend,data=trendEstimates,size=sd)
#dont weight the analysis

#use dev time or voltinim
summary(lm(trend~VoltinismProp,data=trendEstimates))#use this one
summary(lm(trend~devTime,data=trendEstimates))

#german of europe range size
summary(lm(trend~germanRange,data=trendEstimates))
summary(lm(trend~nuEuroGrids,data=trendEstimates,))#use thie one

#flow categorisation
summary(lm(trend~Flow,data=trendEstimates))#negative effect of standing
summary(lm(trend~Flow-1,data=trendEstimates))#running are increasing

#build multiple regression model
lm1 <- lm(trend ~ VoltinismProp + nuEuroGrids + 
                  Flight_start + meanTemp + 
                  DispPot + river + bog + Generalism, data=trendEstimates)
summary(lm1)

#final
mytraits <- c(mytraits,"sbog")
lm1 <- lm(trend ~ smeanTemp + sriver + smedHw + sFlight_start, data=trendEstimates)
summary(lm1)

#extract and plot coefficients
coefDF <- addBest(lm1)
coefDF <- subset(coefDF,param!="(Intercept)")
coefDF$Param <- c("temp pref","river use", "wing length",
                  "flight start","voltinism","bog use")

coefDF$Param <- factor(coefDF$Param,
                       levels=rev(c("temp pref","flight start",
                                "voltinism","wing length",
                                "river use","bog use")))

ggplot(coefDF)+
  geom_crossbar(aes(x=Param,y=Estimate,ymin=lowerQ,ymax=upperQ),width=0.3)+
  coord_flip()+
  theme_classic()+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  ylab("Effect on trend")+xlab("")


ggsave("plots/Trait_effect_sizes.png",width=5,height=4)

#save as csv file
write.csv(coefDF,file="plots/Traits_effects.csv",row.names=FALSE)

#as mixed models including genus
library(lme4)
library(lmerTest)
lme1 <- lmer(trend~ VoltinismProp + nuEuroGrids + 
               Flight_start + meanTemp + 
               DispPot + Habitat + 
               Generalism + (1|Genus), data=trendEstimates)
#tmean is the best trait

r.squaredGLMM(lme1)
#genus explains a lot!!

lme1 <- lmer(trend~  meanTemp + river + (1|Genus), data=trendEstimates)
#river not significant after accounting for phylogeny as a random effect

###change models#########################################################

hist(log(trendEstimates$change))

lm1 <- lm(log10(change) ~ VoltinismProp + 
            Flight_start + meanTemp + 
            DispPot + river + 
            Generalism,
          data=trendEstimates)
summary(lm1)
#only mean temp

#final
lm1 <- lm(change ~ meanTemp, data=trendEstimates)
summary(lm1)

#extract and plot coefficients
coefDF <- addBest(lm1,myresponse="change")
coefDF <- subset(coefDF,param!="(Intercept)")
#coefDF$Param <- c("int","voltinism","range size","flight start",
#                  "temp pref","disp potential","habitat","habitat #breadth")

ggplot(coefDF)+
  geom_crossbar(aes(x=param,y=Estimate,ymin=lowerQ,ymax=upperQ),width=0.3)+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  ylab("Effect on total occupancy change")+xlab("Trait")

#as mixed models including genus
library(lme4)
library(lmerTest)
lme1 <- lmer(change ~ VoltinismProp + 
               Flight_start + meanTemp + 
               DispPot + Habitat + 
               Generalism + (1|Genus), data=trendEstimates)
#tmean is the best trait

r.squaredGLMM(lme1)
#genus explains a lot!!

lme1 <- lmer(trend~  meanTemp + (1|Genus), data=trendEstimates)
summary(lme1)

###taxonomy###################################################################################

#suborder
table(trendEstimates$Suborder)
#Anisoptera  Zygoptera 
#50         27

table(trendEstimates$Suborder,trendEstimates$Trend)

chisq.test(table(trendEstimates$Suborder,trendEstimates$Trend))
#ns
prop.test(c(23,10),c(50,27))
prop.test(c(14,11),c(50,27))

#family
sum(is.na(trendEstimates$Family))
table(trendEstimates$Family)
#most from Libellulidae, Coenagrionidae, and Aeshnidae (69%)

chisq.test(table(trendEstimates$Family,trendEstimates$Trend))
#ns

#genus?
table(trendEstimates$Genus)

ggplot(trendEstimates)+
  geom_boxplot(aes(x=Genus,y=trend),alpha=0.5)+
  geom_point(aes(x=Genus,y=trend),colour="blue")+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=0,linetype="dotted")

ggplot(trendEstimates)+
  geom_boxplot(aes(x=Genus,y=change),alpha=0.5)+
  geom_point(aes(x=Genus,y=change),colour="blue")+
  coord_flip()+
  scale_y_log10()+
  theme_bw()+
  geom_hline(yintercept=1,linetype="dotted")


#as random effects in models
library(lme4)
lme1 <- lmer(trend ~ 1 + (1|Genus) + (1|Family), data=trendEstimates)
lme1 <- lmer(trend ~ 1 + (1|Genus), data=trendEstimates)
lme1 <- lmer(trend ~ 1 + (1|Family), data=trendEstimates)

r.squaredGLMM(lme1)

#genus explains some variation
lme1 <- lmer(change ~ 1 + (1|Genus) + (1|Family), data=trendEstimates)
lme1 <- lmer(change ~ 1 + (1|Genus), data=trendEstimates)
lme1 <- lmer(change ~ 1 + (1|Family), data=trendEstimates)

r.squaredGLMM(lme1)

###phylogeny##############################################################

library(ape)
library(picante)

trendEstimates$Suborder <- factor(trendEstimates$Suborder)
trendEstimates$Family   <- factor(trendEstimates$Family)
trendEstimates$Genus    <- factor(trendEstimates$Genus)
trendEstimates$Species  <- factor(trendEstimates$Species)

frm <- ~Suborder/Family/Genus/Species
tr <- as.phylo.formula(frm, data = trendEstimates)
plot(tr)
cophenetic.phylo(tr)#see assumed distances among species

library(nlme)
row.names(trendEstimates) <- trendEstimates$Species
trendEstimates <- trendEstimates[order(match(trendEstimates$Species,tr$tip.label)),]
all(row.names(trendEstimates)==tr$tip.label)
out2 <- compute.brlen(tr,1)
gls1 <- gls(trend ~ smeanTemp + sriver + smedHw + sFlight_start,
            correlation=corPagel(1,tr,fixed=FALSE),data=trendEstimates)
summary(gls1)


gls1 <- gls(trend ~ smeanTemp + sriver + smedHw + sFlight_start,
            correlation=corPagel(1,out2,fixed=FALSE),data=trendEstimates)
gls2 <- gls(trend ~ smeanTemp + sriver + smedHw + sFlight_start,data=trendEstimates)
anova(gls1,gls2)
#no difference

###trait plots#################################################################

#plot trait vs trends

q1<-ggplot(trendEstimates,aes(x=meanTemp,y=trend))+
  geom_point()+
  theme_classic()+
  ylab("Long-term trend")+
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  xlab(expression("Temperature preference "*~degree*C))+
  stat_smooth(fill="blue",colour="blue",method="lm",alpha=0.2)+
  theme(legend.position="none")

q2<-ggplot(trendEstimates,aes(x=Flight_start,y=trend))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  ylab("Long-term trend")+xlab("Start of flight period (month)")+
  stat_smooth(fill="blue",colour="blue",method="lm",alpha=0.2)+
  theme(legend.position="none")

q3<-ggplot(trendEstimates,aes(x=medHw,y=trend))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  ylab("Long-term trend")+xlab("Wing length (mm)")+
  stat_smooth(fill="blue",colour="blue",method="lm",alpha=0.2)+
  theme(legend.position="none")

q4<-ggplot(trendEstimates,aes(x=VoltinismProp,y=trend))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  ylab("Long-term trend")+xlab("Votinism scale")+
  stat_smooth(fill="blue",colour="blue",method="lm",alpha=0.2)+
  theme(legend.position="none")

trendEstimates$River <- ifelse(trendEstimates$river==1,"Yes","No")
q5<-ggplot(trendEstimates,aes(x=River,y=trend))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  ylab("Long-term trend")+xlab("River use")+
  theme(legend.position="none")

trendEstimates$Bog <- ifelse(trendEstimates$bog==1,"Yes","No")
q6<-ggplot(trendEstimates,aes(x=Bog,y=trend))+
  geom_boxplot()+
  theme_classic()+
  geom_hline(yintercept=0,color="black",linetype="dashed")+
  ylab("Long-term trend")+xlab("Bog use")+
  stat_smooth(method="lm")+
  theme(legend.position="none")


plot_grid(q1,q2,q3,q4,q5,q6,align="v",ncol=2)
ggsave("plots/Trait_plots.png",width=5,height=6)

#plot trait vs change

q1<-ggplot(trendEstimates,aes(x=meanTemp,y=log10(change)))+
  geom_point()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Temperature preference")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

q2<-ggplot(trendEstimates,aes(x=Flight_start,y=log10(change)))+
  geom_point()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Start of flight period")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

q3<-ggplot(trendEstimates,aes(x=DispPot,y=log10(change)))+
  geom_point()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Dispersal potential")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

q4<-ggplot(trendEstimates,aes(x=factor(Generalism),y=log10(change)))+
  geom_boxplot()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Habitat breadth")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

q5<-ggplot(trendEstimates,aes(x=factor(river),y=log10(change)))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  ylab("Total occupancy change")+xlab("River use")+
  theme(legend.position="none")

q6<-ggplot(trendEstimates,aes(x=factor(bog),y=log10(change)))+
  geom_boxplot()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Bog")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

q7<-ggplot(trendEstimates,aes(x=VoltinismProp,y=log10(change)))+
  geom_point()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Votinism")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

q8<-ggplot(trendEstimates,aes(x=Overwintering,y=log10(change)))+
  geom_boxplot()+
  theme_bw()+
  ylab("Total occupancy change")+xlab("Overwintering stage")+
  stat_smooth(method="lm")+
  theme(legend.position="none")

plot_grid(q1,q2,q3,q4,q5,q6,q7,q8,align="v",ncol=2)

###Community-level################################################

#plot for David
out <- readRDS("modelSummary_Odonata_adult_Bav.rds")
out <- getCodeFromFile(out,myfile="out_nuSpecies_adult_Bav_")
out <- subset(out,grepl("psi.fs",out$Param))
out$ParamNu <- as.numeric(sub(".*\\[([^][]+)].*", "\\1", out$Param))
out$Year <- out$ParamNu+1979

#for each species and year get 1000 possible occupancies

outRandom <- ddply(out,.(Year,Species),summarise,
                   vals = rnorm(1000,mean,sd),
                   sim = 1:length(vals))

#get mean occupancy per year
outRandomMean <- ddply(outRandom,.(Year,sim),summarise,mean=mean(vals))
outRandomQ <- ddply(outRandomMean,.(Year),summarise,
                    medianM = median(mean),
                    lowerQ = quantile(mean,0.025),
                    upperQ = quantile(mean,0.975))

ggplot(data=outRandomQ)+
  theme_classic()+xlab("Year")+ylab("Mittlere Vorkommenswahrscheinlichkeit")+
  geom_line(aes(x=Year,y=medianM))+
  geom_ribbon(aes(x=Year,ymin=lowerQ,ymax=upperQ),alpha=0.4)

ggsave("NULplot.tiff",dpi=600)

####cwm#################################################################

load("randomMatrix.RData")

randomMatrixM<-melt(randomMatrix,id=c("Year","Species"))
randomMatrixM<-arrange(randomMatrixM,Species,Year)
randomMatrixM$Species[!randomMatrixM$Species %in% alltraits$Species]
randomTrends <- merge(randomMatrixM,alltraits,by="Species")

#european distribution
tmeansMeans <- ddply(randomTrends,.(Year,variable),summarise,
                     tmean = weighted.mean(nuEuroGrids,value))
tmeansMeans <- ddply(tmeansMeans,.(Year),summarise,
                     my.mean = mean(tmean),
                     lowerCI = quantile(tmean,0.025),
                     upperCI=quantile(tmean,0.975))

tmeansMeans$Year <- tmeansMeans$Year+1979
g1 <- ggplot(tmeansMeans)+
  geom_line(aes(x=Year,y=my.mean))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_classic()+
  ylab("Mean range size")+
  theme(legend.position="none")

####species richness#################################################

load("randomMatrix.RData")

out <- ddply(randomMatrix,.(Year),function(x){
  numcolwise(sum)(x)})

#get mean and 95% CI across communities
out$meanRichness<-apply(out[,2:1001],1,median)
out$lowerRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.025))
out$upperRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.975))

out$Year <- out$Year+1979
g2 <- ggplot(out)+
  geom_line(aes(x=Year,y=meanRichness))+
  geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness),alpha=0.3)+
  theme_classic()+
  ylab("Mean species richness")

###disimmilarity##############################################

load("randomMatrix.RData")

library(vegan)
randomMatrixM<-melt(randomMatrix,id=c("Year","Species"))
randomMatrixM<-arrange(randomMatrixM,Species,Year)

nuSites <- 500
randomMatrixM$predProp <- round(randomMatrixM$value*nuSites)

out <- ddply(randomMatrixM,.(variable),function(x){
  moo2 <- acast(x,Year~Species,value.var="predProp")
  diss <- as.numeric(as.matrix(vegdist(moo2,method="bray",binary=FALSE))[,1])
  year <- names(as.matrix(vegdist(moo2,binary=TRUE))[,1])
  data.frame(year,diss)
})

out$year <- as.numeric(as.character(out$year))+1979
out <- subset(out,year>1980)

#take average across all
out <- ddply(out,.(year),summarise,
             meanD = mean(diss),
             lowerCI = quantile(diss,0.025),
             upperCI = quantile(diss,0.975))

out$Year <-as.numeric(as.character(out$year))
ggplot(out)+
  geom_line(aes(x=Year,y=meanD))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_bw()+
  ylab("Average dissimilarity")

###diversity#######################################################################

data(BCI)
head(BCI)
dim(BCI)
diversity(BCI)

H <- diversity(acast(BCI,Year~Species,value.var="predProp"))

out <- ddply(randomMatrixM,.(variable),function(x){
  moo2 <- acast(x,Year~Species,value.var="predProp")
  div <- diversity(moo2) 
  year <- sort(unique(x$Year))
  data.frame(year,div)
})

out$year <- as.numeric(as.character(out$year))+1979
out <- subset(out,year>1979)

#take average across all
out <- ddply(out,.(year),summarise,
             meanD = mean(div),
             lowerCI = quantile(div,0.025),
             upperCI = quantile(div,0.975))

out$Year <-as.numeric(as.character(out$year))

g3 <- ggplot(out)+
  geom_line(aes(x=Year,y=meanD))+
  geom_ribbon(aes(x=Year,ymin=lowerCI,ymax=upperCI),alpha=0.3)+
  theme_classic()+
  ylab("Mean diversity")

plot_grid(g2,g3,g1,nrow=1)
ggsave("plots/Community_level.png",width=7,height=2)

###func diversity#########################################################

library(FD)

#gowdis computes the Gower dissimilarity from different trait types (continuous, ordinal, nominal, or binary)
trendEstimates$Species <- as.character(trendEstimates$Species)
trendEstimates <- arrange(trendEstimates,Species)
myTraits <- trendEstimates[,allTraits]
row.names(myTraits) <- trendEstimates$Species
ex1 <- gowdis(myTraits)

myMatrix <- dcast(annualDF,Year~Species,value.var="mean")
#Rows are sites and species are columns
ex3 <- dbFD(ex1,myMatrix[,-1])#takes 1 min

#Interesting output
fdDF <- data.frame(Year=1980:2016,
                   funcDiv=ex3$FDiv,
                   funcDisp=ex3$FDis)

qplot(Year,funcDiv,data=fdDF,geom="line")
qplot(Year,funcDisp,data=fdDF,geom="line")

ggplot(fdDF)+
  geom_line(aes(x=Year,y=funcDiv))+
  #geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness),alpha=0.3)+
  theme_classic()+
  ylab("Average functional diversity")


#apply as a function to each iteraction
load("randomMatrix.RData")
randomMatrixM <- melt(randomMatrix,id=c("Year","Species"))
randomMatrixM <- arrange(randomMatrixM,Species,Year)
randomMatrixM$Species[!randomMatrixM$Species %in% alltraits$Species]

fdMeans <- ddply(subset(randomMatrixM,variable %in% 1:10),.(variable),
                 function(x){
              
                     myMatrix <- dcast(x,Year~Species,value.var="value")
                     ex3 <- dbFD(ex1,myMatrix[,-1])
                     
                     #Interesting output
                     fdDF <- data.frame(Year=1980:2016,
                                        funcDiv=ex3$FDiv,
                                        funcDisp=ex3$FDis)
                     return(fdDF)
                     
                     })


g3 <- ggplot(fdDF)+
  geom_line(aes(x=Year,y=funcDisp))+
  #geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness),alpha=0.3)+
  theme_classic()+
  ylab("Functional diversity")

plot_grid(g2,g1,g3,nrow=1)
ggsave("plots/Community_level.png",width=7,height=2)

###end###############################################################