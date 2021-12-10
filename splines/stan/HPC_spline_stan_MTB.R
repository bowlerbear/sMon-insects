#code to work out number of knots for each species
library(tidyverse)

myfolder <- "/data/idiv_ess/Odonata"#HPC
#myfolder <- "/data/dbowler/Odonata/data"#rstudio

### load in the full dataset ####

adultData <- readRDS(paste(myfolder,"adultData.rds",sep="/"))
#load(paste(myfolder,"mtbqsDF.RData",sep="/"))
mtbsDF_extended <- readRDS(paste(myfolder,"MTB_extendedpoints.rds",sep="/"))
names(mtbsDF_extended)[1] <- "MTB"
mtbsDF <- subset(mtbsDF_extended, !is.na(MTB))

# subset  

df <- subset(adultData, Year>=1990  & Year<2017)
#table(df$Month)
df <- subset(df, Month %in% 4:10)

#get species

speciesTaskID <- read.delim(paste(myfolder,"speciesTaskID_adult.txt",sep="/"),as.is=T)
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)]
#myspecies <- "Sympetrum danae" #test case

#get MTB
df$MTB <- sapply(as.character(df$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,(len-1))})

### function to get data for a species ####

all(df$MTB %in% mtbsDF$MTB)

#mtbsDF$CoarseNatur[is.na(mtbsDF$CoarseNatur)] <- mtbsDF$MTB_CoarseNatur[is.na(mtbsDF$CoarseNatur)]
df$CoarseNatur <- mtbsDF$Coarsenaturraum[match(df$MTB,mtbsDF$MTB)]
sum(is.na(df$CoarseNatur))

#mtbsDF$MTB_MidNatur <- gdata::trim(mtbsDF$MTB_MidNatur)
#df$MidNatur <- mtbsDF$MTB_MidNatur[match(df$MTB,mtbsDF$MTB)]
#sum(is.na(df$MidNatur))

#define a visit
df$visit <- paste(df$MTB,df$Date,df$Beobachter,sep="_")

#get occurence matrix  - detection and non-detection
getOccurrenceMatrix<-function(df){
  out <- reshape2::acast(df,visit~Species,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
  out[out>0]<-1
  return(out)
}
occMatrix <- getOccurrenceMatrix(df)

#get list length
getListLength<-function(df){
  df %>% 
    dplyr::group_by(visit,Date,MTB) %>%
    dplyr::summarise(nuSpecies=length(unique(Species)),
                    nuRecords=length(Species)) %>%
    arrange(.,visit)
}
listlengthDF <- getListLength(df)

#check rows of occuMatrix match visits
all(listlengthDF$visit==row.names(occMatrix))

#add on some indices
listlengthDF$Date <- as.Date(listlengthDF$Date)
listlengthDF$Year <- lubridate::year(listlengthDF$Date)
listlengthDF$yday <- lubridate::yday(listlengthDF$Date)
listlengthDF$yearIndex <- as.numeric(factor(listlengthDF$Year))

#add natur raum
listlengthDF$CoarseNaturraum <- mtbsDF$Coarsenaturraum[match(listlengthDF$MTB,mtbsDF$MTB)]
listlengthDF$cnIndex <- as.numeric(factor(listlengthDF$CoarseNaturraum))
subset(listlengthDF,is.na(CoarseNaturraum))

# listlengthDF$MidNaturraum <- mtbqsDF$MTB_MidNatur[match(listlengthDF$MTB,mtbsDF$MTB)]
# listlengthDF$mnIndex <- as.numeric(factor(listlengthDF$MidNaturraum))
# subset(listlengthDF,is.na(MidNaturraum))

#get other effort variables
listlengthDF$singleList <- ifelse(listlengthDF$nuSpecies==1,1,0)
listlengthDF$shortList <- ifelse(listlengthDF$nuSpecies%in%2:3,1,0)
listlengthDF$longList <- ifelse(listlengthDF$nuSpecies>3,1,0)

#order data
listlengthDF <- arrange(listlengthDF,visit)
all(row.names(occMatrix)==listlengthDF$visit)

#add on species data to the listlength object
listlengthDF$Species <- occMatrix[,myspecies]
table(listlengthDF$Species)

### expand to fill mtbq ####

df <- listlengthDF

nrow(df)

#sort site indices
df$siteIndex <- as.numeric(as.factor(paste0(df$MTB,df$Year)))
df <- arrange(df,siteIndex,visit)
nrow(df)

#add on geographic coordinates
df$x <- mtbsDF$x[match(df$MTB, mtbsDF$MTB)]
df$y <- mtbsDF$y[match(df$MTB, mtbsDF$MTB)]

#make values smaller
medX <- median(df$x/10000)
medY <- median(df$y/1000000)
df$x <- df$x/10000 - medX
df$y <- df$y/1000000 - medY

# Aggregate observed data
site_df <- df %>%
  group_by(siteIndex,x,y,MTB,Year) %>%
  summarise(Species = max(Species, na.rm=T))
table(site_df$Species)

#few checks
nrow(site_df)#27606
nrow(site_df)==length(unique(df$siteIndex))#should be TRUE

#qplot(x, y, color = Species, data = site_df, alpha=0.5,facets=~Year)+theme(legend.position="none")

# expand to whole range 
siteInfo_NAs <- mtbsDF_extended
nuSites <- nrow(siteInfo_NAs)
siteInfo_NAs <- siteInfo_NAs %>% slice(rep(1:n(), each = length(unique(df$Year))))
siteInfo_NAs$Year <- rep(sort(unique(df$Year)),nuSites)
siteInfo_NAs$yearIndex <- rep(1:length(unique(df$Year)),nuSites)

#sort species data
siteInfo_NAs$Species <- siteInfo_NAs$SpeciesOrig  <- site_df$Species[match(interaction(siteInfo_NAs$MTB,siteInfo_NAs$Year),
                                                                           interaction(site_df$MTB,site_df$Year))]
siteInfo_NAs$Species[is.na(siteInfo_NAs$Species)] <- 0
table(siteInfo_NAs$Species)

#keep original coordinates
siteInfo_NAs$x_MTB <- siteInfo_NAs$x
siteInfo_NAs$y_MTB <- siteInfo_NAs$y

#sort coordinates
siteInfo_NAs$x <- siteInfo_NAs$x/10000 - medX
siteInfo_NAs$y <- siteInfo_NAs$y/1000000 - medY

#arrange file
siteInfo_NAs$siteIndex <- as.numeric(as.factor(paste0(siteInfo_NAs$MTB,siteInfo_NAs$Year)))
siteInfo_NAs <- arrange(siteInfo_NAs, siteIndex)
#saveRDS(siteInfo_NAs, file="splines/siteInfo_NAs.rds")

### fit model ####

library(rstan)
library(brms)
library(mgcv)

#use make_stancode to get model code

#original model
#mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex)),
#                                           data = siteInfo_NAs, 
#                                           family = bernoulli())


#more wiggly model(v2)
#mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, d = c(2,1), k = c(5,8))),
#                                           data = siteInfo_NAs, 
#                                           family = bernoulli())
 

#v3                                           
#mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, d = c(2,1), 
#                                                           k = c(5,5))),data = siteInfo_NAs, 
#                                           family = bernoulli())


#v4
# mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, d = c(2,1), 
#                                                            bs = c("ts", "cs"),
#                                                            k = c(10,10))),data = siteInfo_NAs, 
#                                                            family = bernoulli())

#v5
# mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, d = c(2,1), 
#                                                            bs = c("ts", "cs"),
#                                                            k = c(15,10))),data = siteInfo_NAs, 
#                                                           family = bernoulli())


#v6
# mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, 
#                                                            d = c(2,1),k = c(8,5))),
#                                            data = siteInfo_NAs,
#                                            family = bernoulli())
                                                           

#v7 - spatial only model
# mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, k=15)),
#                               data = siteInfo_NAs, 
#                               family = bernoulli())


#v8
#mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, k=8)),
#                                           data = siteInfo_NAs, 
#                                           family = bernoulli())
                                           

#v9
#mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, k=12)),
#                                           data = siteInfo_NAs, 
#                                           family = bernoulli())



#v10 include covariates on the detection model
# effort term - short and long list
# yday and yday2
# extended MTB frame
#mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, k=8)),
#                                           data = siteInfo_NAs, 
#                                           family = bernoulli())


#v11 include covariates on the detection model
# effort term - short and long list
# yday and yday2
# extended MTB frame
mydata_complete <- mydata <- make_standata(bf(Species ~ t2(x, y, yearIndex, k=7)),
                                           data = siteInfo_NAs, 
                                           family = bernoulli())


names(mydata_complete) <- sapply(names(mydata_complete), 
                                 function(x) paste0("complete_",x))

# subset code to sites with data -----------------------------------------------------------------

#define sites with data
presentData <- !is.na(siteInfo_NAs$SpeciesOrig)
table(presentData)

#subset remaining data
mydata$Y <- mydata_complete$complete_Y[presentData]
mydata$X <- mydata_complete$complete_X[presentData,]
mydata$Xs <- mydata_complete$complete_Xs[presentData,]
mydata$Zs_1_1 <- mydata_complete$complete_Zs_1_1[presentData,]
mydata$Zs_1_2 <- mydata_complete$complete_Zs_1_2[presentData,]
mydata$Zs_1_3 <- mydata_complete$complete_Zs_1_3[presentData,]
mydata$Zs_1_4 <- mydata_complete$complete_Zs_1_4[presentData,]
mydata$Zs_1_5 <- mydata_complete$complete_Zs_1_5[presentData,]
mydata$Zs_1_6 <- mydata_complete$complete_Zs_1_6[presentData,]
mydata$Zs_1_7 <- mydata_complete$complete_Zs_1_7[presentData,]

# Occupancy states ------------------------------------------

# define a design matrix for site-level occupancy
n_site <- length(unique(df$siteIndex))

#check data lengths match
length(mydata$Y)==n_site 

#lets fit intercept only model 
X_psi <- matrix(c(rep(1, n_site)))
m_psi <- ncol(X_psi)    # m_psi is the number of columns in the site level design matrix

#include naturraum as a fixed effect?? interaction with year
#X_psi <- mydata$X
#m_psi <- ncol(X_psi) 

# Survey data --------------------------

# get observations
df <- arrange(df,siteIndex,visit)
table(df$Species)

# determine number of surveys per site
n_survey <- as.numeric(tapply(df$visit,df$siteIndex,length))
summary(n_survey)
total_surveys <- sum(n_survey)
total_surveys

# define a survey-level design matrix for detection probabilities
#intercept only
X_p <- matrix(c(rep(1, total_surveys)))
m_p <- ncol(X_p)

#or covariate model
#effort terms
X_p_new <- df[,c("singleList","shortList")]
#phenology terms
X_p_new$yday <- as.numeric(scale(df$yday))
X_p_new$yday2 <- as.numeric(scale(df$yday^2))
#combine all
X_p <- as.matrix(cbind(X_p,X_p_new))
m_p <- ncol(X_p)

#save for examining predictions
#saveRDS(unique(df[,c("MTB","Year")]), file="splines/fitDF.rds")

# indices ---------------------------------------------------------

# get start and end indices to extract slices of y for each site
start_idx <- rep(0, n_site)
end_idx <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    site_indices <- which(df$siteIndex == i)
    start_idx[i] <- site_indices[1]
    end_idx[i] <- site_indices[n_survey[i]]
  }
}

# create vector of whether any positive observations were seen at each site
any_seen <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    any_seen[i] <- max(df$Species[start_idx[i]:end_idx[i]])
  }
}

summary(any_seen)
summary(n_survey)
summary(start_idx)
summary(end_idx)

# Bundle data for Stan ----------------------------------------------------

stan_d <- list(n_site = length(unique(df$siteIndex)), 
               m_psi = m_psi, #number of parameters
               X_psi = X_psi, #design matrix for occupancy 
               total_surveys = nrow(df), 
               m_p = m_p, #number of parameters
               X_p = X_p, #design matrix for detection
               site = df$siteIndex, #site indices
               y = df$Species, #observation data
               start_idx = start_idx, #start indices
               end_idx = end_idx, # end indicies
               any_seen = any_seen, #anything seen at a site
               n_survey = n_survey) #nu of surveys at each site

# add on spline elements to the data list
names(mydata)
stan_d <- c(stan_d,mydata)

#add on spline elements for prediction
stan_d <- c(stan_d,mydata_complete)

# Fit model ---------------------------------------------------------------

#select model
m_init <- stan_model(paste(myfolder,
                           'bernoulli-occupancy-long-spline-complete_space_time_v9.stan',sep="/"))

#get cores
# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)


m_fit <- sampling(m_init, 
                  chains = 4,
                  iter = 2000,
                  data = stan_d, 
                  pars = c('beta_p','psi'))#'beta_psi'

saveRDS(summary(m_fit)$summary,file=paste0("m_fit_summary_spacetime_",myspecies,".rds"))
