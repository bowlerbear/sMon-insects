library(tidyverse)
library(ggthemes)

#get list of stan models
modelDirectory <- "model-outputs/=Odonata_stan_spline"
stanFiles <- list.files(modelDirectory) %>% str_subset("m_fit")

#function to apply to each file
readStanModel <- function(file, get="psi"){
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl(get,temp$Param))
  
  #get species name
  temp$File <- file
  temp$Species <- gsub("m_fit_summary_spacetime_","",temp$File)
  temp$Species <- gsub(".rds","",temp$Species)
  
  return(temp)
}

#apply to all files
modelSummaries <- stanFiles %>%
  map_df(readStanModel) %>%
  mutate(siteIndex = parse_number(Param))


#get site and year information
siteInfo_NAs <- readRDS("stan/siteInfo_NAs.rds") %>%#see the HPC_spline_stan_MTB script
                select(!c("Species","SpeciesOrig"))
modelSummaries <- left_join(modelSummaries,siteInfo_NAs, by="siteIndex")


#analyse nationwide time series for each species
nuMTBs <- length(unique(siteInfo_NAs$MTB))
myspecies <- sort(unique(modelSummaries$Species))
annualTS <- modelSummaries %>%
            group_by(Species,Year) %>%
            summarise(total = sum(mean)/nuMTBs)

ggplot(annualTS)+
  geom_line(aes(x=Year, y=total))+
  facet_wrap(~Species)+
  theme_bw()
#very smooth time-series... maybe try with more knots

#plot spatial maps
#just compare the start and the end for each species
modelSummaries_Limits <- filter(modelSummaries, Year %in% c(1990,2016))
ggplot(filter(modelSummaries_Limits, Species %in% myspecies[1:10]))+
  geom_point(aes(x=x_MTB, y=y, colour=mean))+
    facet_grid(Year~Species)+
  scale_color_viridis_c()

ggplot(filter(modelSummaries_Limits, Species %in% myspecies[11:20]))+
  geom_point(aes(x=x_MTB, y=y, colour=mean))+
  facet_grid(Year~Species)+
  scale_color_viridis_c()

ggplot(filter(modelSummaries_Limits, Species %in% myspecies[21:30]))+
  geom_point(aes(x=x_MTB, y=y, colour=mean))+
  facet_grid(Year~Species)+
  scale_color_viridis_c()
