library(lubridate)
library(sparta)

ortho_s<-readRDS("raw-data\\ortho_Sachs.rds")
## filter for data with only years
ortho_s$day<- lubridate:: day(ortho_s$Datum)
ortho_s$month<- lubridate:: month(ortho_s$Datum)

##fliter out year NA
ortho_s<- ortho_s[!is.na(ortho_s$Jahr),]
ortho_s<- ortho_s[!is.na(ortho_s$day),]
ortho_s<- ortho_s[!is.na(ortho_s$month),]
###filter out 01-01
ortho_s<- ortho_s[!((ortho_s$day==1) & (ortho_s$month==1)),]

# reduce the data to the subset from 1990 to 2014
#Data available up to 2017, very low number of records - Data missing?
ortho_s <- subset(ortho_s, Jahr %in% 1990:2014)

# Change Gomphocerinae subfam to the right fam in fam column -> Acrididae
ortho_s$Family<-gsub("Gomphocerinae","Acrididae", ortho_s$Family,fixed=T)

# Kick XY spec.
ortho_s<-ortho_s[!ortho_s$Species=="Tetrix spec.",]
ortho_s<-ortho_s[!ortho_s$Species=="Chorthippus spec.",]

#length(unique(ortho_s$Species)) #59
#length(unique(ortho_s$Family))  #7 - 2 additional families compared with mv Rhaphidophoridae & Mantidae


dataDiagnostics(taxa = ortho_s$Species,
                site = ortho_s$MTB_Q,
                time_period = ortho_s$Jahr)
