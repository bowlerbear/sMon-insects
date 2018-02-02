library(lubridate)
library(sparta)

ortho_mv<-readRDS("raw-data\\ortho_mv_.rds")


## filter for data with only years
ortho_mv$day<- lubridate:: day(ortho_mv$Datum)
ortho_mv$month<- lubridate:: month(ortho_mv$Datum)

###filter out 01-01
ortho_mv<- ortho_mv[!((ortho_mv$day==1) & (ortho_mv$month==1)),]


##fliter out year NA
ortho_mv<- ortho_mv[!is.na(ortho_mv$Jahr),]
ortho_mv<- ortho_mv[!is.na(ortho_mv$day),]
ortho_mv<- ortho_mv[!is.na(ortho_mv$month),]


# reduce the data to the subset from 1990
ortho_mv <- subset(ortho_mv, Jahr %in% 1990:2003)

#combine $ANZAHL and $ANZSPEZ
#different Method a problem??? $METHODE

#length(unique(ortho_mv$Species)) #39
#length(unique(ortho_mv$Family))  #5


dataDiagnostics(taxa = ortho_mv$Species,
                site = ortho_mv$MTB_Q,
                time_period = ortho_mv$Jahr)

