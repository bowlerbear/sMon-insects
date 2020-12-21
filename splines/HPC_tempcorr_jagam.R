#load data
#need to load in y
bugs.data <- readRDS("/data/idiv_ess/Odonata/bugs.data_NAs.rds")

#set model parameters
ni <- 40000   ;   nb <- 20000   ;   nt <- 20   ;   nc <- 3

#fit model
library(rjags)
library(R2WinBUGS)
library(jagsUI)

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
#set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
#n.cores = 3

modelfile="/data/idiv_ess/Odonata/BUGS_sparta_nation_naturraum_spline.txt"

#get inits
zst <- readRDS("/data/idiv_ess/Odonata/zst.rds")
inits <- function(){list(z = zst)}

#specifiy effort
effort = "shortList"
bugs.data$Effort <- bugs.data[[effort]]

#specify parameters to monitor
params <- c("psi")

Sys.time()
#run model
out <- jags(bugs.data, inits=inits, params, modelfile, n.thin=nt,
            n.chains=n.cores, n.burnin=nb,n.iter=ni,parallel=T)
Sys.time()

#save as output file
saveRDS(out$summary,file="outSummary_tempcorr_spline_NAs.rds")
