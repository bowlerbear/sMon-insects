#load data
siteInfo <- readRDS("/data/idiv_ess/Odonata/siteInfo_NAs.rds")

#need to load in y
bugs.data <- list(y = siteInfo$obs,
                  n = length(siteInfo$obs))

#add jagam objects
jags.ready <- readRDS("/data/idiv_ess/Odonata/jags.ready_NAs.rds")
bugs.data$X = jags.ready$jags.data$X
bugs.data$S1 = jags.ready$jags.data$S1
bugs.data$zero = jags.ready$jags.data$zero
                  
#set model parameters
ni <- 50000   ;   nb <- 25000   ;   nt <- 20   ;   nc <- 3

#fit model
library(rjags)
library(R2WinBUGS)
library(jagsUI)

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
#set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 
#n.cores = 3

modelfile="/data/idiv_ess/Odonata/jagam.txt"

#specify parameters to monitor
params <- c("mu")

Sys.time()
#run model
out <- jags(bugs.data, inits=NULL, params, modelfile, n.thin=nt,
            n.chains=n.cores, n.burnin=nb,n.iter=ni,parallel=T)
Sys.time()

#save as output file
saveRDS(out$summary,file="outSummary_static_spline_NAs.rds")
