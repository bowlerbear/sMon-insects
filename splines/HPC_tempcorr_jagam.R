#load data
#run analysis_HPC_nation_spline_mtb

myfolder <- "/data/idiv_ess/Odonata/"
#myfolder <- "splines"

#need to load in y
bugs.data <- readRDS(paste(myfolder,"bugs.data.rds",sep="/"))

#set model parameters
ni <- 4000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

#fit model
library(rjags)
library(jagsUI)

#get core info
n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) 
#n.cores = 3

modelfile=paste(myfolder,"BUGS_sparta_nation_naturraum_spline.txt",sep="/")

#get inits
zst <- readRDS(paste(myfolder,"zst.rds",sep="/"))
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
if(class(out)=="jagsUI"){
  saveRDS(out$summary,file="outSummary_tempcorr_spline.rds")
}else{
  saveRDS(out,file="out_tempcorr_spline.rds")
}
