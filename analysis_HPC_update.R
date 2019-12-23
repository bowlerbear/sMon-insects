########################################################################################
suppressMessages(library(rjags))
suppressMessages(library(R2WinBUGS))
suppressMessages(library(jagsUI))
suppressMessages(library(lubridate))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))


#load the relational table of task ids and species
speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_adult.txt"),as.is=T)
#get task id
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 
#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)] 

#set stage
stage="adult"

#set seed
set.seed(3)

#number of MCMC samples
niterations = 10000

########################################################################################

#fit model
library(rjags)
library(R2WinBUGS)
library(jagsUI)

#JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

#get core info
n.cores = as.integer(Sys.getenv("NSLOTS", "1")) 

###########################################################################################

#choose model file
myfolder <- "/work/bowler/Odonata_adult_nation_naturraum_5000iter/6057269"
myfiles <- list.files(myfolder)
myfile <- myfiles[grepl(myspecies,myfiles)]
out <- readRDS(paste(myfolder,myfile,sep="/"))

Sys.time()

#update model
out <- update(out,n.iter=niterations)

Sys.time()

#save as new output file
saveRDS(out,file=paste0("out_dynamic_nation_naturraum_",stage,"_",myspecies,".rds"))

########################################################################################
