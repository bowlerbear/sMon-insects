#load the relational table of task ids and species
speciesTaskID <- read.delim(paste0("/data/idiv_ess/Odonata/speciesTaskID_adult.txt"),as.is=T)

#get task id
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")) 

#get species for this task
myspecies <- speciesTaskID$Species[which(speciesTaskID$TaskID==task.id)]

print(myspecies)


