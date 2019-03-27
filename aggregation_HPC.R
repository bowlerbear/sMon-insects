#each task on each folder

#get task id
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1")) 

#get species for this task
folderTaskID <- read.delim("/data/idiv_ess/Odonata/folderTaskID.txt",as.is=T)
myfolder <- folderTaskID$Folder[which(folderTaskID$TaskID==task.id)] 

#open folder and get list of files within
folder1 <- paste0("/work/bowler/",myfolder)
folder2 <- list.files(folder1)[1]#should only be one folder in there

#read in each list and pull out the model summary
modelFiles <- list.files(paste(folder1,folder2,sep="/"))
library(plyr)
modelSummary <- ldply(modelFiles, function(x){
  temp <- readRDS(paste(folder1,folder2,x,sep="/"))
  #get summary of parameters
  tempSummary <- data.frame(temp$summary)
  tempSummary$Param <- as.character(row.names(tempSummary))
  tempSummary$File <- x
  return(tempSummary)
})

#save the folder name as a column
saveRDS(modelSummary,file=paste0("modelSummary_",myfolder,".rds"))