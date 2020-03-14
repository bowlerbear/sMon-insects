#get species for this task
myfolder <- "Odonata_adult_nation_state"

#open folder and get list of files within
folder1 <- paste0("/work/bowler/",myfolder)
folder2 <- list.files(folder1)[1]#should only be one folder in there

#read in each list and pull out the model summary
modelFiles <- list.files(paste(folder1,folder2,sep="/"))
library(plyr)
modelSummary <- ldply(modelFiles, function(x){
  temp <- readRDS(paste(folder1,folder2,x,sep="/"))
  #get summary of parameters
  tempSummary <- data.frame(temp)
  tempSummary$Param <- as.character(row.names(tempSummary))
  tempSummary$File <- x
  return(tempSummary)
})

#save the folder name as a column
saveRDS(modelSummary,file=paste0("modelSummary_",myfolder,".rds"))
