#run the modelSummaries_stan

#the script asks whether species about:
#the distribution of changes - skewed or uniform

#initial subsetting
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))
allspecies <- sort(unique(modelSummaries_Limits$Species))

### functions for analysis #####

getGridChange <- function(modelSummaries_Limits){
  
  gridChange <- modelSummaries_Limits %>%
    dplyr::select(c(Year,Species,MTB,mean)) %>%
    tidyr::pivot_wider(.,names_from="Year",values_from="mean") %>%
    dplyr::group_by(Species,MTB) %>%
    janitor::clean_names() %>%
    dplyr::rename(.,Species = species) %>%
    dplyr::summarise(change = log(x2016+0.01/x1990+0.01))
  
    return(gridChange)
  
}


trimChange <- function(gridChange, myquantile = 0.95){
  
  upper <- quantile(gridChange$change, myquantile)
  
  gridChange$change <- ifelse(gridChange$change>upper,upper,gridChange$change)
  
  lower <- quantile(gridChange$change, 1-myquantile)
  
  gridChange$change <- ifelse(gridChange$change<lower,lower,gridChange$change)
  
  return(gridChange)
  
}


sortSpecies <- function(gridChange){
  
  #get mean change of a species
  meanChange <- gridChange %>%
                group_by(Species) %>%
                summarise(meanChange = mean(change)) %>%
                arrange(desc(meanChange))
  
  gridChange$Species <- factor(gridChange$Species,
                               levels = meanChange$Species)
  
  return(gridChange)

}


plotChange <- function(gridChange, type=1){
  
  require(ggridges)
  
  if(type==1){
  
    ggplot(gridChange, aes(x = change, y = Species)) + geom_density_ridges() +
    theme_few() +
    geom_vline(xintercept = 0, linetype = "dashed")
    
  }else if(type==2){
  
    ggplot(gridChange, aes(x = change, y = Species, fill = stat(x))) +  
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    theme_few() +
    theme(axis.text.y = element_text(size = 8)) +
    scale_fill_viridis_c(name = "Occupancy change", option = "C") +
    geom_vline(xintercept = 0, linetype = "dashed")
    
  }
}


summariseChange <- function(modelSummaries_Limits, myspecies = NULL){

  x <- gridChange$change[gridChange$Species==myspecies]
    
  #shape, centre, spread
  
  #shape: skewness and modality/peakedness - shape
  skewness = moments::skewness(x)#direction of skew
  kurtosis = moments::kurtosis(x)#heaviness of tails/peakedness
  modality = as.numeric(diptest::dip.test(x)$statistic)
  
  #centre: median/mean
  mymean = mean(x)
  mymedian = median(x)
  
  #spread - sd and interquartile range and 95% range
  mysd = sqrt(var(x))
  myinterquartile = as.numeric(quantile(x,0.75) - quantile(x,0.25))
  
  #outliers
  myoutliers = length(boxplot.stats(x)$out)
  
  data.frame(Species = myspecies, skewness, kurtosis, mean = mymean, abs_median = abs(mymedian),
             sd = mysd, interquartile = myinterquartile, outliers = myoutliers, modality,
             direction = ifelse(mymedian>0,"winner","loser"))
  
}


#### applying functions #####

summaryData <- allspecies %>%
                map_dfr(function(x){
                  summariseChange(modelSummaries_Limits, myspecies = x)}) %>%
                  pivot_longer(!c(Species,direction),
                               names_to = "change_type",
                               values_to = "value")


ggplot(filter(summaryData, !change_type %in% c("mean","sd"))) +
  geom_violin(aes(x = direction, y = value))+
  facet_wrap( ~ change_type, scales = "free")

#winners: higher range, lower kurtosis, higher modality, fewer outliers, more left-skewed
#simple terms: 
# large magnitude of change,
# changes more variable over range, 
# but stronger peak change,
# more modality
# fewer outliers
# more likely to be left skewed i.e., most changes at the largest changes


#losers: smaller range, higher kurtosis, lower modality, more outliers, slightly more left skewed
#simple terms:
#small magnitude of change
#smaller range of change
#less peakedness
#less modality
#more outliers
#both left and right skewed
