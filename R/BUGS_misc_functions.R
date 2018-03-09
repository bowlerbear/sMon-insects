#load libraries
library(unmarked)
library(rjags)
library(AHMbook)
library(R2WinBUGS)
library(jagsUI)
library(ggplot2)

# JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

# Default parameters
ni <- 6000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

my_ggs_density<-function(D){
  library(ggplot2)
  f <- ggplot(D, aes(x = value))
  f <- f + geom_density(alpha = 0.3,fill="steelblue")+facet_wrap(~Chain)
  return(f)
}

my_ggs_ppmean<-function(D,outcome){
  ppM <- D %>% dplyr::group_by(Iteration) %>% dplyr::summarize(m = mean(value))
  f <- ggplot()+ geom_density(data=ppM,aes(x=m),alpha=0.5,
                              fill="steelblue",
                              colour="steelblue") +  
    geom_vline(xintercept = mean(outcome,na.rm=T))
  return(f)
}