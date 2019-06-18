#powerAnalysis:

#########################################################################################################
#- test the ability of the predicted occupancy models to estimate trends of a given amount

#for each species, how strong would the real trend have to be before we could detect it, given the SE of the occupancies

#given this SE and mean occupancy, we could only detect a change in occupancy of X

#get this for each species??

###########################################################################################################

#YearA, 80 (10) 
#YearB, 100 (20)

#given the SE, what sort of trend would be difficult from zero at 5%

#sample size required to detect a difference in occupancy with a given power

#MacKenzie & Royle 2005;
#Guillera-Arroita, Ridout & Morgan 2010).

#We provide an approximate
#expression to calculate power and derive a closed-formula that
#conveniently allows the number of sites that need to be sampled
#to be determined with just a few simple calculations, while
#accounting for species detectability.

#we provide R code for conducting
#power analysis, both based on the formula and via
#simulations

#G - power to detect a difference in occupancy

#Let R be the proportional difference in occupancy, so that
#w2 = w1 (1 – R), with R > 0 representing a decline and
#R < 0 an increase

#For a given R, power increases as the number of sampling sites increases (Fig. 1a). 
#Power also increases with the number of replicate surveys

#The larger the initial occupancy probability w1, the larger the power to detect a
#given proportional difference R, as this translates into a larger absolute occupancy difference

#MacKenzie (2005) and Field, Tyre & Possingham(2005),

###########################################################################################################

#assume the occupancy estimates are normally distributed
#then the difference in occupancy estimates between 2 time points is also normally distributed
#with mean as the difference in occupancies, and variances as the sum of the variances
#effect size is the difference divided by square root of their variance sum

#power analysis of its difference from zero

#or do it on the logistic scale

# test for significance	
Wald <- abs(psi1$mean - psi2$mean)/sqrt(psi1$sd^2 + psi2$sd^2)>qnorm(1-0.05/2)
return(Wald)

Wald <- abs(psi1$mean - psi2$mean)/sqrt(psi1$sd^2 + psi2$sd^2)>qnorm(1-0.05/2)
return(Wald)

#read in nation map for a species
pwr.norm.test() 

#power analysis to detect the observed difference in occupancy at 2 time points
calcPowerFormula <- function(psi1,psi2,var1,var2,alpha){
  limL<- (+qnorm(1-alpha/2)*(sqrt(var1+var2))-(psi1-psi2))/sqrt(var1+var2)
  limU<- (-qnorm(1-alpha/2)*(sqrt(var1+var2))-(psi1-psi2))/sqrt(var1+var2)
  G <- 1 - pnorm(limL) + pnorm(limU)
  return(G)
}

calcPowerFormula(0.9,0.1,0.1,0.1,0.05)

#power analysis to detect a set difference - proportional loss of 10%?
calcsetPower <- function(psi1,R,var1,var2,alpha){
  psi2 <- psi1*(1-R)
  limL<- (+qnorm(1-alpha/2)*(sqrt(var1+var2))-(psi1-psi2))/sqrt(var1+var2)
  limU<- (-qnorm(1-alpha/2)*(sqrt(var1+var2))-(psi1-psi2))/sqrt(var1+var2)
  G <- 1 - pnorm(limL) + pnorm(limU)
  return(G)
}

###########################################################################################################

#get some data:
exampleDF <- readRDS("model-outputs/modelSummary_dynamic_Odonata_adult_Sa.rds")

subset(exampleDF,Param=="psi.fs[1]")
#0.98984651 0.02881651 
subset(exampleDF,Param=="psi.fs[25]")
#0.99387734 0.01941881

calcPowerFormula(psi1=0.98984651,psi2=0.90,var1=0.02881651^2,var2=0.01941881^2,alpha=0.05)
calcsetPower(psi1=0.98984651,R=0.1,var1=0.02881651^2,var2=0.01941881^2,alpha=0.05)

#power to detect a given amount of change
pwr.norm.test()

###########################################################################################################






