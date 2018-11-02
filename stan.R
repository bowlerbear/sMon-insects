#http://mbjoseph.github.io/2014/11/14/stan_occ.html

#http://mc-stan.org/users/documentation/case-studies/dorazio-royle-occupancy.html

#https://groups.google.com/forum/#!searchin/stan-users/simple$20occupancy$20model$20stan%7Csort:date/stan-users/EpX-wzOj7L4/tJWQ68va9kkJ

#https://groups.google.com/forum/#!searchin/stan-users/occupancy$20model%7Csort:date/stan-users/ZuLmNl4Q9xc/u7ySLj_6OQMJ


#https://www.r-bloggers.com/bayesian-regression-with-stan-part-1-normal-regression/
#http://thestatsgeek.com/2015/06/08/my-first-foray-with-stan/

#dealing with missing values
https://discourse.mc-stan.org/t/how-to-treat-missing-data-na/4498/5


#write model
/**
  * Occupancy model; Single season, no covariates
*
  * Translated from model by Marc Kery, December 2010:
  * http://www.fisheriesstockassessment.com/TikiWiki/tiki-index.php?page=BUGS+Occupancy
  */
    
    data {
      int<lower=0> R;
      int<lower=0> T;
      int<lower=0,upper=1> y[R,T];
    }
  parameters {
    real<lower=0,upper=1> psi1;
    real<lower=0,upper=1> p;
  }
  model {
    // local variables to avoid recomputing log(psi1) and log(1 - psi1)
    real log_psi1;
    real log1m_psi1;
    log_psi1 <- log(psi1);
    log1m_psi1 <- log1m(psi1);
    
    // priors
    psi1 ~ uniform(0,1);
    p ~ uniform(0,1);
    
    // likelihood
    for (r in 1:R) {
      if (sum(y[r]) > 0)
        increment_log_prob(log_psi1 + bernoulli_log(y[r],p));
      else
        increment_log_prob(log_sum_exp(log_psi1 + bernoulli_log(y[r],p),
                                       log1m_psi1));
    }
  }

#simulate data
R <- 250;
T <- 10;
p <- 0.55;
psi <- 0.4;

z <- rbinom(R,1,psi);

y <- matrix(NA,R,T);
for (r in 1:R)
  y[r,] <- rbinom(T,1,z[r] * p);

# to fit in Stan, use this:
# library('rstan')
# fit <- stan('occupancy.stan', data=c("R","T","y"));


#https://groups.google.com/forum/#!searchin/stan-users/occupancy$20model%7Csort:date/stan-users/EpX-wzOj7L4/v6Dgg-QkX9MJ

library(rstan)
library(MASS)

set_cppo("fast")  # for best running speed
# Set values for simulating data
p <- 0.7   #Values for p at intercept (mean of covariate values)
psi <- 0.5 #Values for psi at intercept (mean of covariate value(s))
n.site <- 240 #no. sites
n.visit <- 2 #no. visits (will try n.visit<-c(5,10,20,40) for n.site = 30 and p = 0.1 later)
B.x1 <- 1   #Slope of linear relationship btwn psi and covariate 1
B.x2 <- 0  #Slope of linear relationship btwn psi and covariate 2
A.x1 <- 0   #Slope of linear relationship btwn p and covariate 1
cov.x <- 0    #Correlation between x1 and x2

# Functions for simulating data
expit <-
  function(x){
    exp(x)/(1+exp(x))
  }
logit <-
  function(x){
    log(x/(1-x))
  }

Dat.sim <- function(p0,psi0,n.sts,n.vsts,b.x1,b.x2,a.x1,covx) {
  Y<-matrix(NA,nrow=n.sts,ncol=n.vsts)
  mns<-c(0,0)
  Sgma<-as.matrix(rbind(c(1,covx),c(covx,1)))
  x<-mvrnorm(n.sts,mns,Sgma)
  psi<-expit(logit(psi0)+b.x1*x[,1]+b.x2*x[,2])
  p<-expit(logit(p0)+a.x1*x[,1])
  for(ii in 1:nrow(Y)) {
    
    z<-rbinom(1,1,psi[ii])
    Y[ii,]<-rbinom(n.vsts,1,z*p[ii])
  }
  Dat<-list(Y=Y,psi=psi,p=p,x1=x[,1],x2=x[,2])
  return(Dat)
}

# Simulate data
dat <- Dat.sim(p,psi,n.site,n.visit,B.x1,B.x2,A.x1,cov.x)

R <- nrow(dat$Y)
J <- ncol(dat$Y)
Data <- list(Y=dat$Y,R=R,J=J,x1=dat$x1,x2=dat$x2)

# Set initial values for HMC sampling
inits <- function()
  list(psi_mean=runif(1),p_mean=runif(1),B_x1=rnorm(1),B_x2=rnorm(1),A_x1=rnorm(1))

#Define which parameters to save.
parameters <- c("psi_mean","p_mean","B0","B_x1","B_x2","A0","A_x1")

# Model code
mod_occupancy_simple <- '
data {
int<lower=0> R; // number of sites
int<lower=0> J; // number of visits
int<lower=0,upper=1> Y[R,J]; // detection-nondetection data
vector[R] x1; // 1st covariate
vector[R] x2; // 2nd covariate
}
parameters {
real<lower=0,upper=1> psi_mean;
real<lower=0,upper=1> p_mean;
real B_x1;   
real B_x2;   
real A_x1;   
}
transformed parameters {
real B0; // intercept for linear logistic model for occupancy
real A0; // intercept for linear logistic model for detection
vector[R] psi;
matrix[R,J] p;

B0 <- logit(psi_mean);
A0 <- logit(p_mean);
for (r in 1:R) {
psi[r] <- inv_logit(B0 + B_x1*x1[r] + B_x2*x2[r]);
for (j in 1:J) p[r,j] <- inv_logit(A0 + A_x1*x1[r]);
}
}
model {
// local variables to avoid recomputing log(psi) and log(1-psi)
vector[R] log_psi;
vector[R] log1m_psi;

for (r in 1:R) {
log_psi[r] <- log(psi[r]);
log1m_psi[r] <- log1m(psi[r]);
}
// priors

psi_mean ~ uniform(0,1);
p_mean ~ uniform(0,1);
B_x1 ~ normal(0,10);
B_x2 ~ normal(0,10);
A_x1 ~ normal(0,10);
// likelihood
for (r in 1:R) {
if (sum(Y[r]) > 0)
increment_log_prob(log_psi[r] + bernoulli_log(Y[r],p[r]));
else
increment_log_prob(log_sum_exp(log_psi[r] + bernoulli_log(Y[r],p[r]),log1m_psi[r]));
}
}
'
# Settings for HMC chains
nc <- 6
nb <- 1000
ni <- 2000
nt <- 10
#Run model
out <- stan(model_code = mod_occupancy_simple, data = Data, pars = parameters, init = inits,
            iter = ni, chains = nc, warmup = nb, thin = nt)
Results:
  
  
  TRANSLATING MODEL 'mod_occupancy_simple' FROM Stan CODE TO C++ CODE NOW.
COMPILING THE C++ CODE FOR MODEL 'mod_occupancy_simple' NOW.
cygwin warning:
  MS-DOS style path detected: C:/PROGRA~1/R/R-30~1.1/etc/x64/Makeconf
Preferred POSIX equivalent is: /cygdrive/c/PROGRA~1/R/R-30~1.1/etc/x64/Makeconf
CYGWIN environment variable option "nodosfilewarning" turns off this warning.
Consult the user's guide for more details about POSIX paths:
http://cygwin.com/cygwin-ug-net/using.html#using-pathnames
C:/Program Files/R/R-3.0.1/library/rstan/include//stansrc/stan/agrad/
rev/var_stack.hpp:49:17: warning: 'void stan::agrad::free_memory()' defined but not used [-Wunused-function]
SAMPLING FOR MODEL 'mod_occupancy_simple' NOW (CHAIN 1).
Iteration: 2000 / 2000 [100%]  (Sampling)
Elapsed Time: 2.604 seconds (Warm-up)
2.148 seconds (Sampling)
4.752 seconds (Total)

SAMPLING FOR MODEL 'mod_occupancy_simple' NOW (CHAIN 2).
Iteration: 2000 / 2000 [100%]  (Sampling)
Elapsed Time: 2.554 seconds (Warm-up)
2.452 seconds (Sampling)
5.006 seconds (Total)

SAMPLING FOR MODEL 'mod_occupancy_simple' NOW (CHAIN 3).
Iteration: 2000 / 2000 [100%]  (Sampling)
Elapsed Time: 2.764 seconds (Warm-up)
2.238 seconds (Sampling)
5.002 seconds (Total)

SAMPLING FOR MODEL 'mod_occupancy_simple' NOW (CHAIN 4).
Iteration: 2000 / 2000 [100%]  (Sampling)
Elapsed Time: 2.553 seconds (Warm-up)
2.429 seconds (Sampling)
4.982 seconds (Total)

SAMPLING FOR MODEL 'mod_occupancy_simple' NOW (CHAIN 5).
Iteration: 2000 / 2000 [100%]  (Sampling)
Elapsed Time: 2.585 seconds (Warm-up)
2.431 seconds (Sampling)
5.016 seconds (Total)

SAMPLING FOR MODEL 'mod_occupancy_simple' NOW (CHAIN 6).
Iteration: 2000 / 2000 [100%]  (Sampling)
Elapsed Time: 2.606 seconds (Warm-up)
2.291 seconds (Sampling)
4.897 seconds (Total)

> out
Inference for Stan model: mod_occupancy_simple.
6 chains, each with iter=2000; warmup=1000; thin=10; 
post-warmup draws per chain=100, total post-warmup draws=600.

mean se_mean  sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
psi_mean    0.5     0.0 0.0    0.4    0.5    0.5    0.5    0.6   600    1
p_mean      0.7     0.0 0.0    0.6    0.7    0.7    0.7    0.8   600    1
B0          0.0     0.0 0.2   -0.3   -0.1    0.0    0.1    0.3   600    1
B_x1        1.0     0.0 0.2    0.6    0.8    1.0    1.1    1.4   600    1
B_x2       -0.1     0.0 0.2   -0.4   -0.2   -0.1    0.0    0.2   600    1
A0          1.0     0.0 0.2    0.6    0.8    0.9    1.1    1.4   600    1
A_x1       -0.1     0.0 0.2   -0.4   -0.2   -0.1    0.0    0.3   592    1
lp__     -262.5     0.1 1.5 -266.3 -263.3 -262.2 -261.4 -260.6   600    1

Samples were drawn using NUTS(diag_e) at Sun Jan 05 12:08:04 2014.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).