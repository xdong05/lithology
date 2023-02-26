#Bayesian hierarchical model to analyze the effect of bedrock properties; 

rm(list = ls())

library(rjags)
library(runjags)
library(coda)
library(MASS)
library(parallel)

sss = read.csv("computedSens.csv",header=TRUE)
D   = read.csv("distance.csv",header=F)
D   = as.matrix(D)

# --- Standardize and transform variables ----
sss$porosity     = (sss$porosity-mean(sss$porosity,na.rm=TRUE))/sd(sss$porosity,na.rm=TRUE)
sss$permeability = (sss$permeability-mean(sss$permeability,na.rm=TRUE))/sd(sss$permeability,na.rm=TRUE)
sss$zSoil = (sss$zSoil-mean(sss$zSoil,na.rm=TRUE))/sd(sss$zSoil,na.rm=TRUE)
sss$zReg  = (sss$zReg-mean(sss$zReg,na.rm=TRUE))/sd(sss$zReg,na.rm=TRUE)
sss$sCwd  = log(abs(sss$sCwd))

nSites = nrow(sss)
# SD of ecosystem sensitivity calculated from Stage I analysis; used as "data"
# in the model here;it represents the uncertainty associated with the inferred ecosystem 
# sensitivity, and is used as "weights" in the current model
sd_p   = sss$sCwdSD 
sdM_p  = matrix(ncol = nSites,nrow = nSites)
for (i in 1:nSites){
  for (j in 1:nSites){
    sdM_p[i,j] = sd_p[i]*sd_p[j]
  }
} 


jags.model.sp ="

model{
    # Likielihood of ecosystem sensitivity calculated from first stage analysis
    sCWD ~ dmnorm.vcov(Mu_p, Sigma)
    
    # Replicated data for evaluating model fit and calculating bayesian R2
    sCWD_rep ~ dmnorm.vcov(Mu_p, Sigma)
    
    
    for (i in 1:N){ #loop through N grids
    # Mean model to infer the effect of regolith porosity ('pors'),permeability ('perm'),
    # soil thickness ('zSoil') and regolith thickness ('upReg'),with covariates varying with 
    # regions (four regions in total);
    Mu_p[i] <- b0+alpha[rockType[i]]+b1[region[i]]*pors[i]+b2[region[i]]*perm[i]+b3[region[i]]*upReg[i]+b4[region[i]]*zSoil[i]
    
    # calculate the var-cov matrix to account for the effect of spatial autocorrelation
    Sigma[i,i] <- etasq * sdM_p[i,i]
    for (j in (i+1):N) {
          Sigma[i,j] <- (etasq * sdM_p[i,j]) * exp(-rhosq * pow(D[i,j],2)) 
          Sigma[j,i] <- Sigma[i,j]
    }
     
    
    # imputating missing values in predictors 
    pors[i]  ~ dnorm(0,1)
    perm[i]  ~ dnorm(0,1)
    upReg[i] ~ dnorm(0,1)
    zSoil[i] ~ dnorm(0,1)
    rockType[i] ~ dcat(rPriorProb)
    region[i]   ~ dcat(gPriorProb)
    }
    
    # 15rock types in RegionI;14 rock types in RegionII;11 rock types in Region III;and 
    # 11rock types in Region IV (total = 51)
    for (k in 1:51) { 
    rPriorProb[k] <- 1/51
    }
    
    #four regions 
    for (k in 1:4){ 
    gPriorProb[k] <- 1/4
    }


   for (k in 1:15) { #15 types of rocks in Region I
   alpha[k] ~ dnorm(0,1)
   alpha.star[k] <- alpha[k]-ave.alpha1 #post sweeping 
   }
   ave.alpha1<-mean(alpha[1:15])
   
   for (k in 16:29) {  #14 types of rocks in Region II
   alpha[k] ~ dnorm(0,1)
   alpha.star[k] <- alpha[k]-ave.alpha2 #post sweeping 
   }
   ave.alpha2<-mean(alpha[16:29])
   
   for (k in 30:40) { #11 types of rocks in Region III
   alpha[k] ~ dnorm(0,1)
   alpha.star[k] <- alpha[k]-ave.alpha3 #post sweeping 
   }
   ave.alpha3 <- mean(alpha[30:40])
   
   for (k in 41:51) { #11 types of rocks in Region IV
   alpha[k] ~ dnorm(0,1)
   alpha.star[k] <- alpha[k]-ave.alpha4 #post sweeping 
   }
   ave.alpha4 <- mean(alpha[41:51])
   
   ave.ave.alpha <- (ave.alpha1+ave.alpha2+ave.alpha3+ave.alpha4)/4
   
   # priors
   tau_p2 ~ dgamma(0.1,0.1)
   b0     ~ dnorm(0,1)
   b0.star <- b0 + ave.ave.alpha #post sweeping for the intercept
   
   rhosq ~ dexp(0.5)
   etasq ~ dexp(2)
   
   # hierarhical priors pooling information across regions
   for(k in 1:4){ 
   
   b1[k] ~ dnorm(hyper_mu1,hyper_precis1)
   b2[k] ~ dnorm(hyper_mu2,hyper_precis2)
   b3[k] ~ dnorm(hyper_mu3,hyper_precis3)
   b4[k] ~ dnorm(hyper_mu4,hyper_precis4)
   
   }
   
   hyper_mu1 ~ dnorm(0,1)
   hyper_mu2 ~ dnorm(0,1)
   hyper_mu3 ~ dnorm(0,1)
   hyper_mu4 ~ dnorm(0,1)
   
   hyper_precis1 ~ dgamma(0.1,0.1)
   hyper_precis2 ~ dgamma(0.1,0.1)
   hyper_precis3 ~ dgamma(0.1,0.1)
   hyper_precis4 ~ dgamma(0.1,0.1)
}
"

jd<-list(N = nSites,D = D,sCWD = sss$sCwd,zSoil = sss$zReg,upReg = sss$zReg,
         pors = sss$porosity,perm = sss$permeability,sdM_p = sdM_p,
         region = sss$region,rockType = sss$rockType)


# sensible initial values for part of the parameters to get the model started
for(k in 1:20){
  assign (make.names(k), list(b0 =rnorm(1, mean=0.001, sd=0.005),b1 =rnorm(4, mean=0.001, sd=0.005),
                      b2 =rnorm(4, mean=0.001, sd=0.005),b3 =rnorm(4, mean=0.001, sd=0.005),
                      b4 =rnorm(4, mean=0.001, sd=0.005),rhosq = 0.01,etasq = 0.1)) 
}
myInitials <-list(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20)

# fit the model
jags.out <- run.jags(jags.model.sp,
                     data   = jd,
                     adapt  = 5000,
                     burnin = 3000,
                     sample = 3000,
                     n.chains = 20,
                     thin = 20,inits = myInitials,
                     monitor=c('b1','b2','b3','b4','b0.star','alpha.star','rhosq','etasq','sCWD_rep'), 
                     modules='glm',method="parallel")


# export the all monitored parameters of interest 
outs = as.mcmc.list(jags.out)
r1   = as.matrix(outs)
bs   = as.data.frame(r1)

write.csv(bs,"results.csv")


