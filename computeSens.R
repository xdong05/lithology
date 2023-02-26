rm(list = ls())

library(rjags)
library(runjags)
library(coda)

data30ndvi = read.csv("ndvi.csv",header=TRUE)# input data: NDVI data 
dfff   = read.csv("climatewaterdef.csv",header=TRUE)# input data: CWD data
sss    = read.csv("covariates.csv",header=TRUE) #matrix with all predictors used in StageII model
nSites = nrow(sss)
sss$sCwd   = NA
sss$sCwdSD = NA

# model intials 
for(k in 1:3){
  assign (make.names(k), list(alpha = rnorm(1,mean=0.001, sd=0.005),beta =rnorm(1,mean=0.001, sd=0.005))) 
}
myInitials <-list(X1,X2,X3)

for (k in 1:nSites){#loop through each grid to calculate sensitivity 
  print(k)
  
  dummy  = intersect(which(data30ndvi$lat==sss$lat[k]), which(data30ndvi$lon==sss$lon[k]));
  dummy2 = intersect(which(dfff$lat==sss$lat[k]), which(dfff$lon==sss$lon[k]));
  df2 = data30ndvi[dummy,]
  df3 = df2[order(df2$year),]
  ys  = df3$ndvi_annual   #NDVI time series for a given grid
  xs  = dfff$cwd[dummy2]  #CWD time series for a given grid
  
  nn = intersect(which(!is.na(ys)),which(!is.na(xs)))
  allYears = seq(1,32,by=1)
  if (length(nn) < 8 | length(unique(ys[nn]))< 4 |length(unique(xs[nn])) < 4 ){
    #Remove grids whose time series is too short to infer sensitivity reliably 
    sss$sCwd[k]    <- NA
    sss$sCwdSD[k]  <- NA
    
  }else{
    #detrend CWD and NDVI time series before computing sensitivity
    y_fit <- lm(ys ~ allYears)
    y_detrend = ys - as.numeric(y_fit$fitted.values)  
    x_fit <- lm(xs ~ allYears)
    x_detrend = xs - as.numeric(x_fit$fitted.values)  
    
    jags.model.sp ="
    model{
    for (i in 1:N){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta*x[i]
    }
    
     beta  ~ dnorm(0,1)
     alpha ~ dnorm(0,1)
     tau   ~ dgamma(0.1,0.1)
     
    }
    " 
    jd<-list(y=y_detrend, x=x_detrend, N=nSites)
    jags.out <- run.jags(jags.model.sp,
                         data=jd,
                         adapt = 5000,
                         burnin = 3000,
                         sample = 3000,
                         n.chains = 3,
                         thin = 20,inits = myInitials,
                         monitor=c('beta'), 
                         modules='glm')
    
    outs = as.mcmc.list(jags.out)
    r1 = as.matrix(outs)
    bs = as.data.frame(r1)
    sss$sCwd[k]   = mean(bs[,1])
    sss$sCwdSD[k] = sd(bs[,1])
  }
  
}

write.csv(sss,"computedSens.csv")



