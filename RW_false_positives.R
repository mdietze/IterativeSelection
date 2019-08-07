## Iterative model selection: Random walk vs. DLM
## Goal: determine the rate of false positives when null model is true (RW)
##   subject to early stopping
devtools::install_github("EcoForecast/EF_Activities")
library(ecoforecastR)

## Settings & true parameters
nboot = 10 ## number of pseudodata simulations
nt = 100  ## simulation length
nf = 4    ## number of steps ahead to forecast & evaluate
dt = 4    ## number of steps between refits
nfit = floor(nt/dt)
tau = 1   ## process error; conclusion should be invariant to tau because can always rescale y-axis
sigma = 0.1 ## observation error
models <- list()
models[[1]] = ""
models[[2]] = "~ 1 + X"
fp <- list()

## model selection metrics to store
DICfp = matrix(NA,nboot,nfit)
DICfp.t <- matrix(NA,length(models),nfit)
PLfp <- RMSEfp <- array(NA,c(length(models),nf,nfit,nboot)) ## model, step ahead, time step, replicate

for(b in seq_len(nboot)){
  
  ## Simulate data
  Yb.TRUE = cumsum(rnorm(nt,0,tau))
  Yb = rnorm(nt,Yb.TRUE,sigma)
  
  ## Fit models
  for(i in seq_along(models)){
    fp[[i]] <- iterative_fit_dlm(model=list(fixed=models[[i]],obs="Y",n.iter=10000),
                                  data=data.frame(Y=Yb),"refit",dt=dt,nf=nf)
  }
  
  ## Calculate DICs along the way
  for(i in seq_along(fp)){
    DICfp.t[i,] <- apply(fp[[i]]$dic,1,sum)
  }
  DICfp[b,] = DICfp.t[2,] - DICfp.t[1,] ## only need to store the difference (neg = false positive)
  
  ## Calculate RMSE & Predictive Loss
  for(m in seq_along(fp)){
    for(i in seq_len(nf)){
      for(t in 1:length(tpred)){
        RMSEfp[m,i,t,b] <- sqrt(mean((fp[[m]]$xp[1:t,4,i]-Yb[tpred[1:t]+i])^2,na.rm = TRUE))
      }
      PLfp[m,i,,b] <- sqrt(SFO[[m]]$xp[,5,i]+RMSEfp[m,i,,b]^2)
    }
  }
  
  ## save periodically
  if(b %% 3 == 0){
    print(paste("saving",b))
    save(DICfp,RMSEfp,PLfp,file="RWfp.Rdata")
    print("done")
  }
}

## final save
save(DICfp,RMSEfp,PLfp,file="RWfp.Rdata")
## note: will load up this output into main Rmd for visualization
hist(DICfp[,nfit])
for(i in 1:nf){  ## for all these, < 1 is a false positive
  hist(PLfp[2,i,nfit,]/PLfp[1,i,nfit,],main=paste("PL:",i,"step ahead"))
  hist(RMSEfp[2,i,nfit,]/RMSEfp[1,i,nfit,],main=paste("RMSE:",i,"step ahead"))
}
