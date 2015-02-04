
library(doMC)
registerDoMC(3)
nchains=3

library(hSDM)
library(raster)
library(sp)

data("cfr.env")
cfr.env$cell=1:nrow(cfr.env)

#= Sample the observation sites in the landscape
nsite <- 250
seed=1234
set.seed(seed)

covar=c("min07","smdwin","fert3")

cfr.env[,covar]=scale(cfr.env[,covar])

pts=cfr.env[sample(1:nrow(cfr.env),nsite,replace = F),]
X=as.matrix(pts[,covar])

## number of visits
set.seed(seed)
visits<- rpois(nsite,3)
visits[visits==0] <- 1

#= Ecological process (suitability)
set.seed(seed)
beta.target <- c(-1,1,-1)
logit.theta <- X %*% beta.target
theta <- inv.logit(logit.theta)
set.seed(seed)
Y <- rbinom(nsite,visits,theta)

#= Data-sets
data.obs <- data.frame(presences=Y,trials=visits,X,cell=pts$cell)

## metadata
meta=list(species="simulated species")

#==================================
#== Site-occupancy model
library(coda)

mod1=foreach(ch=1:nchains) %dopar% {
  hSDM.ZIB(presences=data.obs$presences,
                  trials=data.obs$trials,
                  suitability=as.formula(paste0("~",paste0(covar,collapse="+"))),
                  observability=~1,
                  data=data.obs,
                  suitability.pred=na.omit(cfr.env),
                  burnin=1000, mcmc=100, thin=5,
                  beta.start=0,gamma.start=0,
                  mubeta=0, Vbeta=1.0E6,
                  seed=1234, verbose=1,
                  save.p=0,
                  meta=meta)
}
## summarize and write to disk

hSDM.ncWriteOutput(mod1,file="mod1.nc",overwrite=T,meta=list(modelname="Full"),autocor=F)

## reduced model (first two covariates)
mod2=foreach(ch=1:nchains) %dopar% {
  hSDM.ZIB(presences=data.obs$presences,
                   trials=data.obs$trials,
                   suitability=as.formula(paste0("~",paste0(covar[1:2],collapse="+"))),
                   observability=~1,
                   data=data.obs,
                   suitability.pred=na.omit(cfr.env),
                   burnin=1000, mcmc=100, thin=5,
                   beta.start=0,gamma.start=0,
                   mubeta=0, Vbeta=1.0E6,
                   seed=1234, verbose=1,
                   save.p=0)
}

## summarize and write to disk
hSDM.ncWriteOutput(mod2,file="mod2.nc",overwrite=T,meta=list(modelname="Reduced"),autocor=F)


## show contents of file
nc=nc_open("mod1.nc")
print(nc)

## get list of all nc files available for processing
fresults=list.files(".",pattern=".nc$",full=T,recursive=T)

## extract model comparison for all files
eval=hSDM.ncExtract(files=fresults,what="eval")
eval

## Extract parameter summaries (and convergence metrics) for all files
coef=hSDM.ncExtract(files=fresults,what="coef")
coef

## Clean up
file.remove(fresults)
