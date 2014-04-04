
## ----setup, include=FALSE, cache=FALSE, echo=FALSE-----------------------
library(knitr)
opts_chunk$set(results="hide", warning=FALSE, message=FALSE, fig.path="figures/", 
               highlight=TRUE, fig.show="hide", size="small", tidy=FALSE)
Sys.setenv(TEXINPUTS=getwd(),
           BIBINPUTS=getwd(),
           BSTINPUTS=getwd())


## ----Altitude, cache=TRUE------------------------------------------------
# Import altitudinal data
library(raster)
fname <- "http://hsdm.sourceforge.net/altitude.csv"
alt.df <- read.csv(fname,header=TRUE) 
alt.orig <- rasterFromXYZ(alt.df)
plot(alt.orig)

# Center and scale altitudinal data
alt <- scale(alt.orig,center=TRUE,scale=TRUE)
plot(alt)


## ----theta-binomial------------------------------------------------------
# Load hSDM library
library(hSDM)                                
# Target parameters
beta.target <- matrix(c(1,2,-4),ncol=1)
# Matrix of covariates (including the intercept)
X <- cbind(rep(1,ncell(alt)),values(alt),(values(alt))^2)
# Probability of presence as a quadratic function of altitude
logit.theta <- X %*% beta.target
theta <- inv.logit(logit.theta)
# Transform the probability of presence into a raster
theta <- rasterFromXYZ(cbind(coordinates(alt),theta))
# Plot the probability of presence
plot(theta,col=rev(heat.colors(255)))


## ----observations-binomial-----------------------------------------------
# Load dismo library
library(dismo) # For randomPoints() function
# Number of observations
nobs <- 200
# Set seed for repeatability
seed <- 1234
# Sample the observations in the landscape
set.seed(seed)
obs <- randomPoints(alt,nobs)
# Extract altitude data for observations
alt.obs <- extract(alt,obs)
# Compute theta for these observations
X.obs <- cbind(rep(1,nobs),alt.obs,alt.obs^2)
logit.theta.obs <- X.obs %*% beta.target
theta.obs <- inv.logit(logit.theta.obs)
# Simulate observations
trials <- rep(1,nobs)
set.seed(seed)
Y <- rbinom(nobs,trials,theta.obs)
# Group explicative and response variables in a data-frame
data.obs.df <- data.frame(Y,trials=trials,alt=alt.obs)
# Transform observations in a spatial object
library(sp)
data.obs <- SpatialPointsDataFrame(coords=coordinates(obs),data=data.obs.df)
# Plot observations
plot(alt.orig)
points(data.obs[data.obs$Y==1,],pch=16)
points(data.obs[data.obs$Y==0,],pch=1)


## ------------------------------------------------------------------------
data.obs$alt2 <- (data.obs$alt)^2


## ------------------------------------------------------------------------
data.pred <- data.frame(alt=values(alt),alt2=(values(alt))^2)


## ------------------------------------------------------------------------
mod.hSDM.binomial <- hSDM.binomial(presences=data.obs$Y,
                                   trials=data.obs$trials,
                                   suitability=~alt+alt2,
                                   data=data.obs,
                                   suitability.pred=data.pred,
                                   burnin=1000, mcmc=1000, thin=1,
                                   beta.start=0,
                                   mubeta=0, Vbeta=1.0E6,
                                   seed=1234, verbose=1, save.p=1)


## ----results="markup"----------------------------------------------------
summary(mod.hSDM.binomial$mcmc)


## ----results="markup"----------------------------------------------------
#== glm results for comparison
mod.glm <- glm(cbind(Y,trials-Y)~alt+alt2,family="binomial",data=data.obs)
summary(mod.glm)


## ----mcmc-binomial-------------------------------------------------------
plot(mod.hSDM.binomial$mcmc)


## ----results="markup"----------------------------------------------------
str(mod.hSDM.binomial$prob.p.latent)
summary(mod.hSDM.binomial$prob.p.latent)


## ----predictions-binomial------------------------------------------------
# Create a raster for predictions
theta.pred.mean <- raster(theta)
# Create rasters for uncertainty
theta.pred.2.5 <- theta.pred.97.5 <- raster(theta)
# Attribute predicted values to raster cells
theta.pred.mean[] <- apply(mod.hSDM.binomial$prob.p.pred,2,mean)
theta.pred.2.5[] <- apply(mod.hSDM.binomial$prob.p.pred,2,quantile,0.025)
theta.pred.97.5[] <- apply(mod.hSDM.binomial$prob.p.pred,2,quantile,0.975)
# Plot the predicted probability of presence and uncertainty
plot(theta.pred.mean,col=rev(heat.colors(255)))
plot(theta.pred.2.5,col=rev(heat.colors(255)))
plot(theta.pred.97.5,col=rev(heat.colors(255)))


## ----pred-obs-binomial---------------------------------------------------
# Comparing predictions to initial values
plot(theta[],theta.pred.mean[],cex.lab=1.4)
abline(a=0,b=1,col="red",lwd=2)


## ----siteocc, cache=TRUE-------------------------------------------------
# Import altitudinal data
library(raster)
fname <- "http://hsdm.sourceforge.net/altitude.csv"
alt.df <- read.csv(fname,header=TRUE) 
alt.orig <- rasterFromXYZ(alt.df)
plot(alt.orig)

# Center and scale altitudinal data
alt <- scale(alt.orig,center=TRUE,scale=TRUE)
plot(alt)


