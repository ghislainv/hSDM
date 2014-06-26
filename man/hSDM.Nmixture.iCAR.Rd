\name{hSDM.Nmixture.iCAR}
\alias{hSDM.Nmixture.iCAR}

\title{N-mixture model with CAR process}

\description{The \code{hSDM.Nmixture.iCAR} function can be used to model
  species distribution including different processes in a hierarchical
  Bayesian framework: a \eqn{\mathcal{P}oisson}{Poisson} suitability
  process (refering to environmental suitability explaining abundance)
  which takes into account the spatial dependence of the observations,
  and a \eqn{\mathcal{B}inomial}{Binomial} observability process
  (refering to various ecological and methodological issues explaining
  the species detection). The \code{hSDM.Nmixture.iCAR} function calls a
  Gibbs sampler written in C code which uses an adaptive Metropolis
  algorithm to estimate the conditional posterior distribution of
  hierarchical model's parameters.}

\usage{hSDM.Nmixture.iCAR(# Observations
                     counts, observability, site, data.observability,
                     # Habitat
                     suitability, data.suitability,
                     # Spatial structure
                     spatial.entity,
                     n.neighbors, neighbors,
                     # Predictions
                     suitability.pred = NULL, spatial.entity.pred = NULL,
                     # Chains
                     burnin = 5000, mcmc = 10000, thin = 10,
                     # Starting values
                     beta.start,
                     gamma.start,
                     Vrho.start,
                     # Priors
                     mubeta = 0, Vbeta = 1.0E6,
                     mugamma = 0, Vgamma = 1.0E6,
                     priorVrho = "1/Gamma",
                     shape = 0.5, rate = 0.0005,
                     Vrho.max = 1000,
                     # Various
                     seed = 1234, verbose = 1,
                     save.rho = 0, save.p = 0, save.N = 0)}

\arguments{
  
  \item{counts}{A vector indicating the count (or abundance) for each
    observation.}

  \item{observability}{A one-sided formula of the form
    \eqn{\sim w_1+...+w_q}{'~w_1+...+w_q'} with \eqn{q} terms specifying the
    explicative variables for the observability process.}

  \item{site}{A vector indicating the site identifier (from one to the
    total number of sites) for each observation. Several observations
    can occur at one site. A site can be a raster cell for example.}
  
  \item{data.observability}{A data frame containing the model's
    variables for the observability process.}

  \item{suitability}{A one-sided formula of the form
    \eqn{\sim x_1+...+x_p}{'~x_1+...+x_p'} with \eqn{p} terms specifying the
    explicative variables for the suitability process.}

  \item{data.suitability}{A data frame containing the model's variables
    for the suitability process. The number of rows of the data frame
    should be equal to the total number of spatial entities.}

  \item{spatial.entity}{A vector (of length 'nsite') indicating the
    spatial entity identifier for each site. Values must be between 1
    and the total number of spatial entities. Several sites can be found
    in one spatial entity. A spatial entity can be a raster cell for
    example.}

  \item{n.neighbors}{A vector of integers that indicates the number of
    neighbors (adjacent entities) of each spatial
    entity. \code{length(n.neighbors)} indicates the total number of
    spatial entities.}
  
  \item{neighbors}{A vector of integers indicating the neighbors
    (adjacent entities) of each spatial entity. Must be of the form
    c(neighbors of entity 1, neighbors of entity 2, ... , neighbors of
    the last entity). Length of the \code{neighbors} vector should be
    equal to \code{sum(n.neighbors)}.}

  \item{suitability.pred}{An optional data frame in which to look for
    variables with which to predict. If NULL, the data frame
    \code{data.suitability} for observations is used.}

  \item{spatial.entity.pred}{An optional vector indicating the spatial
    entity identifier (from one to the total number of entities) for
    predictions. If NULL, the vector \code{spatial.entity} for
    observations is used.}

  \item{burnin}{The number of burnin iterations for the sampler.}
    
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to
    \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and
    superior or equal to 100 so that the progress bar can be displayed.}
    
  \item{thin}{The thinning interval used in the simulation. The number
    of mcmc iterations must be divisible by this value.}

  \item{beta.start}{Starting values for \eqn{\beta}{beta} parameters of
    the suitability process. This can either be a scalar or a
    \eqn{p}-length vector.}

  \item{gamma.start}{Starting values for \eqn{\beta}{gamma} parameters
    of the observability process. This can either be a scalar or a
    \eqn{q}-length vector.}
  
  \item{Vrho.start}{Positive scalar indicating the starting value for the
    variance of the spatial random effects.}

  \item{mubeta}{Means of the priors for the \eqn{\beta}{beta} parameters
    of the suitability process. \code{mubeta} must be either a scalar or a
    p-length vector. If \code{mubeta} takes a scalar value, then that value will
    serve as the prior mean for all of the betas. The default value is set
    to 0 for an uninformative prior.}

  \item{Vbeta}{Variances of the Normal priors for the \eqn{\beta}{beta}
    parameters of the suitability process. \code{Vbeta} must be either a
    scalar or a p-length vector. If \code{Vbeta} takes a scalar value,
    then that value will serve as the prior variance for all of the
    betas. The default variance is large and set to 1.0E6 for an
    uninformative flat prior.}

  \item{mugamma}{Means of the Normal priors for the \eqn{\gamma}{gamma}
    parameters of the observability process. \code{mugamma} must be either
    a scalar or a p-length vector. If \code{mugamma} takes a scalar value,
    then that value will serve as the prior mean for all of the
    gammas. The default value is set to 0 for an uninformative prior.}

  \item{Vgamma}{Variances of the Normal priors for the
    \eqn{\gamma}{gamma} parameters of the observability
    process. \code{Vgamma} must be either a scalar or a p-length
    vector. If \code{Vgamma} takes a scalar value, then that value will
    serve as the prior variance for all of the gammas. The default
    variance is large and set to 1.0E6 for an uninformative flat prior.}

  \item{priorVrho}{Type of prior for the variance of the spatial random
    effects. Can be set to a fixed positive scalar, or to an inverse-gamma
    distribution ("1/Gamma") with parameters \code{shape} and \code{rate},
    or to a uniform distribution ("Uniform") on the interval
    [0,\code{Vrho.max}]. Default set to "1/Gamma".}
  
  \item{shape}{The shape parameter for the Gamma prior on the precision
    of the spatial random effects. Default value is \code{shape=0.05} for
    uninformative prior.}
  
  \item{rate}{The rate (1/scale) parameter for the Gamma prior on the
    precision of the spatial random effects. Default value is
    \code{rate=0.0005} for uninformative prior.}

  \item{Vrho.max}{Upper bound for the uniform prior of the spatial
    random effect variance. Default set to 1000.}
 
  \item{seed}{The seed for the random number generator. Default set to
  1234.}
  
  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the
    Gibbs sampler.}

  \item{save.rho}{A switch (0,1) which determines whether or not the
    sampled values for rhos are saved. Default is 0: the posterior mean
    is computed and returned in the \code{rho.pred} vector. Be careful,
    setting \code{save.rho} to 1 might require a large amount of
    memory.}

  \item{save.p}{A switch (0,1) which determines whether or not the
    sampled values for predictions are saved. Default is 0: the
    posterior mean is computed and returned in the \code{lambda.pred}
    vector. Be careful, setting \code{save.p} to 1 might require a large
    amount of memory.}

  \item{save.N}{A switch (0,1) which determines whether or not the
    sampled values for the latent count variable N for each observed
    cells are saved. Default is 0: the mean (rounded to the closest
    integer) is computed and returned in the \code{N.pred} vector. Be
    careful, setting \code{save.N} to 1 might require a large amount of
    memory.}

}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_{it} P(y_{it},N_i|...))}{%
      D=-2log(prod_it P(y_it,N_i|...))}, is also
    provided.}

  \item{rho.pred}{If \code{save.rho} is set to 0 (default),
    \code{rho.pred} is the predictive posterior mean of the spatial
    random effect associated to each spatial entity. If \code{save.rho} is
    set to 1, \code{rho.pred} is an \code{mcmc} object with sampled
    values for each spatial random effect associated to each spatial
    entity.}

  \item{lambda.pred}{If \code{save.p} is set to 0 (default),
    \code{lambda.pred} is the predictive posterior mean of the abundance
    associated to the suitability process for each prediction. If
    \code{save.p} is set to 1, \code{lambda.pred} is an \code{mcmc}
    object with sampled values of the abundance associated to the
    suitability process for each prediction.}

  \item{N.pred}{If \code{save.N} is set to 0 (default), \code{N.pred} is
    the posterior mean (rounded to the closest integer) of the latent
    count variable N for each observed cell. If \code{save.N} is set to 1,
    \code{N.pred} is an \code{mcmc} object with sampled values of the
    latent count variable N for each observed cell.}
  
  \item{lambda.latent}{Predictive posterior mean of the abundance
    associated to the suitability process for each observation.}
  
  \item{delta.latent}{Predictive posterior mean of the probability
  associated to the observability process for each observation.}
  
}

\details{

  The model integrates two processes, an ecological process associated
  to the abundance of the species due to habitat suitability and an
  observation process that takes into account the fact that the
  probability of detection of the species is inferior to one. The
  ecological process includes an intrinsic conditional autoregressive
  model (iCAR) model for spatial autocorrelation between observations,
  assuming that the abundance of the species at one site depends on the
  abundance of the species on neighboring sites.

  \bold{Ecological process:}
  \deqn{N_i \sim \mathcal{P}oisson(\lambda_i)}{N_i ~ Poisson(lambda_i)}
  \deqn{log(\lambda_i) = X_i \beta + \rho_i}{log(lambda_i) = X_i *
    beta + rho_i}
  \eqn{\rho_i}{rho_i}: spatial random effect
  
  \bold{Spatial autocorrelation:}
  
  An intrinsic conditional autoregressive model (iCAR) is assumed:
  \deqn{\rho_i \sim \mathcal{N}ormal(\mu_i,V_{\rho} / n_i)}{rho_i ~
    Normal(mu_i,V_rho / n_i)}
  \eqn{\mu_i}{mu_i}: mean of \eqn{\rho_{i'}}{rho_i'} in the
  neighborhood of \eqn{i}.
  
  \eqn{V_{\rho}}{V_rho}: variance of the spatial random effects.
  
  \eqn{n_i}{n_i}: number of neighbors for spatial entity \eqn{i}.
  
  \bold{Observation process:}
  \deqn{y_{it} \sim \mathcal{B}inomial(N_i, \delta_{it})}{y_it ~ Binomial(N_i, delta_it)}
  \deqn{logit(\delta_{it}) = W_{it} \gamma}{logit(delta_it) = W_it * gamma} 

}

\references{

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.
  
  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
    Applications}, 16, 33-50.

  Royle, J. A. (2004) N-mixture models for estimating population size from
  spatially replicated counts. \emph{Biometrics}, 60, 108-115.

}

\author{
  Ghislain Vieilledent \email{ghislain.vieilledent@cirad.fr}
}

\seealso{
  \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
}

\examples{

\dontrun{

#==============================================
# hSDM.Nmixture.iCAR()
# Example with simulated data
#==============================================

#============
#== Preambule
library(hSDM)
library(raster)
library(sp)
library(mvtnorm)

#==================
#== Data simulation

# Set seed for repeatability
seed <- 1234

# Target parameters
beta.target <- matrix(c(0.2,0.5,0.5),ncol=1)
gamma.target <- matrix(c(1),ncol=1) #= This means that theta.2~0.75: 75\% of chance to detect the species
## Uncomment if you want covariates on the observability process
## gamma.target <- matrix(c(1,0.5,0.5),ncol=1)
Vrho.target <- 1 # Spatial Variance

# Landscape
Landscape <- raster(ncol=20,nrow=20,crs='+proj=utm +zone=1')
ncell <- ncell(Landscape)

# Neighbors
neighbors.mat <- adjacent(Landscape, cells=c(1:ncell), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]

# Generate symmetric adjacency matrix, A
A <- matrix(0,ncell,ncell)
index.start <- 1
for (i in 1:ncell) {
    index.end <- index.start+n.neighbors[i]-1
    A[i,adj[c(index.start:index.end)]] <- 1
    index.start <- index.end+1
}

# Spatial effects
d <- 1	# Spatial dependence parameter = 1 for intrinsic CAR
Q <- diag(n.neighbors)-d*A + diag(.0001,ncell) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos
set.seed(seed)
rho <- c(rmvnorm(1,sigma=covrho)) # Spatial Random Effects
rho <- rho-mean(rho) # Centering rhos on zero

# Visited cells
n.visited <- 150 # Compare with 400, 350 and 100 for example
set.seed(seed)
visited.cells <- sort(sample(1:ncell,n.visited,replace=FALSE)) # Draw visited cells at random
notvisited.cells <- c(1:ncell)[-visited.cells]

# Number of observations
nobs <- 300

# Cell vector
set.seed(seed)
cells <- c(visited.cells,sample(visited.cells,nobs-n.visited,replace=TRUE))
coords <- xyFromCell(Landscape,cells) # Get coordinates

# Covariates for "suitability" process
set.seed(seed)
X1 <- rnorm(n=ncell,0,1) # Here, we decided to define one habitat by spatial entity (raster cell)
set.seed(2*seed)
X2 <- rnorm(n=ncell,0,1)
X <- cbind(rep(1,ncell),X1,X2)

# Covariates for "observability" process
W <- cbind(rep(1,nobs))
## Uncomment if you want covariates on the observability process 
## set.seed(3*seed)
## W1 <- rnorm(n=nobs,0,1)
## set.seed(4*seed)
## W2 <- rnorm(n=nobs,0,1)
## W <- cbind(rep(1,nobs),W1,W2)

#== Simulating latent variables

# Suitability
log.lambda <- vector()
log.lambda <- X\%*\%beta.target+rho
lambda <- exp(log.lambda)
set.seed(seed)
N <- rpois(ncell,lambda)

# Observability
logit.theta <- W\%*\%gamma.target
theta <- inv.logit(logit.theta)

#== Simulating response variable
set.seed(seed)
Y <- rbinom(nobs,N[cells],theta)

#== Data-sets
Data.obs <- data.frame(Y,cells)
## Uncomment if you want covariates on the observability process 
## Data.obs <- data.frame(Y,cells,W1,W2)
Data.suit <- data.frame(cells=c(1:ncell),X1,X2)
Data.obs <- SpatialPointsDataFrame(coords=coords,data=Data.obs)
plot(Data.obs)
Data.pred <- Data.suit

#==================================
#== N-mixture model

mod.hSDM.Nmixture.iCAR <- hSDM.Nmixture.iCAR(# Observations
                                   counts=Data.obs$Y, observability=~1,
                                   spatial.entity=Data.obs$cells, data.observability=Data.obs,
                                   # Habitat
                                   suitability=~X1+X2, data.suitability=Data.suit,
                                   # Spatial structure
                                   n.neighbors=n.neighbors, neighbors=adj,
                                   # Predictions
                                   suitability.pred=Data.pred, spatial.entity.pred=Data.pred$cells,
                                   # Chains
                                   burnin=10000, mcmc=10000, thin=10,
                                   # Starting values
                                   beta.start=0,
                                   gamma.start=0,
                                   Vrho.start=10,
                                   # Priors
                                   mubeta=0, Vbeta=1.0E6,
                                   mugamma=0, Vgamma=1.0E6,
                                   priorVrho="1/Gamma",
                                   #priorVrho="Uniform",
                                   #priorVrho=10,
                                   shape=0.5, rate=0.0005,
                                   #Vrho.max=1000,
                                   # Various
                                   seed=1234, verbose=1,
                                   save.rho=1, save.p=0, save.N=1)

#==========
#== Outputs

#= Parameter estimates
summary(mod.hSDM.Nmixture.iCAR$mcmc)

#= MCMC and posteriors
pdf(file="Posteriors_hSDM.Nmixture.iCAR.pdf")
plot(mod.hSDM.Nmixture.iCAR$mcmc)
dev.off()

pdf(file="Posteriors.rho_hSDM.Nmixture.iCAR.pdf")
plot(mod.hSDM.Nmixture.iCAR$rho.pred)
dev.off()

pdf(file="Posteriors.N_hSDM.Nmixture.iCAR.pdf")
plot(mod.hSDM.Nmixture.iCAR$N.pred)
dev.off()

#= Check that Ns are correctly estimated
M <- as.matrix(mod.hSDM.Nmixture.iCAR$N.pred)
N.est <- apply(M,2,mean)
Y.by.cell <- tapply(Data.obs$Y,Data.obs$cells,mean)
pdf(file="Check_N.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(Y.by.cell, N.est) ## More individuals are expected in reality (N > Y) due to detection process
abline(a=0,b=1,col="red")
plot(N[visited.cells], N.est) ## N are well estimated
abline(a=0,b=1,col="red")
cor(N[visited.cells], N.est) ## Close to 1
dev.off()

#= Summary plots

# rho
r.rho <- r.rho.pred <- r.visited <- Landscape
r.rho[] <- rho
rho.pred <- apply(mod.hSDM.Nmixture.iCAR$rho.pred,2,mean)
r.rho.pred[] <- rho.pred
r.visited[] <- 0
r.visited[visited.cells] <- tapply(Data.obs$Y,Data.obs$cells,mean)
# lambda
r.lambda <- Landscape
r.lambda[] <- mod.hSDM.Nmixture.iCAR$lambda.pred

pdf(file="Summary_hSDM.Nmixture.iCAR.pdf")
par(mfrow=c(3,2))
plot(r.rho, main="rho target")
plot(r.visited,main="Visited cells and abundance")
plot(Data.obs[Y>0,],add=TRUE,pch=16,cex=0.5)
plot(r.rho.pred, main="rho estimated")
plot(rho[visited.cells],rho.pred[visited.cells],
     xlab="rho target",
     ylab="rho estimated")
points(rho[notvisited.cells],rho.pred[notvisited.cells],pch=16,col="blue")
legend(x=-4,y=3.5,legend="Unvisited cells",col="blue",pch=16,bty="n")
abline(a=0,b=1,col="red")
plot(r.lambda,main="Predicted abundance")
plot(Data.obs[Y>0,],add=TRUE,pch=16,cex=0.5)
dev.off()

}

}

\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{N-mixture models}
\keyword{spatial correlation}
\keyword{intrinsic CAR model}
\keyword{conditional autoregressive model}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}