\name{hSDM.Nmixture}
\alias{hSDM.Nmixture}

\title{N-mixture model}

\description{The \code{hSDM.Nmixture} function can be used to model
  species distribution including different processes in a hierarchical
  Bayesian framework: a \eqn{\mathcal{P}oisson}{Poisson} suitability
  process (refering to environmental suitability) and a
  \eqn{\mathcal{B}inomial}{Binomial} observability process (refering to
  various ecological and methodological issues explaining the species
  detection). The \code{hSDM.Nmixture} function calls a Gibbs
  sampler written in C code which uses an adaptive Metropolis algorithm
  to estimate the conditional posterior distribution of hierarchical
  model's parameters.}

\usage{hSDM.Nmixture(# Observations
                     counts, observability, spatial.entity, data.observability,
                     # Habitat
                     suitability, data.suitability,
                     # Predictions
                     suitability.pred = NULL,
                     # Chains
                     burnin = 5000, mcmc = 10000, thin = 10,
                     # Starting values
                     beta.start,
                     gamma.start,
                     # Priors
                     mubeta = 0, Vbeta = 1.0E6,
                     mugamma = 0, Vgamma = 1.0E6,
                     # Various
                     seed = 1234, verbose = 1,
                     save.p = 0, save.N = 0)}

\arguments{
  
  \item{counts}{A vector indicating the count (or abundance) for each
    observation.}

  \item{observability}{A one-sided formula of the form
    \eqn{\sim w_1+...+w_q}{'~w_1+...+w_q'} with \eqn{q} terms specifying the
    explicative variables for the observability process.}
  
  \item{spatial.entity}{A vector indicating the spatial entity
    identifier (from one to the total number of entities) for each
    observation. Several observations can occur in one spatial entity. A
    spatial entity can be a raster cell for example.}

  \item{data.observability}{A data frame containing the model's
    variables for the observability process.}

  \item{suitability}{A one-sided formula of the form
    \eqn{\sim x_1+...+x_p}{'~x_1+...+x_p'} with \eqn{p} terms specifying the
    explicative variables for the suitability process.}

  \item{data.suitability}{A data frame containing the model's
    variables for the suitability process.}

  \item{suitability.pred}{An optional data frame in which to look for
    variables with which to predict. If NULL, the observations are
    used.}

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
 
  \item{seed}{The seed for the random number generator. Default set to
  1234.}
  
  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the
    Gibbs sampler.}

  \item{save.p}{A switch (0,1) which determines whether or not the
    sampled values for predictions are saved. Default is 0: the
    posterior mean is computed and returned in the \code{prob.p.pred}
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

  \item{prob.p.pred}{If \code{save.p} is set to 0 (default),
    \code{prob.p.pred} is the predictive posterior mean of the
    probability associated to the suitability process for each
    prediction. If \code{save.p} is set to 1, \code{prob.p.pred} is an
    \code{mcmc} object with sampled values of the probability associated
    to the suitability process for each prediction.}

  \item{N.pred}{If \code{save.N} is set to 0 (default), \code{N.pred} is
    the posterior mean (rounded to the closest integer) of the latent
    count variable N for each observed cell. If \code{save.N} is set to 1,
    \code{N.pred} is an \code{mcmc} object with sampled values of the
    latent count variable N for each observed cell.}
  
  \item{prob.p.latent}{Predictive posterior mean of the probability
    associated to the suitability process for each observation.}
  
  \item{prob.q.latent}{Predictive posterior mean of the probability
  associated to the observability process for each observation.}
  
}

\details{

  The model integrates two processes, an ecological process associated
  to the presence or absence of the species due to habitat suitability
  and an observation process that takes into account the fact that the
  probability of detection of the species is inferior to one.

  \bold{Ecological process:}
  \deqn{N_i \sim \mathcal{P}oisson(\lambda_i)}{N_i ~ Poisson(lambda_i)}
  \deqn{log(\lambda_i) = X_i \beta}{log(lambda_i) = X_i * beta}
  
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
# hSDM.Nmixture()
# Example with simulated data
#==============================================

#============
#== Preambule
library(hSDM)

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

# Number of observations
nobs <- 300

# Cell vector
set.seed(seed)
cells <- sample(c(1:100),nobs,replace=TRUE)
Ident.cells <- as.numeric(as.factor(cells))
Levels.cells <- sort(unique(cells))
ncell <- length(Levels.cells)

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
log.lambda <- X \%*\% beta.target
lambda <- exp(log.lambda)
set.seed(seed)
N <- rpois(ncell,lambda)

# Observability
logit.theta <- W \%*\% gamma.target
theta <- inv.logit(logit.theta)

#== Simulating response variable
set.seed(seed)
Y <- rbinom(nobs,N[Ident.cells],theta)

#== Data-sets
Data.obs <- data.frame(Y,cells)
## Uncomment if you want covariates on the observability process 
## Data.obs <- data.frame(Y,cells,W1,W2)
Data.suit <- data.frame(Levels.cells,X1,X2)
Data.pred <- Data.suit

#==================================
#== N-mixture model

mod.hSDM.Nmixture <- hSDM.Nmixture(# Observations
                                   counts=Data.obs$Y, observability=~1,
                                   spatial.entity=Data.obs$cells, data.observability=Data.obs,
                                   # Habitat
                                   suitability=~X1+X2, data.suitability=Data.suit,
                                   # Predictions
                                   suitability.pred=Data.pred,
                                   # Chains
                                   burnin=10000, mcmc=10000, thin=10,
                                   # Starting values
                                   beta.start=0,
                                   gamma.start=0,
                                   # Priors
                                   mubeta=0, Vbeta=1.0E6,
                                   mugamma=0, Vgamma=1.0E6,
                                   # Various
                                   seed=1234, verbose=1,
                                   save.p=0, save.N=1)

#==========
#== Outputs

#= Parameter estimates
summary(mod.hSDM.Nmixture$mcmc)

#= MCMC and posteriors
pdf(file="Posteriors_hSDM.Nmixture.pdf")
plot(mod.hSDM.Nmixture$mcmc)
dev.off()

pdf(file="Posteriors.N_hSDM.Nmixture.pdf")
plot(mod.hSDM.Nmixture$N.pred)
dev.off()

#= Check that Ns are correctly estimated
M <- as.matrix(mod.hSDM.Nmixture$N.pred)
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

}

}

\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{N-mixture models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
