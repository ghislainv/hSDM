\name{hSDM.Nmixture.K}
\alias{hSDM.Nmixture.K}

\title{N-mixture model with K, the maximal theoretical abundance}

\description{The \code{hSDM.Nmixture.K} function can be used to model
  species distribution including different processes in a hierarchical
  Bayesian framework: a \eqn{\mathcal{P}oisson}{Poisson} suitability
  process (refering to environmental suitability explaining abundance)
  and a \eqn{\mathcal{B}inomial}{Binomial} observability process
  (refering to various ecological and methodological issues explaining
  species detection). The \code{hSDM.Nmixture.K} function calls a Gibbs
  sampler written in C code which uses an adaptive Metropolis algorithm
  to estimate the conditional posterior distribution of hierarchical
  model's parameters. K is the maximal theoretical abundance sensus
  Royle 2004.}

\usage{hSDM.Nmixture.K(# Observations
                       counts, observability, site, data.observability,
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
                       K,
                       seed = 1234, verbose = 1,
                       save.p = 0)}

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

  \item{K}{Maximal theoretical abundance sensus Royle 2004. It
    corresponds to the integer upper index of integration for
    N-mixture. This should be set high enough so that it does not affect
    the parameter estimates. Note that computation time will increase
    with K.}
 
  \item{seed}{The seed for the random number generator. Default set to
  1234.}
  
  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the
    Gibbs sampler.}

  \item{save.p}{A switch (0,1) which determines whether or not the
    sampled values for predictions are saved. Default is 0: the
    posterior mean is computed and returned in the \code{lambda.pred}
    vector. Be careful, setting \code{save.p} to 1 might require a large
    amount of memory.}

}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_{it} P(y_{it},N_i|...))}{%
      D=-2log(prod_it P(y_it,N_i|...))}, is also
    provided.}

  \item{lambda.pred}{If \code{save.p} is set to 0 (default),
    \code{lambda.pred} is the predictive posterior mean of the abundance
    associated to the suitability process for each prediction. If
    \code{save.p} is set to 1, \code{lambda.pred} is an \code{mcmc}
    object with sampled values of the abundance associated to the
    suitability process for each prediction.}
  
  \item{lambda.latent}{Predictive posterior mean of the abundance
    associated to the suitability process for each observation.}
  
  \item{delta.latent}{Predictive posterior mean of the probability
    associated to the observability process for each observation.}
  
}

\details{

  The model integrates two processes, an ecological process associated
  to the abundance of the species due to habitat suitability and an
  observation process that takes into account the fact that the
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
# hSDM.Nmixture.K()
# Example with simulated data
#==============================================

#=================
#== Load libraries
library(hSDM)

#==================
#== Data simulation

# Number of observation sites
nsite <- 200

#= Set seed for repeatability
seed <- 4321

#= Ecological process (suitability)
set.seed(seed)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
beta.target <- c(-1,1,-1) # Target parameters
log.lambda <- X \%*\% beta.target
lambda <- exp(log.lambda)
set.seed(seed)
N <- rpois(nsite,lambda)

#= Number of visits associated to each observation point
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1
# Vector of observation points
sites <- vector()
for (i in 1:nsite) {
    sites <- c(sites,rep(i,visits[i]))
}

#= Observation process (detectability)
nobs <- sum(visits)
set.seed(seed)
w1 <- rnorm(nobs,0,1)
set.seed(2*seed)
w2 <- rnorm(nobs,0,1)
W <- cbind(rep(1,nobs),w1,w2)
gamma.target <- c(-1,1,-1) # Target parameters
logit.delta <- W \%*\% gamma.target
delta <- inv.logit(logit.delta)
set.seed(seed)
Y <- rbinom(nobs,N[sites],delta)

#= Data-sets
data.obs <- data.frame(Y,w1,w2,site=sites)
data.suit <- data.frame(x1,x2)

#================================
#== Parameter inference with hSDM

Start <- Sys.time() # Start the clock
mod.hSDM.Nmixture.K <- hSDM.Nmixture.K(# Observations
                           counts=data.obs$Y,
                           observability=~w1+w2,
                           site=data.obs$site,
                           data.observability=data.obs,
                           # Habitat
                           suitability=~x1+x2,
                           data.suitability=data.suit,
                           # Predictions
                           suitability.pred=NULL,
                           # Chains
                           burnin=5000, mcmc=5000, thin=5,
                           # Starting values
                           beta.start=0,
                           gamma.start=0,
                           # Priors
                           mubeta=0, Vbeta=1.0E6,
                           mugamma=0, Vgamma=1.0E6,
                           # Various
                           K=max(data.obs$Y)*2,
                           seed=1234, verbose=1,
                           save.p=0)
Time.hSDM <- difftime(Sys.time(),Start,units="sec") # Time difference

#= Computation time
Time.hSDM

#==========
#== Outputs

#= Parameter estimates
summary(mod.hSDM.Nmixture.K$mcmc)
pdf(file="Posteriors_hSDM.Nmixture.K.pdf")
plot(mod.hSDM.Nmixture.K$mcmc)
dev.off()

#= Predictions
summary(mod.hSDM.Nmixture.K$lambda.latent)
summary(mod.hSDM.Nmixture.K$delta.latent)
summary(mod.hSDM.Nmixture.K$lambda.pred)
pdf(file="Pred-Init.K.pdf")
plot(lambda,mod.hSDM.Nmixture.K$lambda.pred)
abline(a=0,b=1,col="red")
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
