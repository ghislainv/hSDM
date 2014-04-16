\name{hSDM.siteocc}
\alias{hSDM.siteocc}

\title{Site occupancy model}

\description{The \code{hSDM.siteocc} function can be used to model
  species distribution including different processes in a hierarchical
  Bayesian framework: a \eqn{\mathcal{B}ernoulli}{Bernoulli} suitability
  process (refering to environmental suitability) and a
  \eqn{\mathcal{B}ernoulli}{Bernoulli} observability process (refering
  to various ecological and methodological issues explaining the species
  detection). The \code{hSDM.siteocc} function calls a Gibbs sampler
  written in C code which uses a Metropolis algorithm to estimate the
  conditional posterior distribution of hierarchical model's
  parameters.}

\usage{hSDM.siteocc(# Observations
                     presences, observability, spatial.entity, data.observability,
                     # Habitat
                     suitability, data.suitability,
                     # Predictions
                     suitability.pred = NULL,
                     # Chains
                     burnin = 1000, mcmc = 1000, thin = 1,
                     # Starting values
                     beta.start,
                     gamma.start,
                     # Priors
                     mubeta = 0, Vbeta = 1.0E6,
                     mugamma = 0, Vgamma = 1.0E6,
                     # Various
                     seed = 1234, verbose = 1,
                     save.p = 0)}

\arguments{
  
  \item{presences}{A vector indicating the presence/absence for each
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

  \item{theta.pred}{If \code{save.p} is set to 0 (default),
    \code{theta.pred} is the predictive posterior mean of the
    probability associated to the suitability process for each
    prediction. If \code{save.p} is set to 1, \code{theta.pred} is an
    \code{mcmc} object with sampled values of the probability associated
    to the suitability process for each prediction.}
  
  \item{theta.latent}{Predictive posterior mean of the probability
    associated to the suitability process for each observation.}
  
  \item{delta.latent}{Predictive posterior mean of the probability
    associated to the observability process for each observation.}
  
}

\details{

  The model integrates two processes, an ecological process associated
  to the presence or absence of the species due to habitat suitability
  and an observation process that takes into account the fact that the
  probability of detection of the species is inferior to one.

  \bold{Ecological process:}
  \deqn{z_i \sim \mathcal{B}ernoulli(\theta_i)}{z_i ~ Bernoulli(theta_i)}
  \deqn{logit(\theta_i) = X_i \beta}{logit(theta_i) = X_i * beta}
   
  \bold{Observation process:}
  \deqn{y_i \sim \mathcal{B}ernoulli(z_i * \delta_i)}{y_i ~ Bernoulli(z_i
    * delta_i)}
  \deqn{logit(\delta_i) = W_i \gamma}{logit(delta_i) = W_i * gamma} 

}

\references{

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.
  
  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
    Applications}, 16, 33-50.

  MacKenzie, D. I.; Nichols, J. D.; Lachman, G. B.; Droege, S.; Andrew
  Royle, J. and Langtimm, C. A. (2002) Estimating site occupancy rates when
  detection probabilities are less than one. \emph{Ecology}, 83,
  2248-2255.

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
# hSDM.siteocc()
# Example with simulated data
#==============================================

#============
#== Preambule
library(hSDM)

#==================
#== Data simulation

# Set seed for repeatability
seed <- 1234

# Number of observations
nobs <- 1000

# Target parameters
beta.target <- matrix(c(0.2,0.5,0.5),ncol=1)
gamma.target <- matrix(c(1),ncol=1) #= This means that theta.2~0.75: 75\% of chance to detect the species
## Uncomment if you want covariates on the observability process
## gamma.target <- matrix(c(0.2,0.5,0.5),ncol=1)

# Covariates for "suitability" process
set.seed(seed)
X1 <- rnorm(n=nobs,0,1)
set.seed(2*seed)
X2 <- rnorm(n=nobs,0,1)
X <- cbind(rep(1,nobs),X1,X2)

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
logit.theta.1 <- X \%*\% beta.target
theta.1 <- inv.logit(logit.theta.1)
set.seed(seed)
y.1 <- rbinom(nobs,1,theta.1)

# Observability
set.seed(seed)
trials <- rpois(nobs,5)	# Number of trial associated to each observation
trials[trials==0] <- 1
logit.theta.2 <- W \%*\% gamma.target
theta.2 <- inv.logit(logit.theta.2)
set.seed(seed)
y.2 <- rbinom(nobs,trials,theta.2)

#== Simulating response variable
Y <- y.2*y.1

#== Data-set
Data <- data.frame(Y,trials,X1,X2)
## Uncomment if you want covariates on the observability process 
## Data <- data.frame(Y,trials,X1,X2,W1,W2)

#==================================
#== siteocc model

mod.hSDM.siteocc <- hSDM.siteocc(presences=Data$Y,
                         trials=Data$trials,
                         suitability=~X1+X2,
                         observability=~1, #=~1+W1+W2 if covariates
                         data=Data,
                         suitability.pred=NULL,
                         burnin=1000, mcmc=1000, thin=5,
                         beta.start=0,
                         gamma.start=0,
                         mubeta=0, Vbeta=1.0E6,
                         mugamma=0, Vgamma=1.0E6,
                         seed=1234, verbose=1,
                         save.p=0)

#==========
#== Outputs
pdf(file="Posteriors_hSDM.siteocc.pdf")
plot(mod.hSDM.siteocc$mcmc)
dev.off()
summary(mod.hSDM.siteocc$prob.p.latent)
summary(mod.hSDM.siteocc$prob.q.latent)
summary(mod.hSDM.siteocc$prob.p.pred)

}

}

\keyword{biodiversity}
\keyword{species distribution models}
\keyword{site-occupancy models}
\keyword{hierarchical Bayesian models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
