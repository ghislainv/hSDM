\name{hSDM.site.occupancy}
\alias{hSDM.site.occupancy}

\title{The \code{hSDM.site.occupancy} function can be used to model
  species distribution including different processes in a hierarchical
  Bayesian framework: a Bernoulli suitability process (refering to
  environmental suitability) which takes into account the spatial
  dependence of the observations, and a Binomial observability process
  (refering to various ecological and methodological issues explaining
  the species presence).}

\description{The \code{hSDM.site.occupancy} function calls a
  Gibbs sampler written in C code which uses a Metropolis algorithm to
  estimate the conditional posterior distribution of hierarchical
  model's parameters.}

\usage{hSDM.site.occupancy(presences, trials, suitability,
observability, spatial.entity, data, n.neighbors, neighbors,
suitability.pred=NULL, spatial.entity.pred=NULL, burnin = 5000, mcmc =
10000, thin = 10, beta.start, gamma.start, Vrho.start, mubeta = 0, Vbeta
= 1e+06, mugamma = 0, Vgamma = 1e+06, priorVrho = "1/Gamma", shape =
0.5, rate = 0.0005, Vrho.max=1000, seed = 1234, verbose = 1, save.rho =
0, save.p = 0)}

\arguments{
  
  \item{presences}{A vector indicating the number of successes (or
    presences) for each observation.}

  \item{trials}{A vector indicating the number of trials for each
    observation. \eqn{t_n} should be superior to zero and superior or
    equal to \eqn{y_n}, the number of successes for observation \eqn{n}.}

  \item{suitability}{A one-sided formula of the form '~x1+...+xp' with p
    terms specifying the explicative variables for the suitability
    process of the model.}

  \item{observability}{A one-sided formula of the form '~w1+...+wq' with
    q terms specifying the explicative variables for the observability
    process of the model.}
  
  \item{spatial.entity}{A vector indicating the spatial entity
    identifier (from one to the total number of entities) for each
    observation. Several observations can occur in one spatial entity. A
    spatial entity can be a raster cell for example.}

  \item{data}{A data frame containing the model's variables.}

  \item{n.neighbors}{A vector of integers that indicates the number of
    neighbors (adjacent entities) of each spatial
    entity. \code{length(n.neighbors)} indicates the total number of
    spatial entities.}
  
  \item{neighbors}{A vector of integers indicating the neighbors
    (adjacent entities) of each spatial entity. Must be of the form c(neighbors
    of entity 1, neighbors of entity 2, ... , neighbors of the last
    entity). Length of the \code{neighbors} vector should be equal to
    sum(data$num).}

  \item{suitability.pred}{An optional data frame in which to look for
    variables with which to predict. If NULL, the observations are
    used.}

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

  \item{beta.start}{Starting values for \eqn{\beta}{beta}
    parameters. This can either be a scalar or a p-length vector.}

  \item{gamma.start}{Starting values for \eqn{\beta}{gamma}
    parameters. This can either be a scalar or a q-length vector.}
  
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
  [0,\code{Vrho.max}]. Default to "1/Gamma".}

  \item{shape}{The shape parameter for the Gamma prior on the precision
  of the spatial random effects. Default value is \code{shape=0.05} for
  uninformative prior.}

  \item{rate}{The rate (1/scale) parameter for the Gamma prior on the
    precision of the spatial random effects. Default value is
    \code{rate=0.0005} for uninformative prior.}

  \item{Vrho.max}{Upper bound for the uniform prior of the spatial
    random effect variance. Default to 1000.}
 
  \item{seed}{The seed for the random number generator. Default to
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
    posterior mean is computed and returned in the \code{prob.p.pred}
    vector. Be careful, setting \code{save.p} to 1 might require a large
    amount of memory.}

}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_i P(y_i,n_i|u_i,\beta,\gamma,\rho_i))}{%
      D=-2log(prod_i P(y_i,n_i|u_i,beta,gamma,rho_i))}, is also
    provided.}

  \item{rho.pred}{If \code{save.rho} is set to 0 (default),
    \code{rho.pred} is the predictive posterior mean of the spatial
    random effect associated to each spatial entity. If \code{save.rho} is
    set to 1, \code{rho.pred} is an \code{mcmc} object with sampled
    values for each spatial random effect associated to each spatial
    entity.}

  \item{prob.p.pred}{If \code{save.p} is set to 0 (default),
    \code{prob.p.pred} is the predictive posterior mean of the
    probability associated to the suitability process for each
    prediction. If \code{save.p} is set to 1, \code{prob.p.pred} is an
    \code{mcmc} object with sampled values of the probability associated
    to the suitability process for each prediction.}
  
  \item{prob.p.latent}{Predictive posterior mean of the probability
    associated to the suitability process for each observation.}
  
  \item{prob.q.latent}{Predictive posterior mean of the probability
  associated to the observability process for each observation.}
  
}

\references{

  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
  Applications}, 16, 33-50.

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.

  MacKenzie, D. I.; Nichols, J. D.; Lachman, G. B.; Droege, S.; Andrew
  Royle, J. and Langtimm, C. A. (2002) Estimating site occupancy rates when
  detection probabilities are less than one. \emph{Ecology}, 83, 2248-2255.

  Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive rejection
  Metropolis sampling. \emph{Applied Statistics}, 44, 455-472.

}

\author{
  Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
}

\seealso{
  \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
}

\examples{

\dontrun{

#==============================================
# hSDM.site.occupancy()
# Example with simulated data
#==============================================

#============
#== Preambule
library(mvtnorm)
library(hSDM)

#==================
#== Data simulation

# Set seed for repeatability
set.seed(1234)

# Constants
ncell <- 150 # Number of cells
nobs <- 100*ncell # Number of observation for the *binomial* random variable
trials <- rpois(nobs,5)	# Number of trial associated to each observation
cell <- rep(c(1:ncell),each=nobs/ncell)

# Covariates for "suitability" process
X1 <- rnorm(n=nobs,0,1)
X2 <- rnorm(n=nobs,0,1)
X <- cbind(rep(1,nobs),X1,X2)

# Alteration
U <- runif(n=nobs,min=0,max=1)

# Covariates for "observability" process
W1 <- rnorm(n=nobs,0,1)
W2 <- rnorm(n=nobs,0,1)
W <- cbind(rep(1,nobs),W1,W2)

# Target parameters
beta.target <- matrix(c(0.2,0.1,0.1),ncol=1) # fixed effects
gamma.target <- matrix(c(0.3,0.1,0.1),ncol=1) # fixed effects
Vrho.target <- 10 # Spatial Variance

# Generate symmetric adjacency matrix, A  
A <- matrix(0,ncell,ncell)
A[upper.tri(A,diag=F)] <- rbinom(ncell*(ncell-1)/2,1,.05)
A <- A+t(A) 
n.neighbors <- apply(A,1,sum)
f.adjacent <- function (x) {
  which(x==1)
}
adj <- unlist(apply(A,1,f.adjacent))
  
# Spatial effects
d <- 1	# Spatial dependence parameter = 1 for intrinsic CAR
Q <- diag(n.neighbors)-d*A + diag(.0001,ncell) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos
rho <- c(rmvnorm(1,sigma=covrho)) # Spatial Random Effects
rho <- rho-mean(rho) # Centering rhos on zero

#== Simulating latent variables

# Suitability
logit.theta.1 <- vector()
for (n in 1:nobs) {
  logit.theta.1[n] <- X[n,]\%*\%beta.target+rho[cell[n]]
}
theta.1 <- inv.logit(logit.theta.1)
y.1 <- rbinom(nobs,1,theta.1)

# Alteration
u <- rbinom(nobs,1,U)

# Observability
logit.theta.2 <- W\%*\%gamma.target
theta.2 <- inv.logit(logit.theta.2)
y.2 <- rbinom(nobs,trials,theta.2)

#== Simulating response variable
Y <- y.2*(1-u)*y.1

#== Data-set
Data <- data.frame(Y,trials,cell,X1,X2,W1,W2,U)
str(Data)

#==================================
#== Statistical modelling with hSDM

model <- hSDM.site.occupancy(presences=Data$Y,
                                    trials=Data$trials,
                                    suitability=~X1+X2,
                                    cells=Data$cell,
                                    n.neighbors=n.neighbors,
                                    neighbors=adj,
                                    alteration=Data$U,
                                    observability=~W1+W2,
                                    data=Data, burnin=500,
                                    mcmc=500, thin=1,
                                    beta.start=0,
                                    gamma.start=0,
                                    Vrho.start=1,
                                    priorVrho="1/Gamma",
                                    #priorVrho="Uniform",
                                    #priorVrho=10,
                                    mubeta=0, Vbeta=1.0E6,
                                    mugamma=0, Vgamma=1.0E6,
                                    shape=0.5, rate=0.0005,
                                    Vrho.max=1000,
                                    seed=1234, verbose=1, save.rho=0)

#==========
#== Outputs
summary(model$mcmc)
pdf(file="Posteriors_hSDM.site.occupancy.pdf")
plot(model$mcmc)
dev.off()
summary(model$prob.p.pred)
summary(model$rho.pred)

}

}

\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{spatial correlation}
\keyword{intrinsic CAR model}
\keyword{conditional autoregressive model}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
\keyword{ZIB models}