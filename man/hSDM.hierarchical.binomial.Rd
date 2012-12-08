\name{hSDM.hierarchical.binomial}
\alias{hSDM.hierarchical.binomial}

\title{The \code{hSDM.hierarchical.binomial} function can be used to
  model species distribution including different processes in a
  hierarchical Bayesian framework: (i) a Bernoulli suitability process
  (refering to environmental suitability) which takes into account the
  spatial dependence of the observations, (ii) an alteration process
  (refering to anthropogenic disturbances), and (iii) a Binomial
  observability process (refering to various ecological and
  methodological issues explaining the species presence). The model
  reduces to a ZIB (Zero-Inflated Binomial) model with spatial
  dependence if the alteration is set to zero.}

\description{The \code{hSDM.hierarchical.binomial} function calls a
  Gibbs sampler written in C code which uses a Metropolis algorithm to
  estimate the conditional posterior distribution of hierarchical
  model's parameters.}

\usage{hSDM.hierarchical.binomial(presences, trials, suitability, cells,
n.neighbors, neighbors, alteration, observability, data, burnin = 5000,
mcmc = 10000, thin = 10, beta.start, gamma.start, Vrho.start, mubeta =
0, Vbeta = 1e+06, mugamma = 0, Vgamma = 1e+06, priorVrho = "1/Gamma",
shape = 0.5, rate = 0.0005, Vrho.max=1000, seed = 1234, verbose = 1)}

\arguments{
  
  \item{presences}{A vector indicating the number of successes (or
  presences) for each observation.}

  \item{trials}{A vector indicating the number of trials for each
  observation. \eqn{t_n} should be superior or equal to \eqn{y_n}, the
  number of successes for observation \eqn{n}. If \eqn{t_n=0}, then
  \eqn{y_n=0}.}

  \item{suitability}{A one-sided formula of the form '~x1+...+xp' with p
  terms specifying the explicative variables for the suitability
  process of the model.}

  \item{cells}{A vector indicating the spatial cell identifier (from 1
  to total number of cell) for each observation. Several observations
  can occur in one spatial cell.}

  \item{n.neighbors}{A vector of integers that indicates the number of
    neighbors (adjacent cells) of each spatial cell.}
  
  \item{neighbors}{A vector of integers indicating the neighbors
  (adjacent cells) of each spatial cell. Must be of the form c(neighbors
  of cell 1, neighbors of cell 2, ... , neighbors of the last
  cell). Length of the \code{neighbors} vector should be equal to
  sum(data$num).}
  
  \item{alteration}{A vector indicating the proportion of area in the
    spatial cell which is transformed (by anthropogenic activities for
    example) for each observation. Must be between 0 and 1.}
  
  \item{observability}{A one-sided formula of the form '~x1+...+xq' with
    q terms specifying the explicative variables for the observability
    process of the model.}

  \item{data}{A data frame containing the model's variables.}

  \item{burnin}{The number of burnin iterations for the sampler.}
    
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to
    \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and
    superior or equal to 100 so that the progress bar can be displayed.}
    
  \item{thin}{The thinning interval used in the simulation. The number
    of mcmc iterations must be divisible by this value.}

  \item{beta.start}{Starting values for beta parameters.}

  \item{gamma.start}{Starting values for gamma parameters.}
  
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

}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_i P(y_i,n_i|u_i,\beta,\gamma,\rho_i))}{%
      D=-2log(prod_i P(y_i,n_i|u_i,beta,gamma,rho_i))}, is also
    provided.}

  \item{rho.pred}{Predictive posterior mean of the spatial random effect
    associated to each spatial cell.}
  
  \item{prob.pred.p}{Predictive posterior mean of the probability
    associated to the suitability process for each spatial cell.}
  
  \item{prob.pred.q}{Predictive posterior mean of the probability
  associated to the observability process for each spatial cell.}
  
}

\references{

  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
  Applications}, 16, 33-50.

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.

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
# hSDM.hierarchical.binomial()
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

model <- hSDM.hierarchical.binomial(presences=Data$Y,
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
                                    seed=1234, verbose=1)

#==========
#== Outputs
summary(model$mcmc)
pdf(file="Posteriors_hSDM.hierarchical.binomial.pdf")
plot(model$mcmc)
dev.off()
summary(model$prob.p.pred)
summary(model$rho)

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