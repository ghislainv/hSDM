\name{hSDM.poisson}
\alias{hSDM.poisson}

\title{The \code{hSDM.poisson} function performs a Poisson regression in
  a Bayesian framework.}

\description{The \code{hSDM.poisson} function calls a Gibbs sampler
  written in C code which uses a Metropolis algorithm to estimate the
  conditional posterior distribution of model's parameters.}

\usage{hSDM.poisson(counts, suitability, data, burnin = 5000, mcmc =
10000, thin = 10, beta.start, mubeta = 0, Vbeta = 1e+06, seed = 1234,
verbose = 1)}

\arguments{
  
  \item{counts}{A vector indicating the count (or abundance) for each
  observation.}

  \item{suitability}{A one-sided formula of the form '~x1+...+xp' with p
  terms specifying the explicative covariates for the suitability
  process of the model.}
  
  \item{data}{A data frame containing the model's explicative
  variables.}

  \item{burnin}{The number of burnin iterations for the sampler.}
    
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to
    \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and
    superior or equal to 100 so that the progress bar can be displayed.}
    
  \item{thin}{The thinning interval used in the simulation. The number
    of mcmc iterations must be divisible by this value.} 

  \item{beta.start}{Starting values for beta parameters of the
  suitability process. If \code{beta.start} takes a scalar value, then
  that value will serve for all of the betas.}

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
 
  \item{seed}{The seed for the random number generator. Default to 1234.}
  
  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the
    Gibbs sampler.}

}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_i P(y_i,n_i|\beta))}{%
      D=-2log(prod_i P(y_i,n_i|beta))}, is also
    provided.}
  
  \item{prob.pred.p}{Predictive posterior mean of the probability
    associated to the suitability process for each spatial cell.}
    
}

\references{

  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
  Applications}, 16, 33-50.

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.

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
# hSDM.hierarchical.poisson()
# Example with simulated data
#==============================================

#== Preambule
library(hSDM)

#==================
#== Data simulation

# Set seed for repeatability
set.seed(1234)

# Constants
nobs <- 1000

# Covariates for "suitability" process
X1 <- rnorm(n=nobs,0,1)
X2 <- rnorm(n=nobs,0,1)
X <- cbind(rep(1,nobs),X1,X2)

# Target parameters
beta.target <- matrix(c(0.2,0.1,0.1),ncol=1) # fixed effects

#== Simulating latent variables

# Suitability
log.theta <- X \%*\% beta.target
theta <- exp(log.theta)

#== Simulating response variable
Y <- rpois(n=nobs,lambda=theta)

#== Data-set
Data <- data.frame(Y,X1,X2)
str(Data)
par(mfrow=c(2,1))
plot(Data$X1,log.theta)
plot(Data$X2,log.theta)

#========================
#== Statistical modelling

#== glm resolution to compare
mod.glm <- glm(Y~X1+X2,family="poisson",data=Data)
summary(mod.glm)
beta.hat.glm <- coef(mod.glm)

#== hSDM
mod.hSDM.poisson <- hSDM.poisson(counts=Data$Y,
                                 suitability=~X1+X2,
                                 data=Data, burnin=1000,
                                 mcmc=1000, thin=1,
                                 beta.start=beta.hat.glm,
                                 mubeta=0, Vbeta=1.0E6,
                                 seed=1234, verbose=1)

#==========
#== Outputs
summary(mod.hSDM.poisson$mcmc)
pdf(file="Posteriors_hSDM.poisson.pdf")
plot(mod.hSDM.poisson$mcmc)
dev.off()
summary(mod.hSDM.poisson$prob.p.pred)

}

}

\keyword{Poisson regression}
\keyword{biodiversity}
\keyword{species distribution models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
