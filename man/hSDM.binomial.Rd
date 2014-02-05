\name{hSDM.binomial}
\alias{hSDM.binomial}

\title{Binomial logistic regression model}

\description{The \code{hSDM.binomial} function performs a Binomial
  logistic regression in a Bayesian framework. The function calls
  a Gibbs sampler written in C code which uses an adaptive Metropolis
  algorithm to estimate the conditional posterior distribution of
  model's parameters.}

\usage{hSDM.binomial(presences, trials, suitability, data,
suitability.pred = NULL, burnin = 5000, mcmc = 10000, thin = 10,
beta.start, mubeta = 0, Vbeta = 1e+06, seed = 1234, verbose = 1, save.p
= 0)}

\arguments{
  
  \item{presences}{A vector indicating the number of successes (or
  presences) for each observation.}

  \item{trials}{A vector indicating the number of trials for each
  observation. \eqn{t_n} should be superior or equal to \eqn{y_n}, the
  number of successes for observation \eqn{n}. If \eqn{t_n=0}, then
  \eqn{y_n=0}.}

  \item{suitability}{A one-sided formula of the form '~x1+...+xp' with p
  terms specifying the explicative variables for the suitability process
  of the model.}
  
  \item{data}{A data frame containing the model's explicative
    variables.}

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
    \eqn{D=-2\log(\prod_i P(y_i|\beta,t_i))}{%
      D=-2log(prod_i P(y_i|beta,t_i))}, is also provided.}

  \item{prob.p.pred}{If \code{save.p} is set to 0 (default),
    \code{prob.p.pred} is the predictive posterior mean of the
    probability associated to the suitability process for each
    prediction. If \code{save.p} is set to 1, \code{prob.p.pred} is an
    \code{mcmc} object with sampled values of the probability associated
    to the suitability process for each prediction.}
  
  \item{prob.p.latent}{Predictive posterior mean of the probability
    associated to the suitability process for each observation.}
      
}

\details{

  We model an ecological process where the presence or absence
  of the species is explained by habitat suitability.

  \bold{Ecological process:}
  \deqn{y_i \sim \mathcal{B}inomial(\theta_i,t_i)}{y_i ~ Binomial(theta_i,t_i)}
  \deqn{logit(\theta_i) = X_i \beta}{logit(theta_i) = X_i * beta}
  
}

\references{

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.
  
  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
  Applications}, 16, 33-50.

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
# hSDM.binomial()
# Example with simulated data
#==============================================

#============
#== Preambule
library(hSDM)

#==================
#== Data simulation

# Set seed for repeatability
seed <- 1234

# Constants
nobs <- 1000
set.seed(seed)
trials <- rpois(nobs,5)	# Number of trial associated to each observation

# Covariates for "suitability" process
set.seed(seed)
X1 <- rnorm(n=nobs,0,1)
set.seed(2*seed)
X2 <- rnorm(n=nobs,0,1)
X <- cbind(rep(1,nobs),X1,X2)

# Target parameters
beta.target <- matrix(c(0.2,0.1,0.1),ncol=1) # fixed effects

#== Simulating latent variables

# Suitability
logit.theta <- X \%*\% beta.target
theta <- inv.logit(logit.theta)

#== Simulating response variable
set.seed(seed)
Y <- rbinom(n=nobs,size=trials,prob=theta)

#== Data-set
Data <- data.frame(Y,trials,X1,X2)
str(Data)
par(mfrow=c(2,1))
plot(Data$X1,logit.theta)
plot(Data$X2,logit.theta)

#==================================
#== Statistical modelling

#== glm resolution to compare
mod.glm <- glm(cbind(Y,trials-Y)~X1+X2,family="binomial",data=Data)
summary(mod.glm)
beta.hat.glm <- coef(mod.glm)

#== hSDM
mod.hSDM.binomial <- hSDM.binomial(presences=Data$Y,
                                   trials=Data$trials,
                                   suitability=~X1+X2,
                                   data=Data,
                                   suitability.pred=NULL,
                                   burnin=1000, mcmc=1000, thin=1,
                                   beta.start=beta.hat.glm,
                                   mubeta=0, Vbeta=1.0E6,
                                   seed=1234, verbose=1, save.p=0)

#==========
#== Outputs
summary(mod.hSDM.binomial$mcmc)
pdf(file="Posteriors_hSDM.binomial.pdf")
plot(mod.hSDM.binomial$mcmc)
dev.off()
summary(mod.hSDM.binomial$prob.p.latent)
summary(mod.hSDM.binomial$prob.p.pred)

}

}

\keyword{Binomial logistic regression}
\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
