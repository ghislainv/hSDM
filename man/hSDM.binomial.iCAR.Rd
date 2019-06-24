\name{hSDM.binomial.iCAR}
\alias{hSDM.binomial.iCAR}

\title{Binomial logistic regression model with CAR process}

\description{The \code{hSDM.binomial.iCAR} function performs a Binomial
  logistic regression model in a hierarchical Bayesian framework. The
  suitability process includes a spatial correlation process. The
  spatial correlation is modelled using an intrinsic CAR model. The
  \code{hSDM.binomial.iCAR} function calls a Gibbs sampler written in C
  code which uses an adaptive Metropolis algorithm to estimate the
  conditional posterior distribution of hierarchical model's
  parameters.}

\usage{hSDM.binomial.iCAR(presences, trials, suitability,
spatial.entity, data, n.neighbors, neighbors, suitability.pred=NULL,
spatial.entity.pred=NULL, burnin = 5000, mcmc = 10000, thin = 10,
beta.start, Vrho.start, mubeta = 0, Vbeta = 1e+06, priorVrho =
"1/Gamma", shape = 0.5, rate = 0.0005, Vrho.max=1000, seed = 1234,
verbose = 1, save.rho = 0, save.p = 0)}

\arguments{
  
  \item{presences}{A vector indicating the number of successes (or
    presences) for each observation.}

  \item{trials}{A vector indicating the number of trials for each
    observation. \eqn{t_i} should be superior to zero and superior or
    equal to \eqn{y_i}, the number of successes for observation \eqn{i}.}

  \item{suitability}{A one-sided formula of the form
    \eqn{\sim x_1+...+x_p}{'~x_1+...+x_p'} with \eqn{p} terms specifying the
    explicative variables for the suitability process.}

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
    (adjacent entities) of each spatial entity. Must be of the form
    c(neighbors of entity 1, neighbors of entity 2, ... , neighbors of
    the last entity). Length of the \code{neighbors} vector should be
    equal to \code{sum(n.neighbors)}.}

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

  \item{beta.start}{Starting values for \eqn{\beta}{beta} parameters of
    the suitability process. This can either be a scalar or a
    \eqn{p}-length vector.}
  
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

  \item{priorVrho}{Type of prior for the variance of the spatial random
  effects. Can be set to a fixed positive scalar, or to an inverse-gamma
  distribution ("1/Gamma") with parameters \code{shape} and \code{rate},
  or to a uniform distribution ("Uniform") on the interval
  [0,\code{Vrho.max}]. Default set to "1/Gamma".}

  \item{shape}{The shape parameter for the Gamma prior on the precision
  of the spatial random effects. Default value is \code{shape=0.5} for
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
    posterior mean is computed and returned in the \code{theta.pred}
    vector. Be careful, setting \code{save.p} to 1 might require a large
    amount of memory.}

}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_i P(y_i|...))}{D=-2log(prod_i P(y_i|...))}, is
    also provided.}

  \item{rho.pred}{If \code{save.rho} is set to 0 (default),
    \code{rho.pred} is the predictive posterior mean of the spatial
    random effect associated to each spatial entity. If \code{save.rho} is
    set to 1, \code{rho.pred} is an \code{mcmc} object with sampled
    values for each spatial random effect associated to each spatial
    entity.}

  \item{theta.pred}{If \code{save.p} is set to 0 (default),
    \code{theta.pred} is the predictive posterior mean of the
    probability associated to the suitability process for each
    prediction. If \code{save.p} is set to 1, \code{theta.pred} is an
    \code{mcmc} object with sampled values of the probability associated
    to the suitability process for each prediction.}
  
  \item{theta.latent}{Predictive posterior mean of the probability
    associated to the suitability process for each observation.}
  
}

\details{

  We model an ecological process where the presence or absence of the
  species is explained by habitat suitability. The ecological process
  includes an intrinsic conditional autoregressive (iCAR) model for
  spatial autocorrelation between observations, assuming that the
  probability of presence of the species at one site depends on the
  probability of presence of the species on neighboring sites.

  \bold{Ecological process:}
  \deqn{y_i \sim \mathcal{B}inomial(\theta_i,t_i)}{y_i ~ Binomial(theta_i,t_i)}
  \deqn{logit(\theta_i) = X_i \beta + \rho_{j(i)}}{logit(theta_i) = X_i *
    beta + rho_j(i)}
  \eqn{\rho_j}{rho_j}: spatial random effect
  
  \eqn{j(i)}: index of the spatial entity for observation \eqn{i}.
  
  \bold{Spatial autocorrelation:}
  
  An intrinsic conditional autoregressive model (iCAR) is assumed:
  \deqn{\rho_j \sim \mathcal{N}ormal(\mu_j,V_{\rho} / n_j)}{rho_j ~
    Normal(mu_j,V_rho / n_j)}
  \eqn{\mu_j}{mu_j}: mean of \eqn{\rho_{j'}}{rho_j'} in the
  neighborhood of \eqn{j}.
  
  \eqn{V_{\rho}}{V_rho}: variance of the spatial random effects.
  
  \eqn{n_j}{n_j}: number of neighbors for spatial entity \eqn{j}.
  
}

\references{

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.
   
  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
    Applications}, 16, 33-50.

  Lichstein, J. W.; Simons, T. R.; Shriner, S. A. & Franzreb,
  K. E. (2002) Spatial autocorrelation and autoregressive models in ecology
  \emph{Ecological Monographs}, 72, 445-463.

  Diez, J. M. & Pulliam, H. R. (2007) Hierarchical analysis of species
  distributions and abundance across environmental gradients
  \emph{Ecology}, 88, 3144-3152.

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
# hSDM.binomial.iCAR()
# Example with simulated data
#==============================================

#=================
#== Load libraries                          
library(hSDM)
library(raster)
library(sp)

#===================================
#== Multivariate normal distribution
rmvn <- function(n, mu = 0, V = matrix(1), seed=1234) {
    p <- length(mu)
    if (any(is.na(match(dim(V), p)))) {
        stop("Dimension problem!")
    }
    D <- chol(V)
    set.seed(seed)
    t(matrix(rnorm(n*p),ncol=p)\%*\%D+rep(mu,rep(n,p)))
}

#==================
#== Data simulation

#= Set seed for repeatability
seed <- 1234

#= Landscape
xLand <- 30
yLand <- 30 
Landscape <- raster(ncol=xLand,nrow=yLand,crs='+proj=utm +zone=1')
Landscape[] <- 0
extent(Landscape) <- c(0,xLand,0,yLand)
coords <- coordinates(Landscape)
ncells <- ncell(Landscape)

#= Neighbors
neighbors.mat <- adjacent(Landscape, cells=c(1:ncells), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]

#= Generate symmetric adjacency matrix, A
A <- matrix(0,ncells,ncells)
index.start <- 1
for (i in 1:ncells) {
    index.end <- index.start+n.neighbors[i]-1
    A[i,adj[c(index.start:index.end)]] <- 1
    index.start <- index.end+1
}

#= Spatial effects
Vrho.target <- 5
d <- 1	# Spatial dependence parameter = 1 for intrinsic CAR
Q <- diag(n.neighbors)-d*A + diag(.0001,ncells) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos
set.seed(seed)
rho <- c(rmvn(1,mu=rep(0,ncells),V=covrho,seed=seed)) # Spatial Random Effects
rho <- rho-mean(rho) # Centering rhos on zero

#= Raster and plot spatial effects
r.rho <- rasterFromXYZ(cbind(coords,rho))
plot(r.rho)

#= Sample the observation sites in the landscape
nsite <- 250
set.seed(seed)
x.coord <- runif(nsite,0,xLand)
set.seed(2*seed)
y.coord <- runif(nsite,0,yLand)
sites.sp <- SpatialPoints(coords=cbind(x.coord,y.coord))
cells <- extract(Landscape,sites.sp,cell=TRUE)[,1]

#= Number of visits associated to each observation point
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

#= Ecological process (suitability)
set.seed(seed)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
beta.target <- c(-1,1,-1)
logit.theta <- X \%*\% beta.target + rho[cells]
theta <- inv.logit(logit.theta)
set.seed(seed)
Y <- rbinom(nsite,visits,theta)

#= Relative importance of spatial random effects
RImp <- mean(abs(rho[cells])/abs(X \%*\% beta.target))
RImp

#= Data-sets
data.obs <- data.frame(Y,visits,x1,x2,cell=cells)

#==================================
#== Site-occupancy model

Start <- Sys.time() # Start the clock
mod.hSDM.binomial.iCAR <- hSDM.binomial.iCAR(presences=data.obs$Y,
                                               trials=data.obs$visits,
                                               suitability=~x1+x2,
                                               spatial.entity=data.obs$cell,
                                               data=data.obs,
                                               n.neighbors=n.neighbors,
                                               neighbors=adj,
                                               suitability.pred=NULL,
                                               spatial.entity.pred=NULL,
                                               burnin=5000, mcmc=5000, thin=5,
                                               beta.start=0,
                                               Vrho.start=1,
                                               mubeta=0, Vbeta=1.0E6,
                                               priorVrho="1/Gamma",
                                               shape=0.5, rate=0.0005,
                                               seed=1234, verbose=1,
                                               save.rho=1, save.p=0)
Time.hSDM <- difftime(Sys.time(),Start,units="sec") # Time difference

#= Computation time
Time.hSDM

#==========
#== Outputs

#= Parameter estimates 
summary(mod.hSDM.binomial.iCAR$mcmc)
pdf("Posteriors_hSDM.binomial.iCAR.pdf")
plot(mod.hSDM.binomial.iCAR$mcmc)
dev.off()

#= Predictions
summary(mod.hSDM.binomial.iCAR$theta.latent)
summary(mod.hSDM.binomial.iCAR$theta.pred)
pdf(file="Pred-Init.pdf")
plot(theta,mod.hSDM.binomial.iCAR$theta.pred)
abline(a=0,b=1,col="red")
dev.off()

#= Summary plots for spatial random effects

# rho.pred
rho.pred <- apply(mod.hSDM.binomial.iCAR$rho.pred,2,mean)
r.rho.pred <- rasterFromXYZ(cbind(coords,rho.pred))

# plot
pdf(file="Summary_hSDM.binomial.iCAR.pdf")
par(mfrow=c(2,2))
# rho target
plot(r.rho, main="rho target")
plot(sites.sp,add=TRUE)
# rho estimated
plot(r.rho.pred, main="rho estimated")
# correlation and "shrinkage"
Levels.cells <- sort(unique(cells))
plot(rho[-Levels.cells],rho.pred[-Levels.cells],
     xlim=range(rho),
     ylim=range(rho),
     xlab="rho target",
     ylab="rho estimated")
points(rho[Levels.cells],rho.pred[Levels.cells],pch=16,col="blue")
legend(x=-3,y=4,legend="Visited cells",col="blue",pch=16,bty="n")
abline(a=0,b=1,col="red")
dev.off()

}

}

\keyword{Binomial logistic regression}
\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{spatial correlation}
\keyword{intrinsic CAR model}
\keyword{conditional autoregressive model}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}