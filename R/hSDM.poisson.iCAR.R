####################################################################
##
## hSDM.poisson.iCAR.R
##
####################################################################
##
## Original code by Ghislain Vieilledent, October 2011
## CIRAD UR B&SEF
## ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
##
####################################################################
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2011 Ghislain Vieilledent
## 
####################################################################
##
## Revisions: 
## - G. Vieilledent, on November 19th 2012
##
####################################################################


hSDM.poisson.iCAR <- function (counts,
                               suitability,
                               cells,
                               n.neighbors,
                               neighbors,  
                               data, burnin=5000,
                               mcmc=10000, thin=10, 
                               beta.start,
                               Vrho.start,
                               mubeta=0, Vbeta=1.0E6,
                               priorVrho="1/Gamma",
                               shape=0.5, rate=0.0005,
                               Vrho.max=1000,
                               seed=1234, verbose=1)

{   
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
   
  #======== 
  # Form response, covariate matrices and model parameters
  #========

  #= Response
  Y <- counts
  nobs <- length(Y)
  #= Suitability
  mf.suit <- model.frame(formula=suitability,data=data)
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  #= Spatial correlation
  cells <- cells
  ncell <- length(unique(cells))
  n.neighbors <- n.neighbors
  neighbors <- neighbors
  #= Model parameters
  np <- ncol(X)
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin

  #========== 
  # Check data
  #==========
  check.Y.poisson(Y)
  check.X(X,nobs)
  check.cells(cells,nobs)
  check.neighbors(n.neighbors,ncell,neighbors)
  
  #========
  # Initial starting values for M-H
  #========
  beta.start <- form.beta.start(beta.start,np)
  rho.start <- rep(0,ncell) # Starting values for spatial random effects set to zero.
  Vrho.start <- check.Vrho.start(Vrho.start)
  
  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mubeta,np)
  Vbeta <- check.Vbeta(Vbeta,np)
  check.ig.prior(shape,rate)
  Vrho.max <- check.Vrho.max(Vrho.max)
  priorVrho <- form.priorVrho(priorVrho)

  #========
  # Parameters to save
  #========
  beta <- rep(beta.start,nsamp)
  rho_pred <- rho.start
  Vrho <- rep(Vrho.start,nsamp)
  prob_p_pred <- rep(0,nobs)
  Deviance <- rep(0,nsamp)

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("hSDM_poisson_iCAR",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn), ## Number of iterations, burning and samples
               nobs=as.integer(nobs),
               ncell=as.integer(ncell),
               np=as.integer(np),
               Y_vect=as.integer(c(Y)),
               X_vect=as.double(c(X)),
               #= Spatial correlation
               C_vect=as.integer(c(cells)-1), # Cells range is 1,...,ncell in R. Must start at 0 for C. Don't forget the "-1" term. 
               nNeigh=as.integer(c(n.neighbors)),
               Neigh_vect=as.integer(c(neighbors-1)), # Cells range is 1,...,ncell in R. Must start at 0 for C. Don't forget the "-1" term. 
               #= Starting values for M-H
               beta_start=as.double(c(beta.start)),
               rho_start=as.double(c(rho.start)),
               #= Parameters to save
               beta.nonconst=as.double(beta), ## Fixed parameters of the regression
               rho_pred.nonconst=as.double(rho_pred),
               Vrho.nonconst=as.double(Vrho),
               #= Defining priors
               mubeta=as.double(c(mubeta)), Vbeta=as.double(c(Vbeta)),
               priorVrho=as.double(priorVrho),
               shape=as.double(shape), rate=as.double(rate),
               Vrho.max=as.double(Vrho.max),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               prob_p_pred.nonconst=as.double(prob_p_pred), ## Predictive posterior mean
               #= Seed
               seed=as.integer(seed), 
               #= Verbose
               verbose=as.integer(verbose),
               PACKAGE="hSDM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np+2)
  names.fixed <- paste("beta.",colnames(X),sep="")
  colnames(Matrix) <- c(names.fixed,"Vrho","Deviance")
  
  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[14]],ncol=np)
  Matrix[,ncol(Matrix)-1] <- Sample[[16]]
  Matrix[,ncol(Matrix)] <- Sample[[23]]

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)
  
  #= Output
  return (list(mcmc=MCMC,rho.pred=Sample[[15]],prob.p.pred=Sample[[24]]))

}

#===================================================================
# END
#===================================================================
