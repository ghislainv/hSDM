####################################################################
##
## hSDM.hierarchical.poisson.R
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


hSDM.hierarchical.poisson <- function (counts,
                                       suitability,
                                       cells,
                                       n.neighbors,
                                       neighbors,
                                       alteration,
                                       observability,    
                                       data, burnin=5000,
                                       mcmc=10000, thin=10,
                                       beta.start,
                                       gamma.start,
                                       Vrho.start,
                                       mubeta=0, Vbeta=1.0E6,
                                       mugamma=0, Vgamma=1.0E6,
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
  #= Alteration
  U <- alteration
  #= Suitability
  mf.suit <- model.frame(formula=suitability,data=data)
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  #= Observability
  mf.obs <- model.frame(formula=observability,data=data)
  W <- model.matrix(attr(mf.obs,"terms"),data=mf.obs)
  #= Spatial correlation
  cells <- cells
  ncell <- length(unique(cells))
  n.neighbors <- n.neighbors
  neighbors <- neighbors
  #= Model parameters
  np <- ncol(X)
  nq <- ncol(W)
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin

  #========== 
  # Check data
  #==========
  check.Y.poisson(Y)
  check.X(X,nobs)
  check.U(U,nobs)
  check.W(W,nobs)
  check.cells(cells,nobs)
  check.neighbors(n.neighbors,ncell,neighbors)

  #========
  # Initial starting values for M-H
  #========
  beta.start <- form.beta.start(beta.start,np)
  gamma.start <- form.gamma.start(gamma.start,nq)
  rho.start <- rep(0,ncell) # Starting values for spatial random effects set to zero.
  Vrho.start <- check.Vrho.start(Vrho.start)
  
  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mubeta,np)
  Vbeta <- check.Vbeta(Vbeta,np)
  mugamma <- check.mugamma(mugamma,nq)
  Vgamma <- check.Vgamma(Vgamma,nq)
  check.ig.prior(shape,rate)
  Vrho.max <- check.Vrho.max(Vrho.max)
  priorVrho <- form.priorVrho(priorVrho)

  #========
  # Parameters to save
  #========
  beta <- rep(beta.start,nsamp)
  gamma <- rep(gamma.start,nsamp)
  rho_pred <- rho.start
  Vrho <- rep(Vrho.start,nsamp)
  prob_p_pred <- rep(0,nobs)
  prob_q_pred <- rep(0,nobs)
  Deviance <- rep(0,nsamp)

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("hSDM_hierarchical_poisson",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn), ## Number of iterations, burning and samples
               nobs=as.integer(nobs),
               ncell=as.integer(ncell),
               np=as.integer(np),
               nq=as.integer(nq),
               Y_vect=as.integer(c(Y)),
               X_vect=as.double(c(X)),
               W_vect=as.double(c(W)),
               U_vect=as.double(c(U)),
               #= Spatial correlation
               C_vect=as.integer(c(cells)-1), # Cells range is 1,...,ncell in R. Must start at 0 for C. Don't forget the "-1" term. 
               nNeigh=as.integer(c(n.neighbors)),
               Neigh_vect=as.integer(c(neighbors-1)), # Cells range is 1,...,ncell in R. Must start at 0 for C. Don't forget the "-1" term.
               #= Starting values for M-H
               beta_start=as.double(c(beta.start)),
               gamma_start=as.double(c(gamma.start)),
               rho_start=as.double(c(rho.start)),
               #= Parameters to save
               beta.nonconst=as.double(beta), ## Fixed parameters of the regression
               gamma.nonconst=as.double(gamma),
               rho_pred.nonconst=as.double(rho_pred), 
               Vrho.nonconst=as.double(Vrho), 
               #= Defining priors
               mubeta=as.double(c(mubeta)), Vbeta=as.double(c(Vbeta)),
               mugamma=as.double(c(mugamma)), Vgamma=as.double(c(Vgamma)),
               priorVrho=as.double(priorVrho),
               shape=as.double(shape), rate=as.double(rate),
               Vrho.max=as.double(Vrho.max),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               prob_p_pred.nonconst=as.double(prob_p_pred), ## Predictive posterior mean
               prob_q_pred.nonconst=as.double(prob_q_pred), ## Predictive posterior mean
               #= Seed
               seed=as.integer(seed), 
               #= Verbose
               verbose=as.integer(verbose),
               PACKAGE="hSDM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np+nq+2)
  names.fixed <- c(paste("beta.",colnames(X),sep=""),paste("gamma.",colnames(W),sep=""))
  colnames(Matrix) <- c(names.fixed,"Vrho","Deviance")
  
  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[18]],ncol=np)
  Matrix[,c((np+1):(np+nq))] <- matrix(Sample[[19]],ncol=nq)
  Matrix[,ncol(Matrix)-1] <- Sample[[21]]
  Matrix[,ncol(Matrix)] <- Sample[[30]]

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)
  
  #= Output
  return (list(mcmc=MCMC,rho.pred=Sample[[20]],prob.p.pred=Sample[[31]],prob.q.pred=Sample[[32]]))

}

#===================================================================
# END
#===================================================================
