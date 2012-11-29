###################################################################
##
## hSDM.binomial.R
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
## - G. Vieilledent, on November 15th 2012
##
####################################################################


hSDM.binomial <- function (presences, trials,
                           suitability,
                           data, burnin=5000,
                           mcmc=10000, thin=10, 
                           beta.start,
                           mubeta=0, Vbeta=1.0E6,
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
  Y <- presences
  nobs <- length(Y)
  T <- trials
  #= Suitability
  mf.suit <- model.frame(formula=suitability,data=data)
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  #= Model parameters
  np <- ncol(X)
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin

  #========== 
  # Check data
  #==========
  check.T.binomial(T,nobs)
  check.Y.binomial(Y,T)
  check.X(X,nobs)
  
  #========
  # Initial starting values for M-H
  #========
  beta.start <- form.beta.start(beta.start,np)

  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mubeta,np)
  Vbeta <- check.Vbeta(Vbeta,np)

  #========
  # Parameters to save
  #========
  beta <- rep(beta.start,nsamp)
  prob_p_pred <- rep(0,nobs)
  Deviance <- rep(0,nsamp)

  #========
  # call C++ code to draw sample
  #========
  Sample <- .C("hSDM_binomial",
               #= Constants and data
               ngibbs=as.integer(ngibbs), nthin=as.integer(nthin), nburn=as.integer(nburn), ## Number of iterations, burning and samples
               nobs=as.integer(nobs),
               np=as.integer(np),
               Y_vect=as.integer(c(Y)),
               T_vect=as.integer(c(T)),
               X_vect=as.double(c(X)),
               #= Starting values for M-H
               beta_start=as.double(c(beta.start)),
               #= Parameters to save
               beta.nonconst=as.double(beta), ## Fixed parameters of the regression
               #= Defining priors
               mubeta=as.double(c(mubeta)), Vbeta=as.double(c(Vbeta)),
               #= Diagnostic
               Deviance.nonconst=as.double(Deviance),
               prob_p_pred.nonconst=as.double(prob_p_pred), ## Predictive posterior mean
               #= Seed
               seed=as.integer(seed), 
               #= Verbose
               verbose=as.integer(verbose),
               PACKAGE="hSDM")
 
  #= Matrix of MCMC samples
  Matrix <- matrix(NA,nrow=nsamp,ncol=np+1)
  names.fixed <- paste("beta.",colnames(X),sep="")
  colnames(Matrix) <- c(names.fixed,"Deviance")
  
  #= Filling-in the matrix
  Matrix[,c(1:np)] <- matrix(Sample[[10]],ncol=np)
  Matrix[,ncol(Matrix)] <- Sample[[13]]

  #= Transform Sample list in an MCMC object
  MCMC <- mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)
  
  #= Output
  return (list(mcmc=MCMC,prob.p.pred=Sample[[14]]))

}

#===================================================================
# END
#===================================================================
