////////////////////////////////////////////////////////////////////
//
// hSDM.poisson.c
//
////////////////////////////////////////////////////////////////////
//
// Original code by Ghislain Vieilledent, October 2011
// CIRAD UR B&SEF
// ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
//
////////////////////////////////////////////////////////////////////
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2011 Ghislain Vieilledent
// 
////////////////////////////////////////////////////////////////////


// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// R libraries
#include <R.h> // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <Rmath.h> // for dnorm and dbinom distributions
// GSL libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// My own functions
#include "useful.h"


/* ********************************************************************* */
/* dens_par */

struct dens_par {
  /* Data */
  int NOBS;
  int *Y;
  /* Suitability */
  int NP;
  int pos_beta;
  double **X;
  double *mubeta, *Vbeta;
  double *beta_run;  
};


/* ************************************************************ */
/* betadens */

static double betadens (double beta_k, void *dens_data) {
  // Pointer to the structure: d
  struct dens_par *d;
  d=dens_data;
  // Indicating the rank of the parameter of interest
  int k=d->pos_beta;
  // logLikelihood
  double logL=0.0;
  for (int n=0; n<d->NOBS; n++) {
    /* lambda */
    double Xpart_lambda=0.0;
    for (int p=0; p<d->NP; p++) {
      if (p!=k) {
        Xpart_lambda+=d->X[n][p]*d->beta_run[p];
      }
    }
    Xpart_lambda+=d->X[n][k]*beta_k;
    double lambda=exp(Xpart_lambda);
    /* log Likelihood */
    logL+=dpois(d->Y[n],lambda,1);
  }
  // logPosterior=logL+logPrior
  double logP=logL+dnorm(beta_k,d->mubeta[k],sqrt(d->Vbeta[k]),1);
  return logP;
}


/* ************************************************************ */
/* Gibbs sampler function */

void hSDM_poisson (
    
    // Constants and data
    const int *ngibbs, int *nthin, int *nburn, // Number of iterations, burning and samples
    const int *nobs, // Number of observations
    const int *np, // Number of fixed effects for lambda
    const int *Y_vect, // Count
    const double *X_vect, // Suitability covariates
    // Predictions
    const int *npred, // Number of predictions
    const double *X_pred_vect, // Suitability covariates for predictions
    // Starting values for M-H
    const double *beta_start,
    // Parameters to save
    double *beta_vect,
    // Defining priors
    const double *mubeta, double *Vbeta,
    // Diagnostic
    double *Deviance,
    double *lambda_latent, // Latent proba of suitability (length NOBS) 
    double *lambda_pred, // Proba of suitability for predictions (length NPRED)
    // Seeds
    const int *seed,
    // Verbose
    const int *verbose,
    // Save p
    const int *save_p
  
) {
  
  ////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Defining and initializing objects
  
  ////////////////////////////////////////
  // Initialize random number generator //
  gsl_rng *r=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed[0]);
  
  ///////////////////////////
  // Redefining constants //
  const int NGIBBS=ngibbs[0];
  const int NTHIN=nthin[0];
  const int NBURN=nburn[0];
  const int NSAMP=(NGIBBS-NBURN)/NTHIN;
  const int NOBS=nobs[0];
  const int NP=np[0];
  const int NPRED=npred[0];
  
  ///////////////////////////////////
  // Declaring some useful objects //
  double *lambda_run=malloc(NOBS*sizeof(double));
  for (int n=0; n<NOBS; n++) {
    lambda_run[n]=0.0;
  }
  double *lambda_pred_run=malloc(NPRED*sizeof(double));
  for (int m=0; m<NPRED; m++) {
    lambda_pred_run[m]=0.0;
  }
  
  //////////////////////////////////////////////////////////
  // Set up and initialize structure for density function //
  struct dens_par dens_data;
  
  /* Data */
  dens_data.NOBS=NOBS;
  // Y
  dens_data.Y=malloc(NOBS*sizeof(int));
  for (int n=0; n<NOBS; n++) {
    dens_data.Y[n]=Y_vect[n];
  }
  
  /* Suitability process */
  dens_data.NP=NP;
  dens_data.pos_beta=0;
  dens_data.X=malloc(NOBS*sizeof(double*));
  for (int n=0; n<NOBS; n++) {
    dens_data.X[n]=malloc(NP*sizeof(double));
    for (int p=0; p<NP; p++) {
      dens_data.X[n][p]=X_vect[p*NOBS+n];
    }
  }
  dens_data.mubeta=malloc(NP*sizeof(double));
  dens_data.Vbeta=malloc(NP*sizeof(double));
  for (int p=0; p<NP; p++) {
    dens_data.mubeta[p]=mubeta[p];
    dens_data.Vbeta[p]=Vbeta[p];
  }
  dens_data.beta_run=malloc(NP*sizeof(double));
  for (int p=0; p<NP; p++) {
    dens_data.beta_run[p]=beta_start[p];
  }
  
  /* Predictions */
  // X_pred
  double **X_pred=malloc(NPRED*sizeof(double*));
  for (int m=0; m<NPRED; m++) {
    X_pred[m]=malloc(NP*sizeof(double));
    for (int p=0; p<NP; p++) {
      X_pred[m][p]=X_pred_vect[p*NPRED+m];
    }
  }
  
  ////////////////////////////////////////////////////////////
  // Proposal variance and acceptance for adaptive sampling //
  
  // beta
  double *sigmap_beta = malloc(NP*sizeof(double));
  int *nA_beta = malloc(NP*sizeof(int));
  double *Ar_beta = malloc(NP*sizeof(double)); // Acceptance rate 
  for (int p=0; p<NP; p++) {
    nA_beta[p]=0;
    sigmap_beta[p]=1.0;
    Ar_beta[p]=0.0;
  }
  
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  //R_ProcessEvents(); for windows
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Gibbs sampler
  
  for (int g=0;g<NGIBBS;g++) {
    
    ////////////////////////////////////////////////
    // beta
    
    for (int p=0; p<NP; p++) {
      dens_data.pos_beta=p; // Specifying the rank of the parameter of interest
      double x_now=dens_data.beta_run[p];
      double x_prop=x_now+gsl_ran_gaussian_ziggurat(r, sigmap_beta[p]);
      double p_now=betadens(x_now, &dens_data);
      double p_prop=betadens(x_prop, &dens_data);
      double ratio=exp(p_prop-p_now); // ratio
      double z=gsl_rng_uniform(r);	
      // Actualization
      if (z < ratio) {
        dens_data.beta_run[p]=x_prop;
        nA_beta[p]++;
      }
    }
    
    
    //////////////////////////////////////////////////
    // Deviance
    
    // logLikelihood
    double logL=0.0;
    for (int n=0; n<NOBS; n++) {
      /* lambda */
      double Xpart_lambda=0.0;
      for (int p=0; p<NP; p++) {
        Xpart_lambda+=dens_data.X[n][p]*dens_data.beta_run[p];
      }
      lambda_run[n]=exp(Xpart_lambda);
      /* log Likelihood */
      logL+=dpois(dens_data.Y[n],lambda_run[n],1);
    }
    
    // Deviance
    double Deviance_run=-2*logL;
    
    
    //////////////////////////////////////////////////
    // Predictions
    for (int m=0; m<NPRED; m++) {
      /* lambda_pred_run */
      double Xpart_lambda_pred=0.0;
      for (int p=0; p<NP; p++) {
        Xpart_lambda_pred+=X_pred[m][p]*dens_data.beta_run[p];
      }
      lambda_pred_run[m]=exp(Xpart_lambda_pred);
    }
    
    
    //////////////////////////////////////////////////
    // Output
    if (((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for (int p=0; p<NP; p++) {
        beta_vect[p*NSAMP+(isamp-1)]=dens_data.beta_run[p];
      }
      Deviance[isamp-1]=Deviance_run;
      for (int n=0; n<NOBS; n++) {
        lambda_latent[n]+=lambda_run[n]/NSAMP; // We compute the mean of NSAMP values
      }
      // lambda
      if (save_p[0]==0) { // We compute the mean of NSAMP values
        for (int m=0; m<NPRED; m++) {
          lambda_pred[m]+=lambda_pred_run[m]/NSAMP; 
        }
      }
      if (save_p[0]==1) { // The NSAMP sampled values for lambda are saved
        for (int m=0; m<NPRED; m++) {
          lambda_pred[m*NSAMP+(isamp-1)]=lambda_pred_run[m]; 
        }
      }
    }
    
    
    ///////////////////////////////////////////////////////
    // Adaptive sampling (on the burnin period)
    const double ropt=0.44;
    int DIV=0;
    if (NGIBBS >=1000) DIV=100;
    else DIV=NGIBBS/10;
    /* During the burnin period */
    if((g+1)%DIV==0 && (g+1)<=NBURN){
      for (int p=0; p<NP; p++) {
        Ar_beta[p]=((double) nA_beta[p])/DIV;
        if(Ar_beta[p]>=ropt) sigmap_beta[p]=sigmap_beta[p]*(2-(1-Ar_beta[p])/(1-ropt));
        else sigmap_beta[p]=sigmap_beta[p]/(2-Ar_beta[p]/ropt);
        nA_beta[p]=0.0; // We reinitialize the number of acceptance to zero
      }
    }
    /* After the burnin period */
    if((g+1)%DIV==0 && (g+1)>NBURN){
      for (int p=0; p<NP; p++) {
        Ar_beta[p]=((double) nA_beta[p])/DIV;
        nA_beta[p]=0.0; // We reinitialize the number of acceptance to zero
      }
    }
    
    
    //////////////////////////////////////////////////
    // Progress bar
    double Perc=100*(g+1)/(NGIBBS);
    if (((g+1)%(NGIBBS/100))==0 && verbose[0]==1) {  
      Rprintf("*");
      R_FlushConsole();
      //R_ProcessEvents(); for windows
      if(((g+1)%(NGIBBS/10))==0){
        double mAr_beta=0; // Mean acceptance rate
        for (int p=0; p<NP; p++) {
          mAr_beta+=Ar_beta[p]/NP;
        }
        Rprintf(":%.1f%%, mean accept. rates= beta:%.3f\n",Perc,mAr_beta);
        R_FlushConsole();
        //R_ProcessEvents(); for windows
      }
    } 
    
    
    //////////////////////////////////////////////////
    // User interrupt
    R_CheckUserInterrupt(); // allow user interrupt 	    
    
  } // Gibbs sampler
  
  ///////////////
  // Delete memory allocation (see malloc())
  /* Data */
  free(dens_data.Y);
  /* Suitability */
  for (int n=0; n<NOBS; n++) {
    free(dens_data.X[n]);
  }
  free(dens_data.X);
  free(dens_data.mubeta);
  free(dens_data.Vbeta);
  free(dens_data.beta_run);
  free(lambda_run);
  /* Predictions */
  for (int m=0; m<NPRED; m++) {
    free(X_pred[m]);
  }
  free(X_pred);
  free(lambda_pred_run);
  /* Adaptive MH */
  free(sigmap_beta);
  free(nA_beta);
  free(Ar_beta);
  /* Random seed */
  gsl_rng_free(r);
  
} // end hSDM function


////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
