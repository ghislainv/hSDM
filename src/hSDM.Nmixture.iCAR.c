////////////////////////////////////////////////////////////////////
//
// hSDM.Nmixture.iCAR.c
//
////////////////////////////////////////////////////////////////////
//
// Original code by Ghislain Vieilledent, November 2013
// CIRAD UR B&SEF
// ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
//
////////////////////////////////////////////////////////////////////
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 3. See the package LICENSE file for more
// information.
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
  /* Sites */
  int NSITE;
  int *IdSiteforObs;
  int *nObsSite;
  int **ListObsBySite;
  /* Latent variable */
  int *N_run;
  int pos_N;
  /* Spatial cells */
  int *IdCellforSite;
  int *nSiteCell;
  int **ListSiteByCell;
  /* Spatial correlation */
  int *nNeigh;
  int **Neigh;
  int pos_rho;
  double *rho_run;
  double shape, rate;
  double Vrho_run;
  /* Suitability */
  int NP;
  int pos_beta;
  double **X;
  double *mubeta, *Vbeta;
  double *beta_run;
  /* Observability */
  int NQ;
  int pos_gamma;
  double **W;
  double *mugamma, *Vgamma;
  double *gamma_run;    
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
  for (int i=0; i<d->NSITE; i++) {
    /* lambda */
    double Xpart_lambda=0.0;
    for (int p=0; p<d->NP; p++) {
      if (p!=k) {
        Xpart_lambda+=d->X[i][p]*d->beta_run[p];
      }
    }
    Xpart_lambda+=d->X[i][k]*beta_k;
    double lambda=exp(Xpart_lambda+d->rho_run[d->IdCellforSite[i]]);
    /* log Likelihood */
    logL+=dpois(d->N_run[i],lambda,1);
  }
  // logPosterior=logL+logPrior
  double logP=logL+dnorm(beta_k,d->mubeta[k],sqrt(d->Vbeta[k]),1);
  return logP;
}

/* ************************************************************ */
/* gammadens */

static double gammadens (double gamma_k, void *dens_data) {
  // Pointer to the structure: d 
  struct dens_par *d;
  d=dens_data;
  // Indicating the rank of the parameter of interest
  int k=d->pos_gamma; //
  // logLikelihood
  double logL=0.0;
  for (int n=0; n<d->NOBS; n++) {
    /* delta */
    double logit_delta=0.0;
    for (int q=0; q<d->NQ; q++) {
      if (q!=k) {
        logit_delta+=d->W[n][q]*d->gamma_run[q];
      }
    }
    logit_delta+=d->W[n][k]*gamma_k;
    double delta=invlogit(logit_delta);
    /* log Likelihood */
    logL+=dbinom(d->Y[n],d->N_run[d->IdSiteforObs[n]],delta,1);
  }
  // logPosterior=logL+logPrior
  double logP=logL+dnorm(gamma_k,d->mugamma[k],sqrt(d->Vgamma[k]),1); 
  return logP;
}

/* ************************************************************ */
/* Ndens */ 

static double Ndens (int N_i, void *dens_data) {
  // Pointer to the structure: d 
  struct dens_par *d;
  d=dens_data;
  // Indicating the rank of the site
  int i=d->pos_N; //
  // logLikelihood
  double logL=0;
  for (int m=0; m<d->nObsSite[i]; m++) {
    int wo=d->ListObsBySite[i][m]; // which observation
    /* delta */
    double logit_delta=0.0;
    for (int q=0; q<d->NQ; q++) {
      logit_delta+=d->W[wo][q]*d->gamma_run[q];
    }
    double delta=invlogit(logit_delta);
    /* log Likelihood */
    logL+=dbinom(d->Y[wo],N_i,delta,1);
  }
  // logPosterior=logL+logPrior
  /* lambda */
  double Xpart_lambda=0.0;
  for (int p=0; p<d->NP; p++) {
    Xpart_lambda+=d->X[i][p]*d->beta_run[p];
  }
  double lambda=exp(Xpart_lambda+d->rho_run[d->IdCellforSite[i]]);
  double logP=logL+dpois(N_i,lambda,1); 
  return logP;
}

/* ************************************************************ */
/* rhodens_visited */

static double rhodens_visited (double rho_j, void *dens_data) {
  // Pointer to the structure: d 
  struct dens_par *d;
  d=dens_data;
  // Indicating the rank of the cell
  int j=d->pos_rho; //
  // logLikelihood
  double logL=0;
  // Loop on sites in the cell
  for (int i=0; i<d->nSiteCell[j]; i++) {
    int ws=d->ListSiteByCell[j][i]; // which site
    /* lambda */
    double Xpart_lambda=0.0;
    for (int p=0; p<d->NP; p++) {
      Xpart_lambda+=d->X[ws][p]*d->beta_run[p];
    }
    double lambda=exp(Xpart_lambda+rho_j);
    /* log Likelihood */
    logL+=dpois(d->N_run[ws],lambda,1);
  }
  // logPosterior=logL+logPrior
  int nNeighbors=d->nNeigh[j];
  double sumNeighbors=0.0;
  for (int m=0;m<nNeighbors;m++) {
    sumNeighbors+=d->rho_run[d->Neigh[j][m]];
  }
  double meanNeighbors=sumNeighbors/nNeighbors;
  double logP=logL+dnorm(rho_j,meanNeighbors,sqrt(d->Vrho_run/nNeighbors),1); 
  return logP;
}

/* ************************************************************ */
/* rhodens_unvisited */

static double rhodens_unvisited (const gsl_rng *r, void *dens_data) {
  // Pointer to the structure: d 
  struct dens_par *d;
  d=dens_data;
  // Indicating the rank of the site
  int j=d->pos_rho; //
  // Draw directly in the posterior distribution
  int nNeighbors=d->nNeigh[j];
  double sumNeighbors=0.0;
  for (int m=0;m<nNeighbors;m++) {
    sumNeighbors+=d->rho_run[d->Neigh[j][m]];
  }
  double meanNeighbors=sumNeighbors/nNeighbors;
  double sample=meanNeighbors+gsl_ran_gaussian_ziggurat(r,sqrt(d->Vrho_run/nNeighbors));
  return sample;
}

/* ************************************************************ */
/* Gibbs sampler function */

void hSDM_Nmixture_iCAR (
    
    // Constants and data
    const int *ngibbs, int *nthin, int *nburn, // Number of iterations, burning and samples
    const int *nobs, int *nsite, int *ncell, // Number of observations, sites and cells
    const int *np, // Number of fixed effects for lambda
    const int *nq, // Number of fixed effects for delta
    const int *Y_vect, // Number of successes (presences)
    const double *W_vect, // Observability covariates (nobs x nq)
    const double *X_vect, // Suitability covariates (nsite x np)
    // Sites
    const int *S_vect, // Site Id (nobs)
    // Spatial correlation
    const int *C_vect, // Cell Id (nsite)
    const int *nNeigh, // Number of neighbors for each cell
    const int *Neigh_vect, // Vector of neighbors sorted by cell
    // Predictions
    const int *npred, // Number of predictions
    const double *X_pred_vect, // Suitability covariates for predictions (npred x np)
    const int *C_pred_vect, // Cell Id for predictions (npred)
    // Starting values for M-H
    const double *beta_start,
    const double *gamma_start,
    const double *rho_start,
    const int *N_start,
    // Parameters
    double *beta_vect,
    double *gamma_vect,
    double *rho_pred,
    double *Vrho,
    int *N_pred,
    // Defining priors
    const double *mubeta, double *Vbeta,
    const double *mugamma, double *Vgamma,
    const double *priorVrho,
    const double *shape, double *rate,
    const double *Vrho_max,
    // Diagnostic
    double *Deviance,
    double *lambda_latent, // Latent proba of suitability (length NSITE)
    double *delta_latent, // Latent proba of observability (length NOBS)
    double *lambda_pred, // Proba of suitability for predictions (length NPRED)
    // Seeds
    const int *seed,
    // Verbose
    const int *verbose,
    // Save rho, p and N
    const int *save_rho,
    const int *save_p,
    const int *save_N
  
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
  const int NSITE=nsite[0];
  const int NCELL=ncell[0];
  const int NP=np[0];
  const int NQ=nq[0];
  const int NPRED=npred[0];
  
  ///////////////////////////////////
  // Declaring some useful objects //
  double *lambda_run=malloc(NSITE*sizeof(double));
  for (int i=0; i<NSITE; i++) {
    lambda_run[i]=0.0;
  }
  double *delta_run=malloc(NOBS*sizeof(double));
  for (int n=0; n<NOBS; n++) {
    delta_run[n]=0.0;
  }
  double *lambda_pred_run=malloc(NPRED*sizeof(double));
  for (int m=0; m<NPRED; m++) {
    lambda_pred_run[m]=0.0;
  }
  double *N_pred_double=malloc(NSITE*sizeof(double));
  for (int i=0; i<NSITE; i++) {
    N_pred_double[i]=0.0;
  }
  
  //////////////////////////////////////////////////////////
  // Set up and initialize structure for density function //
  struct dens_par dens_data;
  
  /* Data */
  dens_data.NOBS=NOBS;
  dens_data.NSITE=NSITE;
  // Y
  dens_data.Y=malloc(NOBS*sizeof(int));
  for (int n=0; n<NOBS; n++) {
    dens_data.Y[n]=Y_vect[n];
  }
  
  /* Sites */
  // IdSiteforObs
  dens_data.IdSiteforObs=malloc(NOBS*sizeof(int));
  for (int n=0; n<NOBS; n++) {
    dens_data.IdSiteforObs[n]=S_vect[n];
  }
  // nObsSite
  dens_data.nObsSite=malloc(NSITE*sizeof(int));
  for (int i=0; i<NSITE; i++) {
    dens_data.nObsSite[i]=0;
    for (int n=0; n<NOBS; n++) {
      if (dens_data.IdSiteforObs[n]==i) {
        dens_data.nObsSite[i]++;
      }
    }
  }
  // ListObsBySite
  dens_data.ListObsBySite=malloc(NSITE*sizeof(int*));
  for (int i=0; i<NSITE; i++) {
    dens_data.ListObsBySite[i]=malloc(dens_data.nObsSite[i]*sizeof(int));
    int repSite=0;
    for (int n=0; n<NOBS; n++) {
      if (dens_data.IdSiteforObs[n]==i) {
        dens_data.ListObsBySite[i][repSite]=n;
        repSite++;
      }
    }
  }
  
  /* Latent variable */
  // N_run
  dens_data.N_run=malloc(NSITE*sizeof(int));
  for (int i=0; i<NSITE; i++) {
    dens_data.N_run[i]=N_start[i];
  }
  dens_data.pos_N=0;
  
  /* Spatial correlation */
  // IdCellforSite
  dens_data.IdCellforSite=malloc(NSITE*sizeof(int));
  for (int i=0; i<NSITE; i++) {
    dens_data.IdCellforSite[i]=C_vect[i];
  }
  // nSiteCell
  dens_data.nSiteCell=malloc(NCELL*sizeof(int));
  for (int j=0; j<NCELL; j++) {
    dens_data.nSiteCell[j]=0;
    for (int i=0; i<NSITE; i++) {
      if (dens_data.IdCellforSite[i]==j) {
        dens_data.nSiteCell[j]++;
      }
    }
  }
  // ListSiteByCell
  dens_data.ListSiteByCell=malloc(NCELL*sizeof(int*));
  for (int j=0; j<NCELL; j++) {
    dens_data.ListSiteByCell[j]=malloc(dens_data.nSiteCell[j]*sizeof(int));
    int repCell=0;
    for (int i=0; i<NSITE; i++) {
      if (dens_data.IdCellforSite[i]==j) {
        dens_data.ListSiteByCell[j][repCell]=i;
        repCell++;
      }
    }
  }
  // Number of neighbors by cell
  dens_data.nNeigh=malloc(NCELL*sizeof(int));
  for (int j=0; j<NCELL; j++) {
    dens_data.nNeigh[j]=nNeigh[j];
  }
  // Neighbor identifiers by cell
  int posNeigh=0;
  dens_data.Neigh=malloc(NCELL*sizeof(int*));
  for (int j=0; j<NCELL; j++) {
    dens_data.Neigh[j]=malloc(nNeigh[j]*sizeof(int));
    for (int m=0; m<nNeigh[j]; m++) {
      dens_data.Neigh[j][m]=Neigh_vect[posNeigh+m];
    }
    posNeigh+=nNeigh[j];
  }
  dens_data.pos_rho=0;
  dens_data.rho_run=malloc(NCELL*sizeof(double));
  for (int j=0; j<NCELL; j++) {
    dens_data.rho_run[j]=rho_start[j];
  }
  dens_data.shape=shape[0];
  dens_data.rate=rate[0];
  dens_data.Vrho_run=Vrho[0];
  
  /* Suitability process */
  dens_data.NP=NP;
  dens_data.pos_beta=0;
  dens_data.X=malloc(NSITE*sizeof(double*));
  for (int i=0; i<NSITE; i++) {
    dens_data.X[i]=malloc(NP*sizeof(double));
    for (int p=0; p<NP; p++) {
      dens_data.X[i][p]=X_vect[p*NSITE+i];
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
  
  /* Observability process */
  dens_data.NQ=NQ;
  dens_data.pos_gamma=0;
  dens_data.W=malloc(NOBS*sizeof(double*));
  for (int n=0; n<NOBS; n++) {
    dens_data.W[n]=malloc(NQ*sizeof(double));
    for (int q=0; q<NQ; q++) {
      dens_data.W[n][q]=W_vect[q*NOBS+n];
    }
  }
  dens_data.mugamma=malloc(NQ*sizeof(double));
  dens_data.Vgamma=malloc(NQ*sizeof(double));
  for (int q=0; q<NQ; q++) {
    dens_data.mugamma[q]=mugamma[q];
    dens_data.Vgamma[q]=Vgamma[q];
  }
  dens_data.gamma_run=malloc(NQ*sizeof(double));
  for (int q=0; q<NQ; q++) {
    dens_data.gamma_run[q]=gamma_start[q];
  }
  
  /* Visited cell or not */
  int *viscell = malloc(NCELL*sizeof(int));
  for (int j=0; j<NCELL; j++) {
    viscell[j]=0;
  }
  for (int i=0; i<NSITE; i++) {
    viscell[dens_data.IdCellforSite[i]]++;
  }
  int NVISCELL=0;
  for (int j=0; j<NCELL; j++) {
    if (viscell[j]>0) {
      NVISCELL++;
    }
  }
  
  /* Predictions */
  // IdCell_pred
  int *IdCell_pred=malloc(NPRED*sizeof(int));
  for (int m=0; m<NPRED; m++) {
    IdCell_pred[m]=C_pred_vect[m];
  }
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
  
  // gamma
  double *sigmap_gamma = malloc(NQ*sizeof(double));
  int *nA_gamma = malloc(NQ*sizeof(int));
  double *Ar_gamma = malloc(NQ*sizeof(double)); // Acceptance rate 
  for (int q=0; q<NQ; q++) {
    nA_gamma[q]=0;
    sigmap_gamma[q]=1.0;
    Ar_gamma[q]=0.0;
  }
  
  // rho
  double *sigmap_rho = malloc(NCELL*sizeof(double));
  int *nA_rho = malloc(NCELL*sizeof(int));
  double *Ar_rho = malloc(NCELL*sizeof(double)); // Acceptance rate 
  for (int i=0; i<NCELL; i++) {
    nA_rho[i]=0;
    sigmap_rho[i]=1.0;
    Ar_rho[i]=0.0;
  }
  
  // N
  int *nA_N = malloc(NSITE*sizeof(int));
  double *Ar_N = malloc(NSITE*sizeof(double)); // Acceptance rate 
  for (int i=0; i<NSITE; i++) {
    nA_N[i]=0;
    Ar_N[i]=0.0;
  }
  
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  //R_ProcessEvents(); for windows
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Gibbs sampler
  
  for (int g=0; g<NGIBBS; g++) {
    
    
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
    
    
    ////////////////////////////////////////////////
    // gamma
    
    for (int q=0; q<NQ; q++) {
      dens_data.pos_gamma=q; // Specifying the rank of the parameter of interest
      double x_now=dens_data.gamma_run[q];
      double x_prop=x_now+gsl_ran_gaussian_ziggurat(r, sigmap_gamma[q]);
      double p_now=gammadens(x_now, &dens_data);
      double p_prop=gammadens(x_prop, &dens_data);
      double ratio=exp(p_prop-p_now); // ratio
      double z=gsl_rng_uniform(r);
      // Actualization
      if (z < ratio) {
        dens_data.gamma_run[q]=x_prop;
        nA_gamma[q]++;
      }
    }
    
    
    ////////////////////////////////////////////////
    // rho
    
    /* Sampling rho_run[j] */
    for (int j=0; j<NCELL; j++) {
      dens_data.pos_rho=j; // Specifying the rank of the parameter of interest
      if (viscell[j]>0) {
        double x_now=dens_data.rho_run[j];
        double x_prop=x_now+gsl_ran_gaussian_ziggurat(r, sigmap_rho[j]);
        double p_now=rhodens_visited(x_now, &dens_data);
        double p_prop=rhodens_visited(x_prop, &dens_data);
        double ratio=exp(p_prop-p_now); // ratio
        double z=gsl_rng_uniform(r);
        // Actualization
        if (z < ratio) {
          dens_data.rho_run[j]=x_prop;
          nA_rho[j]++;
        }
      }
      else {
        dens_data.rho_run[j]=rhodens_unvisited(r, &dens_data);
      }
    }
    
    /* Centering rho_run[j] */
    double rho_sum=0.0;
    for (int j=0; j<NCELL; j++) {
      rho_sum+=dens_data.rho_run[j];
    }
    double rho_bar=rho_sum/NCELL;
    for (int j=0; j<NCELL; j++) {
      dens_data.rho_run[j]=dens_data.rho_run[j]-rho_bar;
    }
    
    
    ////////////////////////////////////////////////
    // Vrho
    
    if (priorVrho[0]>0.0) { // fixed value for Vrho
      dens_data.Vrho_run=priorVrho[0];
    }
    else {
      double Sum=0.0;
      for (int j=0; j<NCELL; j++) {
        double Sum_neigh=0.0;
        double nNeigh=dens_data.nNeigh[j];
        double rho_run=dens_data.rho_run[j];
        for (int m=0; m<nNeigh; m++) {
          Sum_neigh += dens_data.rho_run[dens_data.Neigh[j][m]];
        }
        Sum += rho_run*(nNeigh*rho_run-Sum_neigh);
      }
      if (priorVrho[0]==-1.0) { // prior = 1/Gamma(shape,rate)
        double Shape=shape[0]+0.5*(NCELL-1);
        double Rate=rate[0]+0.5*Sum;
        dens_data.Vrho_run=Rate/gsl_ran_gamma(r, Shape, 1.0);
      }
      if (priorVrho[0]==-2.0) { // prior = Uniform(0,Vrho_max)
        double Shape=0.5*NCELL-1;
        double Rate=0.5*Sum;
        dens_data.Vrho_run=1/myrtgamma_left_gsl(r, Shape, Rate, 1/Vrho_max[0]);
      }
    }
    
    
    ////////////////////////////////////////////////
    // N
    
    for (int i=0; i<NSITE; i++) {
      dens_data.pos_N=i; // Specifying the rank of the parameter of interest
      int x_now=dens_data.N_run[i];
      if (x_now==0) {
        double s=gsl_rng_uniform(r);
        if (s < 0.5) {
          //dens_data.N_run[i]=x_now;
          dens_data.N_run[i]=0;
        }
        else {
          // Proposal
          //int x_prop=x_now+1;
          int x_prop=1;
          // Ratio
          double p_now=Ndens(x_now, &dens_data);
          double p_prop=Ndens(x_prop, &dens_data);
          double ratio=exp(p_prop-p_now);
          // Actualization
          double z=gsl_rng_uniform(r);
          if (z < ratio) {
            dens_data.N_run[i]=x_prop;
            nA_N[i]++;
          }
        }
      }
      else {
        // Proposal
        double s=gsl_rng_uniform(r);
        int x_prop=0;
        if (s < 0.5) x_prop=x_now-1;
        else x_prop=x_now+1;
        // Ratio
        double p_now=Ndens(x_now, &dens_data);
        double p_prop=Ndens(x_prop, &dens_data);
        double ratio=exp(p_prop-p_now);
        // Actualization
        double z=gsl_rng_uniform(r);
        if (z < ratio) {
          dens_data.N_run[i]=x_prop;
          nA_N[i]++;
        }
      }
    }
    
    
    //////////////////////////////////////////////////
    // Deviance
    
    // logLikelihood
    double logL1=0.0;
    for (int n=0; n<NOBS; n++) {
      int ws=dens_data.IdSiteforObs[n]; // which site
      /* delta */
      double logit_delta=0.0;
      for (int q=0; q<NQ; q++) {
        logit_delta+=dens_data.W[n][q]*dens_data.gamma_run[q];
      }
      delta_run[n]=invlogit(logit_delta);
      /* log Likelihood */
      logL1+=dbinom(dens_data.Y[n],dens_data.N_run[ws],delta_run[n],1);
    }
    double logL2=0.0;
    for (int i=0; i<NSITE; i++) {
      int wc=dens_data.IdCellforSite[i]; // which cell
      /* lambda */
      double Xpart_lambda=0.0;
      for (int p=0; p<NP; p++) {
        Xpart_lambda+=dens_data.X[i][p]*dens_data.beta_run[p];
      }
      lambda_run[i]=exp(Xpart_lambda+dens_data.rho_run[wc]);
      logL2+=dpois(dens_data.N_run[i],lambda_run[i],1);
    }
    double logL=logL1+logL2;
    
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
      lambda_pred_run[m]=exp(Xpart_lambda_pred+dens_data.rho_run[IdCell_pred[m]]);
    }
    
    
    //////////////////////////////////////////////////
    // Output
    if (((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      // beta
      for (int p=0; p<NP; p++) {
        beta_vect[p*NSAMP+(isamp-1)]=dens_data.beta_run[p];
      }
      // gamma
      for (int q=0; q<NQ; q++) {
        gamma_vect[q*NSAMP+(isamp-1)]=dens_data.gamma_run[q];
      }
      // Deviance
      Deviance[isamp-1]=Deviance_run;
      for (int i=0; i<NSITE; i++) {
        lambda_latent[i]+=lambda_run[i]/NSAMP; // We compute the mean of NSAMP values
      }
      for (int n=0; n<NOBS; n++) {
        delta_latent[n]+=delta_run[n]/NSAMP; // We compute the mean of NSAMP values
      }
      // rho
      if (save_rho[0]==0) { // We compute the mean of NSAMP values
        for (int j=0; j<NCELL; j++) {
          rho_pred[j]+=dens_data.rho_run[j]/NSAMP; 
        }
      }
      if (save_rho[0]==1) { // The NSAMP sampled values for rhos are saved
        for (int j=0; j<NCELL; j++) {
          rho_pred[j*NSAMP+(isamp-1)]=dens_data.rho_run[j]; 
        }
      }
      // lambda_pred
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
      // Vrho
      Vrho[isamp-1]=dens_data.Vrho_run;
      // N
      if (save_N[0]==0) { // We compute the mean of NSAMP values
        for (int i=0; i<NSITE; i++) {
          N_pred_double[i]+= ((double) dens_data.N_run[i])/NSAMP;
        }
      }
      if (save_N[0]==1) { // The NSAMP sampled values for lambda are saved
        for (int i=0; i<NSITE; i++) {
          N_pred[i*NSAMP+(isamp-1)]=dens_data.N_run[i];
        }
      }
    }
    
    
    ///////////////////////////////////////////////////////
    // Adaptive sampling (on the burnin period)
    const double ropt=0.234;
    int DIV=0;
    if (NGIBBS >=1000) DIV=100;
    else DIV=NGIBBS/10;
    /* During the burnin period */
    if ((g+1)%DIV==0 && (g+1)<=NBURN) {
      // beta
      for (int p=0; p<NP; p++) {
        Ar_beta[p]=((double) nA_beta[p])/DIV;
        if(Ar_beta[p]>=ropt) sigmap_beta[p]=sigmap_beta[p]*(2-(1-Ar_beta[p])/(1-ropt));
        else sigmap_beta[p]=sigmap_beta[p]/(2-Ar_beta[p]/ropt);
        nA_beta[p]=0.0; // We reinitialize the number of acceptance to zero
      }
      // gamma
      for (int q=0; q<NQ; q++) {
        Ar_gamma[q]=((double) nA_gamma[q])/DIV;
        if(Ar_gamma[q]>=ropt) sigmap_gamma[q]=sigmap_gamma[q]*(2-(1-Ar_gamma[q])/(1-ropt));
        else sigmap_gamma[q]=sigmap_gamma[q]/(2-Ar_gamma[q]/ropt);
        nA_gamma[q]=0.0; // We reinitialize the number of acceptance to zero
      }
      // rho
      for (int j=0; j<NCELL; j++) {
        if (viscell[j]>0) {
          Ar_rho[j]=((double) nA_rho[j])/DIV;
          if(Ar_rho[j]>=ropt) sigmap_rho[j]=sigmap_rho[j]*(2-(1-Ar_rho[j])/(1-ropt));
          else sigmap_rho[j]=sigmap_rho[j]/(2-Ar_rho[j]/ropt);
          nA_rho[j]=0.0; // We reinitialize the number of acceptance to zero
        }
      }
      // N
      for (int i=0; i<NSITE; i++) {
        Ar_N[i]=((double) nA_N[i])/DIV;
        nA_N[i]=0.0; // We reinitialize the number of acceptance to zero
      }
    }
    /* After the burnin period */
    if ((g+1)%DIV==0 && (g+1)>NBURN) {
      // beta
      for (int p=0; p<NP; p++) {
        Ar_beta[p]=((double) nA_beta[p])/DIV;
        nA_beta[p]=0.0; // We reinitialize the number of acceptance to zero
      }
      // gamma
      for (int q=0; q<NQ; q++) {
        Ar_gamma[q]=((double) nA_gamma[q])/DIV;
        nA_gamma[q]=0.0; // We reinitialize the number of acceptance to zero
      }
      // rho
      for (int j=0; j<NCELL; j++) {
        if (viscell[j]>0) {
          Ar_rho[j]=((double) nA_rho[j])/DIV;
          nA_rho[j]=0.0; // We reinitialize the number of acceptance to zero
        }
      }
      // N
      for (int i=0; i<NSITE; i++) {
        Ar_N[i]=((double) nA_N[i])/DIV;
        nA_N[i]=0.0; // We reinitialize the number of acceptance to zero
      }
    }
    
    
    //////////////////////////////////////////////////
    // Progress bar
    double Perc=100*(g+1)/(NGIBBS);
    if (((g+1)%(NGIBBS/100))==0 && verbose[0]==1) {
      Rprintf("*");
      R_FlushConsole();
      //R_ProcessEvents(); for windows
      if (((g+1)%(NGIBBS/10))==0) {
        double mAr_beta=0; // Mean acceptance rate
        double mAr_gamma=0;
        double mAr_rho=0;
        double mAr_N=0;
        // beta
        for (int p=0; p<NP; p++) {
          mAr_beta+=Ar_beta[p]/NP;
        }
        // gamma
        for (int q=0; q<NQ; q++) {
          mAr_gamma+=Ar_gamma[q]/NQ;
        }
        // rho
        for (int j=0; j<NCELL; j++) {
          if (viscell[j]>0) {
            mAr_rho+=Ar_rho[j]/NVISCELL;
          }
        }
        // N
        for (int i=0; i<NSITE; i++) {
          mAr_N+=Ar_N[i]/NSITE;
        }
        Rprintf(":%.1f%%, mean accept. rates= beta:%.3f, gamma:%.3f, rho:%.3f, N:%.3f\n",Perc,mAr_beta,mAr_gamma,mAr_rho,mAr_N);
        R_FlushConsole();
        //R_ProcessEvents(); for windows
      }
    }
    
    
    //////////////////////////////////////////////////
    // User interrupt
    R_CheckUserInterrupt(); // allow user interrupt 	    
    
  } // Gibbs sampler
  
  //////////////////////////
  // Rounding N_pred if save.N==0
  if (save_N[0]==0) {
    for (int i=0; i<NSITE; i++) {
      N_pred[i]= (int)(N_pred_double[i] < 0 ? (N_pred_double[i]-0.5):(N_pred_double[i]+0.5));
    }
  }
  
  ///////////////
  // Delete memory allocation (see malloc())
  /* Obs */
  free(dens_data.Y);
  /* Site */
  free(dens_data.IdSiteforObs);
  free(dens_data.nObsSite);
  for (int i=0; i<NSITE; i++) {
    free(dens_data.ListObsBySite[i]);
  }
  free(dens_data.ListObsBySite);
  /* Latent variable */
  free(dens_data.N_run);
  free(N_pred_double);
  /* Spatial correlation */
  free(dens_data.IdCellforSite);
  free(dens_data.nSiteCell);
  for (int j=0; j<NCELL; j++) {
    free(dens_data.ListSiteByCell[j]);
  }
  free(dens_data.ListSiteByCell);
  free(dens_data.nNeigh);
  for (int j=0; j<NCELL; j++) {
    free(dens_data.Neigh[j]);
  }
  free(dens_data.Neigh);
  free(dens_data.rho_run);
  free(viscell);
  /* Suitability */
  for (int i=0; i<NSITE; i++) {
    free(dens_data.X[i]);
  }
  free(dens_data.X);
  free(dens_data.mubeta);
  free(dens_data.Vbeta);
  free(dens_data.beta_run);
  free(lambda_run);
  /* Observability */
  for (int n=0; n<NOBS; n++) {
    free(dens_data.W[n]);
  }
  free(dens_data.W);
  free(dens_data.mugamma);
  free(dens_data.Vgamma);
  free(dens_data.gamma_run);
  free(delta_run);
  /* Predictions */
  free(IdCell_pred);
  for (int m=0; m<NPRED; m++) {
    free(X_pred[m]);
  }
  free(X_pred);
  free(lambda_pred_run);
  /* Adaptive MH */
  free(sigmap_beta);
  free(nA_beta);
  free(Ar_beta);
  free(sigmap_gamma);
  free(nA_gamma);
  free(Ar_gamma);
  free(sigmap_rho);
  free(nA_rho);
  free(Ar_rho);
  free(nA_N);
  free(Ar_N);
  /* Random seed */
  gsl_rng_free(r);
  
} // end hSDM function


////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
