////////////////////////////////////////////////////////////////////
//
// hSDM.poisson.iCAR.c
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
//
// Revisions: 
// - G. Vieilledent, on November 15th 2012
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
// My own functions
#include "useful.h"


/* ********************************************************************* */
/* dens_par */

struct dens_par {
    /* Data */
    int NOBS;
    int NCELL;
    int *Y;
    int *IdCell;
    int *nObsCell;
    int **PosCell;
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
	/* prob_p */
	double Xpart_prob_p=0.0;
	for (int p=0; p<d->NP; p++) {
	    if (p!=k) {
		Xpart_prob_p+=d->X[n][p]*d->beta_run[p];
	    }
	}
	Xpart_prob_p+=d->X[n][k]*beta_k;
	double prob_p=exp(Xpart_prob_p+d->rho_run[d->IdCell[n]]);
	/* log Likelihood */
	logL+=dpois(d->Y[n],prob_p,1);
    }
    // logPosterior=logL+logPrior
    double logP=logL+dnorm(beta_k,d->mubeta[k],sqrt(d->Vbeta[k]),1);
    return logP;
}


/* ************************************************************ */
/* rhodens */

static double rhodens (double rho_i, void *dens_data) {
    // Pointer to the structure: d 
    struct dens_par *d;
    d=dens_data;
    // Indicating the rank of the parameter of interest
    int i=d->pos_rho; //
    // logLikelihood
    double logL=0;
    for (int m=0; m<d->nObsCell[i]; m++) {
	int w=d->PosCell[i][m]; // which observation
	/* prob_p */
	double Xpart_prob_p=0.0;
	for (int p=0; p<d->NP; p++) {
	    Xpart_prob_p+=d->X[w][p]*d->beta_run[p];
	}
	double prob_p=exp(Xpart_prob_p+rho_i);
	/* log Likelihood */
	logL+=dpois(d->Y[w],prob_p,1);
    }
    // logPosterior=logL+logPrior
    int nNeighbors=d->nNeigh[i];
    double sumNeighbors=0.0;
    for (int m=0;m<nNeighbors;m++) {
	sumNeighbors+=d->rho_run[d->Neigh[i][m]];
    }
    double meanNeighbors=sumNeighbors/nNeighbors;
    double logP=logL+dnorm(rho_i,meanNeighbors,sqrt(d->Vrho_run/nNeighbors),1); 
    return logP;
}


/* ************************************************************ */
/* Gibbs sampler function */

void hSDM_poisson_iCAR (
	
    // Constants and data
    const int *ngibbs, int *nthin, int *nburn, // Number of iterations, burning and samples
    const int *nobs, // Number of observations
    const int *ncell, // Number of cells
    const int *np, // Number of fixed effects for prob_p
    const int *Y_vect, // Count
    const double *X_vect, // Suitability covariates
    // Spatial correlation
    const int *C_vect, // Cell Id
    const int *nNeigh, // Number of neighbors for each cell
    const int *Neigh_vect, // Vector of neighbors sorted by cell  
    // Starting values for M-H
    const double *beta_start,
    const double *rho_start,
    // Parameters to save
    double *beta_vect,
    double *rho_pred,
    double *Vrho,
    // Defining priors
    const double *mubeta, double *Vbeta,
    const double *priorVrho,
    const double *shape, double *rate,
    const double *Vrho_max,
    // Diagnostic
    double *Deviance,
    double *prob_p_pred, // Proba of suitability
    // Seeds
    const int *seed,
    // Verbose
    const int *verbose
    
    ) {
	
    ////////////////////////////////////////////////////////////////////////////////
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Defining and initializing objects

    ////////////////////////////////////////
    // Initialize random number generator // 
    srand(seed[0]);                       
    
    ///////////////////////////
    // Redefining constants //
    const int NGIBBS=ngibbs[0];
    const int NTHIN=nthin[0];
    const int NBURN=nburn[0];
    const int NSAMP=(NGIBBS-NBURN)/NTHIN;
    const int NOBS=nobs[0];
    const int NCELL=ncell[0];
    const int NP=np[0];

    ///////////////////////////////////
    // Declaring some useful objects //
    double *prob_p=malloc(NOBS*sizeof(double));
    for (int n=0; n<NOBS; n++) {
	prob_p[n]=0.0;
    }

    //////////////////////////////////////////////////////////
    // Set up and initialize structure for density function //
    struct dens_par dens_data;

    /* Data */
    dens_data.NOBS=NOBS;
    dens_data.NCELL=NCELL;
    // Y
    dens_data.Y=malloc(NOBS*sizeof(int));
    for (int n=0; n<NOBS; n++) {
	dens_data.Y[n]=Y_vect[n];
    }

    /* Spatial correlation */
    // IdCell
    dens_data.IdCell=malloc(NOBS*sizeof(int));
    for (int n=0; n<NOBS; n++) {
	dens_data.IdCell[n]=C_vect[n];
    }
    // nObsCell
    dens_data.nObsCell=malloc(NCELL*sizeof(int));
    for (int i=0; i<NCELL; i++) {
	dens_data.nObsCell[i]=0;
	for (int n=0; n<NOBS; n++) {
	    if (dens_data.IdCell[n]==i) {
		dens_data.nObsCell[i]++;
	    }
	}
    }
    // PosCell
    dens_data.PosCell=malloc(NCELL*sizeof(int*));
    for (int i=0; i<NCELL; i++) {
	dens_data.PosCell[i]=malloc(dens_data.nObsCell[i]*sizeof(int));
	int repCell=0;
	for (int n=0; n<NOBS; n++) {
	    if (dens_data.IdCell[n]==i) {
		dens_data.PosCell[i][repCell]=n;
		repCell++;
	    }
	}
    }
    // Number of neighbors by cell
    dens_data.nNeigh=malloc(NCELL*sizeof(int));
    for (int i=0; i<NCELL; i++) {
	dens_data.nNeigh[i]=nNeigh[i];
    }
    // Neighbor identifiers by cell
    int posNeigh=0;
    dens_data.Neigh=malloc(NCELL*sizeof(int*));
    for (int i=0; i<NCELL; i++) {
	dens_data.Neigh[i]=malloc(nNeigh[i]*sizeof(int));
        for (int m=0; m<nNeigh[i]; m++) {
	    dens_data.Neigh[i][m]=Neigh_vect[posNeigh+m];
	}
	posNeigh+=nNeigh[i];
    }
    dens_data.pos_rho=0;
    dens_data.rho_run=malloc(NCELL*sizeof(double));
    for (int i=0; i<NCELL; i++) {
	dens_data.rho_run[i]=rho_start[i];
    }
    dens_data.shape=shape[0];
    dens_data.rate=rate[0];
    dens_data.Vrho_run=Vrho[0];

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

    ////////////////////////////////////////////////////////////
    // Proposal variance and acceptance for adaptive sampling //

    // beta
    double *sigmap_beta = malloc(NP*sizeof(double));
    int *nA_beta = malloc(NP*sizeof(int));
    for (int p=0; p<NP; p++) {
	nA_beta[p]=0;
	sigmap_beta[p]=1.0;
    }
    double *Ar_beta = malloc(NP*sizeof(double)); // Acceptance rate 
    for (int p=0; p<NP; p++) {
	Ar_beta[p]=0.0;
    }

    // rho
    double *sigmap_rho = malloc(NCELL*sizeof(double));
    int *nA_rho = malloc(NCELL*sizeof(int));
    for (int i=0; i<NCELL; i++) {
	nA_rho[i]=0;
	sigmap_rho[i]=1.0;
    }
    double *Ar_rho = malloc(NCELL*sizeof(double)); // Acceptance rate 
    for (int i=0; i<NCELL; i++) {
	Ar_rho[i]=0.0;
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
	    double x_prop=myrnorm(x_now,sigmap_beta[p]);
	    double p_now=betadens(x_now, &dens_data);
	    double p_prop=betadens(x_prop, &dens_data);
	    double r=exp(p_prop-p_now); // ratio
	    double z=myrunif();	
	    // Actualization
	    if (z < r) {
		dens_data.beta_run[p]=x_prop;
		nA_beta[p]++;
	    }
	}


	////////////////////////////////////////////////
	// rho
	
        /* Sampling rho_run[i] */
	for (int i=0; i<NCELL; i++) {
	    dens_data.pos_rho=i; // Specifying the rank of the parameter of interest
	    double x_now=dens_data.rho_run[i];
	    double x_prop=myrnorm(x_now,sigmap_rho[i]);
	    double p_now=rhodens(x_now, &dens_data);
	    double p_prop=rhodens(x_prop, &dens_data);
	    double r=exp(p_prop-p_now); // ratio
	    double z=myrunif();
	    // Actualization
	    if (z < r) {
		dens_data.rho_run[i]=x_prop;
		nA_rho[i]++;
	    }
	}

	/* Centering rho_run[i] */
	double rho_sum=0.0;
	for (int i=0; i<NCELL; i++) {
	    rho_sum+=dens_data.rho_run[i];
	}
	double rho_bar=rho_sum/NCELL;
	for (int i=0; i<NCELL; i++) {
	    dens_data.rho_run[i]=dens_data.rho_run[i]-rho_bar;
	}


	////////////////////////////////////////////////
	// Vrho

	if (priorVrho[0]>0.0) {      // fixed value for Vrho
	    dens_data.Vrho_run=priorVrho[0];
	}
	else {
	    double Sum=0.0;
	    for (int i=0; i<NCELL; i++) {
		double Sum_neigh=0.0;
		double nNeigh=dens_data.nNeigh[i];
		double rho_run=dens_data.rho_run[i];
		for (int m=0; m<nNeigh; m++) {
		    Sum_neigh += dens_data.rho_run[dens_data.Neigh[i][m]];
		}
		Sum += rho_run*(nNeigh*rho_run-Sum_neigh);
	    }
	    if (priorVrho[0]==-1.0) { // prior = 1/Gamma(shape,rate)
		double Shape=shape[0]+0.5*(NCELL-1);
		double Rate=rate[0]+0.5*Sum;
		dens_data.Vrho_run=Rate/myrgamma1(Shape);
	    }
	    if (priorVrho[0]==-2.0) { // prior = Uniform(0,Vrho_max)
		double Shape=0.5*NCELL-1;
		double Rate=0.5*Sum;
		dens_data.Vrho_run=1/myrtgamma_left(Shape,Rate,1/Vrho_max[0]);
	    }
	}


	//////////////////////////////////////////////////
	// Deviance

	// logLikelihood
	double logL=0.0;
	for (int n=0; n<NOBS; n++) {
	    /* prob_p */
	    double Xpart_prob_p=0.0;
	    for (int p=0; p<NP; p++) {
		Xpart_prob_p+=dens_data.X[n][p]*dens_data.beta_run[p];
	    }
	    prob_p[n]=exp(Xpart_prob_p+dens_data.rho_run[dens_data.IdCell[n]]);
	    /* log Likelihood */
	    logL+=dpois(dens_data.Y[n],prob_p[n],1);
	}

	// Deviance
	double Deviance_run=-2*logL;


	//////////////////////////////////////////////////
	// Output
	if (((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
	    int isamp=((g+1)-NBURN)/(NTHIN);
	    for (int p=0; p<NP; p++) {
		beta_vect[p*NSAMP+(isamp-1)]=dens_data.beta_run[p];
	    }
	    Deviance[isamp-1]=Deviance_run;
	    for (int n=0; n<NOBS; n++) {
		prob_p_pred[n]+=prob_p[n]/NSAMP; // We compute the mean of NSAMP values
	    }
	    for (int i=0; i<NCELL; i++) {
		rho_pred[i]+=dens_data.rho_run[i]/NSAMP; // We compute the mean of NSAMP values
	    }
	    Vrho[isamp-1]=dens_data.Vrho_run;
	}


	///////////////////////////////////////////////////////
	// Adaptive sampling (on the burnin period)
	const double ropt=0.24;
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
	    for (int i=0; i<NCELL; i++) {
		Ar_rho[i]=((double) nA_rho[i])/DIV;
		if(Ar_rho[i]>=ropt) sigmap_rho[i]=sigmap_rho[i]*(2-(1-Ar_rho[i])/(1-ropt));
		else sigmap_rho[i]=sigmap_rho[i]/(2-Ar_rho[i]/ropt);
		nA_rho[i]=0.0; // We reinitialize the number of acceptance to zero
	    }
	}
        /* After the burnin period */
	if((g+1)%DIV==0 && (g+1)>NBURN){
	    for (int p=0; p<NP; p++) {
		Ar_beta[p]=((double) nA_beta[p])/DIV;
		nA_beta[p]=0.0; // We reinitialize the number of acceptance to zero
	    }
	    for (int i=0; i<NCELL; i++) {
		Ar_rho[i]=((double) nA_rho[i])/DIV;
		nA_rho[i]=0.0; // We reinitialize the number of acceptance to zero
	    }
	}

    
	//////////////////////////////////////////////////
	// Progress bar
	double Perc=100*(g+1)/(NGIBBS);
	if (((g+1)%(NGIBBS/100))==0 && (*verbose==1)) {  
	    Rprintf("*");
	    R_FlushConsole();
	    //R_ProcessEvents(); for windows
	    if(((g+1)%(NGIBBS/10))==0){
		double mAr_beta=0; // Mean acceptance rate
		double mAr_rho=0;
		for (int p=0; p<NP; p++) {
		    mAr_beta+=Ar_beta[p]/NP;
		}
		for (int i=0; i<NCELL; i++) {
		    mAr_rho+=Ar_rho[i]/NCELL;
		}
		Rprintf(":%.1f%%, mean accept. rates= beta:%.3f, rho:%.3f\n",Perc,mAr_beta,mAr_rho);
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
    free(dens_data.IdCell);
    free(dens_data.nObsCell);
    for (int i=0; i<NCELL; i++) {
	free(dens_data.PosCell[i]);
    }
    free(dens_data.PosCell);
    /* Spatial correlation */
    free(dens_data.nNeigh);
    for (int i=0; i<NCELL; i++) {
    	free(dens_data.Neigh[i]);
    }
    free(dens_data.Neigh);
    free(dens_data.rho_run);
    /* Suitability */
    for (int n=0; n<NOBS; n++) {
    	free(dens_data.X[n]);
    }
    free(dens_data.X);
    free(dens_data.mubeta);
    free(dens_data.Vbeta);
    free(dens_data.beta_run);
    free(prob_p);
    /* Adaptive MH */
    free(sigmap_beta);
    free(nA_beta);
    free(Ar_beta);
    free(sigmap_rho);
    free(nA_rho);
    free(Ar_rho);

} // end hSDM function


////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
