////////////////////////////////////////////////////////////////////
//
// hSDM.Nmixture.c
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
    int *IdSite;
    int *nObsSite;
    int **PosSite;
    /* Latent variable */
    int *N_run;
    int pos_N;
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
	double lambda=exp(Xpart_lambda);
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
	logL+=dbinom(d->Y[n],d->N_run[d->IdSite[n]],delta,1);
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
	int w=d->PosSite[i][m]; // which observation
	/* delta */
	double logit_delta=0.0;
	for (int q=0; q<d->NQ; q++) {
	    logit_delta+=d->W[w][q]*d->gamma_run[q];
	}
	double delta=invlogit(logit_delta);
	/* log Likelihood */
	logL+=dbinom(d->Y[w],N_i,delta,1);
    }
    // logPosterior=logL+logPrior
    /* lambda */
    double Xpart_lambda=0.0;
    for (int p=0; p<d->NP; p++) {
	Xpart_lambda+=d->X[i][p]*d->beta_run[p];
    }
    double lambda=exp(Xpart_lambda);
    double logP=logL+dpois(N_i,lambda,1); 
    return logP;
}

/* ************************************************************ */
/* Gibbs sampler function */

void hSDM_Nmixture (
	
    // Constants and data
    const int *ngibbs, int *nthin, int *nburn, // Number of iterations, burning and samples
    const int *nobs, int *nsite, // Number of observations and sites
    const int *np, // Number of fixed effects for lambda
    const int *nq, // Number of fixed effects for delta
    const int *Y_vect, // Number of successes (presences)
    const double *W_vect, // Observability covariates (nobs x nq)
    const double *X_vect, // Suitability covariates (nsite x np)
    // Spatial sites
    const int *S_vect, // Site Id
    // Predictions
    const int *npred, // Number of predictions
    const double *X_pred_vect, // Suitability covariates for predictions
    // Starting values for M-H
    const double *beta_start,
    const double *gamma_start,
    const int *N_start,
    // Parameters
    double *beta_vect,
    double *gamma_vect,
    int *N_pred,
    // Defining priors
    const double *mubeta, double *Vbeta,
    const double *mugamma, double *Vgamma,
    // Diagnostic
    double *Deviance,
    double *lambda_latent, // Latent proba of suitability (length NSITE)
    double *delta_latent, // Latent proba of observability (length NOBS)
    double *lambda_pred, // Proba of suitability for predictions (length NPRED)
    // Seeds
    const int *seed,
    // Verbose
    const int *verbose,
    // Save p and N
    const int *save_p,
    const int *save_N

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
    const int NSITE=nsite[0];
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
    
    /* Latent variable */
    // N_run
    dens_data.N_run=malloc(NSITE*sizeof(int));
    for (int i=0; i<NSITE; i++) {
    	dens_data.N_run[i]=N_start[i];
    }
    dens_data.pos_N=0;

    /* Sites */
    // IdSite
    dens_data.IdSite=malloc(NOBS*sizeof(int));
    for (int n=0; n<NOBS; n++) {
    	dens_data.IdSite[n]=S_vect[n];
    }
    // nObsSite
    dens_data.nObsSite=malloc(NSITE*sizeof(int));
    for (int i=0; i<NSITE; i++) {
    	dens_data.nObsSite[i]=0;
    	for (int n=0; n<NOBS; n++) {
    	    if (dens_data.IdSite[n]==i) {
    		dens_data.nObsSite[i]++;
    	    }
    	}
    }
    // PosSite
    dens_data.PosSite=malloc(NSITE*sizeof(int*));
    for (int i=0; i<NSITE; i++) {
    	dens_data.PosSite[i]=malloc(dens_data.nObsSite[i]*sizeof(int));
    	int repSite=0;
    	for (int n=0; n<NOBS; n++) {
    	    if (dens_data.IdSite[n]==i) {
    		dens_data.PosSite[i][repSite]=n;
    		repSite++;
    	    }
    	}
    }

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

    // gamma
    double *sigmap_gamma = malloc(NQ*sizeof(double));
    int *nA_gamma = malloc(NQ*sizeof(int));
    double *Ar_gamma = malloc(NQ*sizeof(double)); // Acceptance rate 
    for (int q=0; q<NQ; q++) {
	nA_gamma[q]=0;
	sigmap_gamma[q]=1.0;
	Ar_gamma[q]=0.0;
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
    	// gamma
	
    	for (int q=0; q<NQ; q++) {
    	    dens_data.pos_gamma=q; // Specifying the rank of the parameter of interest
    	    double x_now=dens_data.gamma_run[q];
    	    double x_prop=myrnorm(x_now,sigmap_gamma[q]);
    	    double p_now=gammadens(x_now, &dens_data);
    	    double p_prop=gammadens(x_prop, &dens_data);
    	    double r=exp(p_prop-p_now); // ratio
    	    double z=myrunif();
    	    // Actualization
    	    if (z < r) {
    		dens_data.gamma_run[q]=x_prop;
    		nA_gamma[q]++;
    	    }
    	}


    	////////////////////////////////////////////////
    	// N
	
    	for (int i=0; i<NSITE; i++) {
    	    dens_data.pos_N=i; // Specifying the rank of the parameter of interest
    	    int x_now=dens_data.N_run[i];
	    if (x_now==0) {
    		double s=myrunif();
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
    		    double r=exp(p_prop-p_now);
    		    // Actualization
    		    double z=myrunif();
    		    if (z < r) {
    			dens_data.N_run[i]=x_prop;
    			nA_N[i]++;
    		    }
    		}
    	    }
    	    else {
    		// Proposal
    		double s=myrunif();
    		int x_prop=0;
    		if (s < 0.5) x_prop=x_now-1;
    		else x_prop=x_now+1;
    		// Ratio
    		double p_now=Ndens(x_now, &dens_data);
    		double p_prop=Ndens(x_prop, &dens_data);
    		double r=exp(p_prop-p_now);
    		// Actualization
    		double z=myrunif();
    		if (z < r) {
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
    	    /* delta */
    	    double logit_delta=0.0;
    	    for (int q=0; q<NQ; q++) {
    		logit_delta+=dens_data.W[n][q]*dens_data.gamma_run[q];
    	    }
    	    delta_run[n]=invlogit(logit_delta);
    	    /* log Likelihood */
    	    logL1+=dbinom(dens_data.Y[n],dens_data.N_run[dens_data.IdSite[n]],delta_run[n],1);
    	}
    	double logL2=0.0;
    	for (int i=0; i<NSITE; i++) {
    	    /* lambda */
    	    double Xpart_lambda=0.0;
    	    for (int p=0; p<NP; p++) {
    		Xpart_lambda+=dens_data.X[i][p]*dens_data.beta_run[p];
    	    }
    	    lambda_run[i]=exp(Xpart_lambda);
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
    	    lambda_pred_run[m]=exp(Xpart_lambda_pred);
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
    		double mAr_N=0;
    	    	// beta
    	    	for (int p=0; p<NP; p++) {
    	    	    mAr_beta+=Ar_beta[p]/NP;
    	    	}
    	    	// gamma
    	    	for (int q=0; q<NQ; q++) {
    	    	    mAr_gamma+=Ar_gamma[q]/NQ;
    	    	}
    		// N
    	    	for (int i=0; i<NSITE; i++) {
    		    mAr_N+=Ar_N[i]/NSITE;
    	    	}
    	    	Rprintf(":%.1f%%, mean accept. rates= beta:%.3f, gamma:%.3f, N:%.3f\n",Perc,mAr_beta,mAr_gamma,mAr_N);
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
    /* Data */
    free(dens_data.Y);
    free(dens_data.N_run);
    free(N_pred_double);
    free(dens_data.IdSite);
    free(dens_data.nObsSite);
    for (int i=0; i<NSITE; i++) {
    	free(dens_data.PosSite[i]);
    }
    free(dens_data.PosSite);
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
    free(nA_N);
    free(Ar_N);

} // end hSDM function


////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
