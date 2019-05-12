// ==============================================================================
// author          :Ghislain Vieilledent
// email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
// web             :https://ecology.ghislainv.fr
// license         :GPLv3
// ==============================================================================

// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// GSL libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/*****************************************************************/
/* General functions */
/*****************************************************************/

/******************/
/* Function logit */
double logit (double x) {
  return log(x) - log(1-x);
}

/*********************/
/* Function invlogit */
double invlogit (double x) {
  if (x > 0) {
    return 1 / (1 + exp(-x));
  }
  else {
    return exp(x) / (1 + exp(x));
  }
}


/*****************************************************************/
/* Random draws with GSL */
/*****************************************************************/


/************************/
/* left truncated gamma */
/* Philippe, A. Simulation of right and left truncated gamma distributions by mixtures Statistics and Computing, 1997, 7, 173-181 */

/* the following function returns a random number from TG^+(a,b,1)            */
/* a=shape parameter, b=rate parameter                                        */
/* when a is an integer                                                       */
/* See Devroye, L. (1985)  Non-Uniform Random Variate Generation              */
/* Springer-Verlag, New-York.                                                 */
double integer_gsl(const gsl_rng *r, double a, double b)                             
{
  double u,x;
  double *wl,*wlc;
  int i;
  wl=(double *)malloc((int)(a+1)*sizeof(double));
  wlc=(double *)malloc((int)(a+1)*sizeof(double));
  wl[1]=1.0;
  wlc[1]=1.0;
  for(i=2; i<= (int)a ;i++ )
  {
    wl[i]=wl[i-1]*(a-i+1)/b;
    wlc[i]=wlc[i-1]+wl[i];
  };
  for(i=1; i<= (int) a; i++)
  {
    wlc[i]=wlc[i]/wlc[(int)a];
  };
  u=gsl_rng_uniform(r);
  i=1;
  while(u>wlc[i]){i=i+1;};
  x=gsl_ran_gamma(r, (double)i, 1.0)/b+1.0;        
  free(wl);
  free(wlc);               
  return(x);
}

/* the following function returns a random number from TG^+(a,b,1)            */  
/* a=shape parameter, b=rate parameter                                        */ 
double inter_le_gsl(const gsl_rng *r, double a, double b)
{
  double test=0,x,y,M;
  if (a<1.0)
  {
    M=1.0;
    while (test == 0)
    {
      x=1-(1/b)*log(1-gsl_rng_uniform(r));
      y=1/pow(x,1-a);
      if (gsl_rng_uniform(r)< y/M) test=1.0;
    };
  }
  else
  {
    if (a<b) 
    {
      M=exp(floor(a)-a);
      while (test == 0)
      {
        x=integer_gsl(r, floor(a), b*floor(a)/a);
        y=pow(x, a-floor(a))*exp(-x*b*(1-floor(a)/a));
        if (gsl_rng_uniform(r)< y/M) test=1.0;
      };
    }
    else
    {
      M=exp(floor(a)-a)*pow(a/b,a-floor(a));
      while (test == 0)
      {
        x=integer_gsl(r, floor(a), b+floor(a)-a);
        y=pow(x, a-floor(a))*exp(-x*(-floor(a)+a));
        if (gsl_rng_uniform(r)< y/M) test=1.0;
      };
    };
  };
  return(x);
}

/* the following function returns a random number from TG^+(a,b,t)            */
/* a=shape parameter, b=rate parameter                                        */
double myrtgamma_left_gsl(const gsl_rng *r, double a, double b, double t)                                 
{
  return(inter_le_gsl(r, a, b*t)*t);
}

