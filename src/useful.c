#include <stdlib.h> // For rand() and RAND_MAX
#include <math.h>


/*****************************************************************/
/* General functions */
/*****************************************************************/

/**********************
/* Function logit */
double logit (double x) {
    double Result=log(x)-log(1-x);
    return Result;
}

/**********************
/* Function invlogit */
double invlogit (double x) {
    double Result=1/(1+exp(-x));
    return Result;
}


/*****************************************************************/
/* Random draws */
/*****************************************************************/

/**********************
/* myrunif() */
double myrunif() {
    return ((double)rand() + 0.5)/((double)RAND_MAX + 1.0);
}

/**********************
/* myrgamma1(shape,rate=1), function modified from the Scythe C++ library */
double myrgamma1 (double alpha) {

    double accept;
    int test;
    double u, v, w, x, y, z, b, c;
    
    // Implement Best's (1978) simulator
    b = alpha - 1;
    c = 3 * alpha - 0.75;
    test = 0;
    while (test == 0) {
	u = myrunif();
	v = myrunif();
	
	w = u * (1 - u);
	y = sqrt (c / w) * (u - .5);
	x = b + y;

	if (x > 0) {
            z = 64 * pow (v, 2) * pow (w, 3);
            if (z <= (1 - (2 * pow (y, 2) / x))) {
		test = 1;
		accept = x;
            } else if ((2 * (b * log (x / b) - y)) >= log (z)) {
		test = 1;
		accept = x;
            } else {
		test = 0;
            }
	}
    }
    
    return accept;
}

/**********************
/* rnorm1, Knuth's 2nd volume of TAOCP 3rd edition page 122 */
double rnorm1 () {
  double v1,v2,s;

  do {
    v1 = 2.0 * ((double) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((double) rand()/RAND_MAX) - 1;

    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );

  if (s == 0.0)
    return 0.0;
  else
    return (v1*sqrt(-2.0 * log(s) / s));
}

/* myrnorm */
double myrnorm (double mean, double sd) {
    return (mean + rnorm1() * sd);
}

/**********************
/* left truncated gamma
/* Philippe, A. Simulation of right and left truncated gamma distributions by mixtures Statistics and Computing, 1997, 7, 173-181

/* the following function returns a random number from TG^+(a,b,1)            */
/* a=shape parameter, b=rate parameter                                        */
/* when a is an interger                                                      */
/* See Devroye, L. (1985)  Non-Uniform Random Variate Generation              */
/* Springer-Verlag, New-York.                                                 */
double integer(double a, double b)                             
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
        u=myrunif();
        i=1;
        while(u>wlc[i]){i=i+1;};
        x=myrgamma1( (double) i)/b+1.0;        
        free(wl);
        free(wlc);               
        return(x);
}

/* the following function returns a random number from TG^+(a,b,1)            */  
/* a=shape parameter, b=rate parameter                                        */ 
double inter_le(double a, double b)
{
        double test=0,u,x,y,M;
        if (a<1.0)
            {
            M=1.0;
            while (test == 0)
                         {
                         x=1-(1/b)*log(1-myrunif());
                         y=1/pow(x,1-a);
                         if (myrunif()< y/M) test=1.0;
                         };
            }
            else
            {
            if (a<b) 
                      {
                       M=exp(floor(a)-a);
                       while (test == 0)
                            {
                            x=integer(floor(a), b*floor(a)/a);
                            y=pow(x, a-floor(a))*exp(-x*b*(1-floor(a)/a));
                            if (myrunif()< y/M) test=1.0;
                            };
                       }
                       else
                       {
                        M=exp(floor(a)-a)*pow(a/b,a-floor(a));
                        while (test == 0)
                           {
                            x=integer(floor(a), b+floor(a)-a);
                            y=pow(x, a-floor(a))*exp(-x*(-floor(a)+a));
                            if (myrunif()< y/M) test=1.0;
                            };
		       };
	   };
         return(x);
}

/* the following function returns a random number from TG^+(a,b,t)            */
/* a=shape parameter, b=rate parameter                                        */
double myrtgamma_left(double a, double b, double t)                                 
{
        return(inter_le(a, b*t)*t);
}
