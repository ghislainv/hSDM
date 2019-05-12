#include <gsl/gsl_rng.h>

// Prototype of useful functions
double logit (double x);
double invlogit (double x);
double myrunif();
double myrgamma1(double alpha);
double myrnorm(double mean, double sd);
double myrtgamma_left(double a, double b, double t);
// New with GSL
double myrtgamma_left_gsl(const gsl_rng *r, double a, double b, double t);
