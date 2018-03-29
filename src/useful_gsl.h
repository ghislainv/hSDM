#include <gsl/gsl_rng.h>

// Prototype of useful functions
double logit (double x);
double invlogit (double x);
double my_rtgamma_left(const gsl_rng *r, double a, double b, double t);
