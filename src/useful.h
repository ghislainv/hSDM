#include <gsl/gsl_rng.h>

// Prototype of useful functions
double logit (double x);
double invlogit (double x);
// New with GSL
double myrtgamma_left_gsl(const gsl_rng *r, double a, double b, double t);
