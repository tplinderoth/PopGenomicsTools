/*
* StatDist.h
*
* Functions for calculations from statistical distributions
* Many functions implemented following 'Numerical Recipes in C 2nd Ed.'
*/

# ifndef STATDIST_H_
# define STATDIST_H_

#define ITMAX 100 // maximum number of iterations
#define EPS 3.0e-7 // relative accuracy
#define FPMIN 1.0e-30 // Number near the smallest representable floating-point number

namespace statdist {
void nrerror(const char*); // error handling
double gammln(double xx); // log(gamma(xx))
void gser(double *gamser, double a, double x, double *gln); // incomplete Gamma distribution function P(a,x) evaluated by series representation
void gcf(double *gammcf, double a, double x, double *gln); // incomplete Gamma quantile function Q(a,x) evaluated by continued fraction representation
double gammp(double a, double x); // incomplete gamma distribution function P(a,x)
double pchisq (double x, int df); // Distribution function for chi-square
}

#endif /* STATDIST_H_ */
