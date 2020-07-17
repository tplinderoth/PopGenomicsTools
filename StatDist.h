/*
* StatDist.h
*
* Functions for calculations from statistical distributions
* Many functions implemented following 'Numerical Recipes in C 2nd Ed.'
* and Thorfinn's modifications
*/

# ifndef STATDIST_H_
# define STATDIST_H_

#define ITMAX 1000 // maximum number of iterations
//#define EPS 3.0e-7 // relative accuracy
//#define FPMIN 1.0e-30 // Number near the smallest representable floating-point number

namespace statdist {
void nrerror(const char*); // error handling
double gammln(double xx); // log(gamma(xx))
double gser(const double a, const double x); // incomplete Gamma distribution function P(a,x) evaluated by series representation
double gcf(const double a, const double x); // incomplete Gamma quantile function Q(a,x) evaluated by continued fraction representation
double gammp(double a, double x); // incomplete gamma distribution function P(a,x)
double pchisq (double x, double df); // Distribution function for chi-square
}

#endif /* STATDIST_H_ */
