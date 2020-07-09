#include "StatDist.h"
#include <math.h>
#include <iostream>

void statdist::nrerror (const char* msg) {
	std::cerr << msg << "\n";
}

double statdist::gammln(double xx) {
	// calculates gamma function as log(gamma(xx)) for xx > 0

	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (int j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void statdist::gser(double *gamser, double a, double x, double *gln) {
	/*
	* Returns the incomplete gamma function P(a,x) evaluated by its series
	* representation as gamser. Also returns ln(gamma(a)) as gln.
	*/

	void nerror(char error_text[]);
	int n;
	double sum,del,ap;

	if (x <= 0.0) {
		if (x < 0.0) nrerror("ERROR: x less than 0 in routine statdist::gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1; n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("ERROR: Max number of iterations exceeded in routine statdist::gser");
		return;
	}

}

void statdist::gcf(double *gammcf, double a, double x, double *gln) {
	/*
	* Returns the incomplete gamma function Q(a,x) (complement of P(a,x))
	* evaluated by its continued fraction representation as gammcf.
	* Also return log(gamma(a)) as gln.
	*/

	int i;
	double an,b,c,d,del,h;

	//Set up for evaluating continued fraction by modified Lentz's method with b0=0
	*gln=gammln(a);
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;

	// iterate to convergence
	for(i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c = b + an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}

	if (i > ITMAX) nrerror("ERROR: Max number of iterations exceeded in routine statdist::gcf");
	*gammcf = exp(-x+a*log(x)-(*gln))*h; // put factors in front
}

double statdist::gammp(double a, double x) {
	// incomplete gamma function P(a,x)

	double gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0) {
		nrerror("ERROR: Invalid arguments in routine statdist::gammp");
		return -999;
	}

	if (x < (a+1.0)) {
		// use series representation
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		// use the continued fraction representation and return its complement
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

double statdist::pchisq (double x, int df) {
	/*
	* Distribution function for Chi-Square distribution
	* x = statistic value
	* df = degrees of freedom
	*/

	if (x < 0 || df < 0) {
		nrerror("ERROR: Invalid arguments to routine statdist::pchisq");
		return -999;
	}

	return df != 0 ? gammp(df/2.0, x/2.0) : 0.0;
}
