#include "StatDist.h"
#include <math.h>
#include <iostream>
#include <limits>

const double EPS = std::numeric_limits<double>::epsilon();
const double FPMIN = std::numeric_limits<double>::min()/EPS;

void statdist::nrerror (const char* msg) {
	std::cerr << msg << "\n";
}

double statdist::gammln(double xx) {
	// calculates gamma function as log(gamma(xx)) for xx > 0

	int j;
	double x,tmp,y,ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912,
		14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
		.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
		-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
		.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};

	if (xx <= 0) throw("bad arg in statdist::gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

double statdist::gser(const double a, const double x) {
	/*
	* Returns the incomplete gamma function P(a,x) evaluated by its series
	* representation as gamser. Also returns ln(gamma(a)) as gln.
	*/

	void nerror(char error_text[]);
	int n;
	double sum,del,ap;
	double gln=gammln(a);

	if (x <= 0.0) {
		if (x < 0.0) throw("x less than 0 in routine statdist::gser");
		return 0.0;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1; n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				return sum*exp(-x+a*log(x)-gln);
			}
		}
		throw("Max number of iterations exceeded in routine statdist::gser");
	}

}

double statdist::gcf(const double a, const double x) {
	/*
	* Returns the incomplete gamma function Q(a,x) (complement of P(a,x))
	* evaluated by its continued fraction representation as gammcf.
	* Also return log(gamma(a)) as gln.
	*/

	int i;
	double an,b,c,d,del,h;

	//Set up for evaluating continued fraction by modified Lentz's method with b0=0
	double gln=gammln(a);
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
		if (fabs(del-1.0) <= EPS) break;
	}

	if (i > ITMAX) throw("ERROR: Max number of iterations exceeded in routine statdist::gcf");
	return exp(-x+a*log(x)-gln)*h; // put factors in front
}

double statdist::gammp(double a, double x) {
	// incomplete gamma function P(a,x)

	if (x < 0.0 || a <= 0.0) {
		nrerror("ERROR: Invalid arguments in routine statdist::gammp");
		return -999;
	}

	if (x < (a+1.0)) {
		// use series representation
		try {
			return gser(a,x);
		}
		catch (std::string msg) {
			std::cerr << msg << "\n";
		}
	} else {
		// use the continued fraction representation and return its complement
		try {
			return 1.0 - gcf(a,x);
		}
		catch (std::string msg) {
			std::cerr << msg << "\n";
		}
	}

	return -999;
}

double statdist::pchisq (double x, double df) {
	/*
	* Distribution function for Chi-Square distribution
	* x = statistic value
	* df = degrees of freedom
	*/

	if (x < 0 || df < 0) {
		nrerror("ERROR: Invalid arguments to routine statdist::pchisq");
		return -999;
	}

	return df != 0 ? gammp(0.5*df, 0.5*x) : 0.0;
}
