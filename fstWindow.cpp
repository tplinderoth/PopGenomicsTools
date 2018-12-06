/*
* fstWindow.cpp - calculate sliding windows for Fst from ANGSD variance components
*
* arguments:
* (1) variance components file
* (2) window size (number sites)
* (3) step size (number sites)
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

struct siteComponents {
	unsigned int pos;
	double a;
	double b;
};

void info (unsigned int winsize, unsigned int stepsize) {
	std::cout << "\nUsage:\n"
	<< "fstWindow [ANGSD fst variance component file] [window size (number sites)] [step size (number sites)]\n"
	<< "default window size: " << winsize << "\n"
	<< "default step size: " << stepsize << "\n\n"
	<< "Output:\n"
	<< "(1) chromosome\n"
	<< "(2) window start\n"
	<< "(3) window end\n"
	<< "(4) window midpoint position\n"
	<< "(5) Fst\n"
	<< "(6) Number sites in window\n\n";
}

int parseArgs (int argc, char** argv, std::fstream &fstfile, unsigned int &windowsz, unsigned int &stepsz) {

	int rv = 0;
	if (argc < 2) {
		std::cerr << "Must supply fst file\n";
		return -1;
	}

	fstfile.open(argv[1]);
	if (! fstfile) {
		std::cerr << "Unable to open Fst variance components file " << argv[1] << "\n";
		return -1;
	}

	if (argc > 2) {
		windowsz = atoi(argv[2]);
		if (windowsz <= 0) {
			std::cerr << "Window size must be a positive integer\n";
			return -1;
		}
	}

	if (argc > 3) {
		stepsz = atoi(argv[3]);
		if (stepsz <= 0) {
			std::cerr << "Step size must be a positive integer\n";
		}
	}

	return rv;
}

std::vector<siteComponents>::iterator calcWindow (std::vector<siteComponents>* sitevar, std::string* chr, const unsigned int &winsize, const unsigned int &step, unsigned int *nsites) {
	// find midpoint of window
	unsigned int startpos = (*sitevar)[0].pos;
	unsigned int endpos = (*sitevar)[*nsites-1].pos;
	unsigned int mid = (startpos + endpos)/2;

	// calculate fst
	double asum = 0;
	double bsum = 0;
	unsigned int i=0;

	for(i=0; i<*nsites; ++i) {
		asum += (*sitevar)[i].a;
		bsum += (*sitevar)[i].b;
	}

	double fst = bsum != 0.0 ? asum/bsum : 0.0;

	// print window info
	std::cout << *chr << "\t" << startpos << "\t" << endpos << "\t" << mid << "\t" << fst << "\t" << *nsites << "\n";

	// prepare window for new values
	static unsigned int nnew;
	if (*nsites == winsize) {
		// same chromosome
		nnew = winsize - step;
		for (i=0; i<nnew; ++i) {
			(*sitevar)[i].pos = (*sitevar)[step+i].pos;
			(*sitevar)[i].a = (*sitevar)[step+i].a;
			(*sitevar)[i].b = (*sitevar)[step+i].b;
		}
	} else {
		// new chromosome
		nnew = 0;
	}

	*nsites = nnew;
	return (sitevar->begin() + nnew);
}

int calcFst (std::fstream &varcomp, unsigned int winsize, unsigned int stepsize) {

	std::vector<siteComponents> sitevar;
	sitevar.resize(winsize);
	std::vector<siteComponents>::iterator it = sitevar.begin();

	std::stringstream ss;
	std::string(sitedata);
	std::string(chr);
	std::string(oldchr);

	oldchr = "";
	unsigned int nsites = 0;

	getline(varcomp, sitedata);

	while (!sitedata.empty()) {
		// place data for current site into window
		ss.str(std::string());
		ss.clear();
		ss.str(sitedata);
		ss >> chr;

		if (chr != oldchr && nsites > 0) {
			// if new chromosome calculate Fst for current window
			it = calcWindow(&sitevar, &oldchr, winsize, stepsize, &nsites);
		} else if (nsites == winsize) {
			// if window is full calculate Fst for current window
			it = calcWindow(&sitevar, &chr, winsize, stepsize, &nsites);
		}

		// update window
		ss >> it->pos >> it->a >> it->b;
		++nsites;
		++it;
		oldchr = chr;

		getline(varcomp, sitedata);
	}

	// do last window calculations
	if (nsites > (winsize-stepsize) && nsites <= winsize) {
		it = calcWindow(&sitevar, &chr, winsize, stepsize, &nsites);
	}

	return 0;
}


int main(int argc, char** argv) {
	int rv = 0;
	std::fstream varf;
	unsigned int winsize = 1;
	unsigned int stepsize = 1;

	if (argc < 2) {
		info(winsize, stepsize);
		return 0;
	}

	if ((rv=parseArgs(argc, argv, varf, winsize, stepsize)) != 0)
		return rv;

	rv = calcFst(varf, winsize, stepsize);

	varf.close();

	return rv;
}
