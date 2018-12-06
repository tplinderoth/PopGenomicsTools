/*
* hetWindow.cpp - calculate heterozygosity for an individual in sliding windows
*
* arguments:
* (1) file with genotypes coded as 0,1,2: (1) chr, (2) position, (3) genotype
* (2) window size (number sites)
* (3) step size (number sites)
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>

typedef std::pair<unsigned int, int> genotype;

void info (unsigned int winsize, unsigned int stepsize) {
	std::cout << "\nUsage:\n"
	<< "hetWindow [genotypes file] [window size (number sites)] [step size (number sites)]\n"
	<< "default window size: " << winsize << "\n"
	<< "default step size: " << stepsize << "\n\n"
	<< "Output:\n"
	<< "(1) chromosome\n"
	<< "(2) window start\n"
	<< "(3) window end\n"
	<< "(4) window midpoint position\n"
	<< "(5) heterozygosity\n"
	<< "(6) Number sites in window\n\n";
}

int parseArgs (int argc, char** argv, std::fstream &genofile, unsigned int &windowsz, unsigned int &stepsz) {

	int rv = 0;
	if (argc < 2) {
		std::cerr << "Must supply genotype file\n";
		return -1;
	}

	genofile.open(argv[1]);
	if (! genofile) {
		std::cerr << "Unable to open genotypes file " << argv[1] << "\n";
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

std::vector<genotype>::iterator calcWindow (std::vector<genotype>* sitegeno, std::string* chr, const unsigned int &winsize, const unsigned int &step, unsigned int *nsites) {
	// find midpoint of window
	unsigned int startpos = (*sitegeno)[0].first;
	unsigned int lastpos = (*sitegeno)[*nsites-1].first;
	unsigned int mid = (startpos + lastpos)/2;

	// calculate heterozygosity
	unsigned int i=0;
	unsigned int nonmissing=0;
	unsigned int nhet=0;

	for(i=0; i<*nsites; ++i) {
		if ((*sitegeno)[i].second >= 0) {
			++nonmissing;
			if ((*sitegeno)[i].second == 1) ++nhet;
		}
	}

	double h = nonmissing != 0 ? (double)nhet/nonmissing : 0.0;

	// print window info
	std::cout << *chr << "\t" << startpos << "\t" << lastpos << "\t"  << mid << "\t" << h << "\t" << *nsites << "\n";

	// prepare window for new values
	static unsigned int nnew;
	if (*nsites == winsize) {
		// same chromosome
		nnew = winsize - step;
		for (i=0; i<nnew; ++i) {
			(*sitegeno)[i].first = (*sitegeno)[step+i].first;
			(*sitegeno)[i].second = (*sitegeno)[step+i].second;
		}
	} else {
		// new chromosome
		nnew = 0;
	}

	*nsites = nnew;
	return (sitegeno->begin() + nnew);
}

int calcHeterozygosity (std::fstream &genofile, unsigned int winsize, unsigned int stepsize) {

	std::vector<genotype> sitegeno;
	sitegeno.resize(winsize);
	std::vector<genotype>::iterator it = sitegeno.begin();

	std::stringstream ss;
	std::string(sitedata);
	std::string(chr);
	std::string(oldchr);

	oldchr = "";
	unsigned int nsites = 0;

	getline(genofile, sitedata);

	while (!sitedata.empty()) {
		// place data for current site into window
		ss.str(std::string());
		ss.clear();
		ss.str(sitedata);
		ss >> chr;

		if (chr != oldchr && nsites > 0) {
			// if new chromosome calculate Fst for current window
			it = calcWindow(&sitegeno, &oldchr, winsize, stepsize, &nsites);
		} else if (nsites == winsize) {
			// if window is full calculate Fst for current window
			it = calcWindow(&sitegeno, &chr, winsize, stepsize, &nsites);
		}

		// update window
		ss >> it->first >> it->second;
		++nsites;
		++it;
		oldchr = chr;

		getline(genofile, sitedata);
	}

	// do last window calculations
	if (nsites > (winsize-stepsize) && nsites <= winsize) {
		it = calcWindow(&sitegeno, &chr, winsize, stepsize, &nsites);
	}

	return 0;
}


int main(int argc, char** argv) {
	int rv = 0;
	std::fstream genofile;
	unsigned int winsize = 1;
	unsigned int stepsize = 1;

	if (argc < 2) {
		info(winsize, stepsize);
		return 0;
	}

	if ((rv=parseArgs(argc, argv, genofile, winsize, stepsize)) != 0)
		return rv;

	rv = calcHeterozygosity(genofile, winsize, stepsize);

	genofile.close();

	return rv;
}
