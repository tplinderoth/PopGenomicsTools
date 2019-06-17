/*
* ihsWindows.cpp
* Extracts largest |iHS| score and average |iHS| in nonoverlapping windows
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>

void info (unsigned int winsize, double ihs_cutoff) {
	std::cout << "\nUsage:\n"
	<< "ihsWindow [selscan normalized iHS *.norm file] [options]\n"
	<< "\nAssumes iHS locus ID in format chr*_position\n"
	<< "\nOptions\n"
	<< "-winsize INT Window size (bp) [" << winsize << "]\n"
	<< "-cutoff FLOAT Determine fraction of sites with |iHS| > cutoff [" << ihs_cutoff << "]\n"
	<< "\nOutput:\n"
	<< "(1) chromosome\n"
	<< "(2) window start\n"
	<< "(3) window end\n"
	<< "(4) most extreme iHS score\n"
	<< "(5) extreme iHS position\n"
	<< "(6) proportion |iHS| > cutoff\n"
	<< "(7) Number SNPs in window\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &ihsfile, unsigned int &windowsz, double &ihs_cutoff) {

        int rv = 0;
        if (argc < 2) {
                std::cerr << "Must supply iHS file\n";
		info(windowsz, ihs_cutoff);
                return 1;
        }

        ihsfile.open(argv[1], std::ios_base::in | std::ios_base::binary);
        if (! ihsfile) {
                std::cerr << "Unable to open iHS file " << argv[1] << "\n";
                return -1;
        }

	int argpos = 2;
	while (argpos < argc) {
        	if (strcmp(argv[argpos], "-winsize") == 0) {
                	windowsz = atoi(argv[argpos+1]);
                	if (windowsz <= 0) {
                        	std::cerr << "Window size must be a positive integer\n";
                        	return -1;
                	}
        	}
		else if (strcmp(argv[argpos],"-cutoff") == 0) {
			ihs_cutoff = atof(argv[argpos+1]);
			if (ihs_cutoff < 0) {
				std::cerr << "|iHS| cutoff must be >= zero\n";
				return -1;
			}
		} else {
			std::cerr << "Unknown argument " << argv[argpos] << "\n";
			return -1;
		}

		argpos += 2;
	}

        return rv;
}

std::string& extractChr (std::string &locus) {
	static std::string chr;
	chr.clear();
	for (unsigned int i=0; i<locus.size(); ++i) {
		if (locus[i] == '_') {
			break;
		} else {
			chr.push_back(locus[i]);
		}
	}

	return chr;
}

void updateMax (std::vector<double> &siteinfo, unsigned int pos, double* maxinfo) {
	maxinfo[0] = fabs(siteinfo[4]);
	maxinfo[1] = siteinfo[4];
	maxinfo[2] = pos;
}

void printWindow(std::string &chr, const unsigned int &winstart, const unsigned int &winend, double nbig,
	const unsigned int &nsites, double* maxihs) {
	double pbig = nbig/nsites;
	std::cout << chr << "\t" << winstart << "\t" << winend << "\t"
	<< maxihs[1] << "\t" << static_cast<unsigned int>(maxihs[2]) << "\t" << pbig << "\t" << nsites << "\n";
}

int calciHSWindows(std::ifstream &ihsstream, const unsigned int winsize, const double ihs_cutoff) {
	std::string sitestr;
	std::stringstream ss;
	std::vector<double> sitevec;
	sitevec.resize(6); // ihs norm should have 8 fields
	std::vector<double>::iterator it;
	std::string locus_id;
	std::string chr;
	unsigned int pos;
	unsigned int winstart = 0, winend = 0;
	double maxihs [3]; // 0=|iHS|, 1=iHS, 2=position
	unsigned int nsites = 0;
	int nbig = 0;

	while(getline(ihsstream, sitestr)) {
		// extract site information
		ss.clear();
		ss.str(sitestr);
		ss >> locus_id;
		std::string newchr = extractChr(locus_id);
		ss >> pos;
		it=sitevec.begin();
		while (ss >> *it) {
			++it;
		}

		// check conditions for window printing
		if ((!chr.empty() && newchr != chr) || (nsites > 0 && pos-winstart >= winsize)) {
			if (nsites == 1) {
				winend = winstart+winsize-1;
			}
			printWindow(chr, winstart, winend, nbig, nsites, maxihs);

			nsites = 0;
		}

		//update window
		double absihs = fabs(sitevec[4]);
		if (nsites == 0) {
			winstart = pos;
			nbig = 0;
			updateMax(sitevec, pos, maxihs);
		}

		if (absihs > maxihs[0]) {
			updateMax(sitevec, pos, maxihs);
		}

		if (absihs > ihs_cutoff) {
			++nbig;
		}

		winend = pos;
		chr = newchr;
		++nsites;
	}

	// print last window
	printWindow(chr, winstart, winend, nbig, nsites, maxihs);

	return 0;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::ifstream ihsfile;
	unsigned int winsize = 100000;
	double ihs_cutoff = 3;

	if (!(rv=parseArgs(argc, argv, ihsfile, winsize, ihs_cutoff))) {
		calciHSWindows(ihsfile, winsize, ihs_cutoff);
	}

	if (ihsfile.is_open()) ihsfile.close();

	return rv;
}
