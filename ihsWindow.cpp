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

void info (unsigned int winsize) {
	std::cout << "\nUsage:\n"
	<< "ihsWindow [selscan normalized iHS *.norm file] [window size (bp)]\n"
	<< "default window size: " << winsize << "\n"
	<< "Assumes iHS locus ID in format chr*_position\n"
	<< "\nOutput:\n"
	<< "(1) chromosome\n"
	<< "(2) window start\n"
	<< "(3) window end\n"
	<< "(4) most extreme iHS score\n"
	<< "(5) extreme iHS position\n"
	<< "(6) average |iHS|\n"
	<< "(7) Number SNPs in window\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &ihsfile, unsigned int &windowsz) {

        int rv = 0;
        if (argc < 2) {
                std::cerr << "Must supply iHS file\n";
		info(windowsz);
                return 1;
        }

        ihsfile.open(argv[1], std::ios_base::in | std::ios_base::binary);
        if (! ihsfile) {
                std::cerr << "Unable to open iHS file " << argv[1] << "\n";
                return -1;
        }

        if (argc > 2) {
                windowsz = atoi(argv[2]);
                if (windowsz <= 0) {
                        std::cerr << "Window size must be a positive integer\n";
                        return -1;
                }
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

void printWindow(std::string &chr, const unsigned int &winstart, const unsigned int &winend, const double &sumihs,
	const unsigned int &nsites, double* maxihs) {
	double avgihs = sumihs/nsites;
	std::cout << chr << "\t" << winstart << "\t" << winend << "\t"
	<< maxihs[1] << "\t" << static_cast<unsigned int>(maxihs[2]) << "\t" << avgihs << "\t" << nsites << "\n";
}

int calciHSWindows(std::ifstream &ihsstream, const unsigned int winsize) {
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
	double sumihs = 0;
	unsigned int nsites = 0;

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
			printWindow(chr, winstart, winend, sumihs, nsites, maxihs);

			nsites = 0;
		}

		//update window
		double absihs = fabs(sitevec[4]);
		if (nsites == 0) {
			winstart = pos;
			sumihs = 0;
			updateMax(sitevec, pos, maxihs);
		}

		if (absihs > maxihs[0]) {
			updateMax(sitevec, pos, maxihs);
		}

		sumihs += absihs;
		winend = pos;
		chr = newchr;
		++nsites;
	}

	// print last window
	printWindow(chr, winstart, winend, sumihs, nsites, maxihs);

	return 0;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::ifstream ihsfile;
	unsigned int winsize = 100000;

	if (!(rv=parseArgs(argc, argv, ihsfile, winsize))) {
		calciHSWindows(ihsfile, winsize);
	}

	if (ihsfile.is_open()) ihsfile.close();

	return rv;
}
