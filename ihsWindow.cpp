/*
* ihsWindows.cpp
* Extracts largest |iHS| score and calculates proportion of scores in windows greater than a cutoff
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
#include <map>

void info (unsigned int winsize, double ihs_cutoff) {
	std::cout << "\nUsage:\n"
	<< "ihsWindow [selscan normalized iHS *.norm file] [options]\n"
	<< "\nAssumes iHS locus ID in format chr*_position\n"
	<< "\nOptions:\n"
	<< "-winsize INT Window size (bp) [" << winsize << "]\n"
	<< "-cutoff FLOAT Determine fraction of sites with |iHS| > cutoff [" << ihs_cutoff << "]\n"
	<< "-chrlen FILE TSV-file with columns (1) chr (2) chromosome length (bp), and each row is a different chromosome\n"
	<< "\nOutput:\n"
	<< "(1) chromosome\n"
	<< "(2) window start\n"
	<< "(3) window stop\n"
	<< "(4) most extreme iHS score\n"
	<< "(5) extreme iHS position\n"
	<< "(6) proportion |iHS| > cutoff\n"
	<< "(7) Number SNPs in window\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &ihsfile, unsigned int &windowsz, double &ihs_cutoff, std::ifstream &chrfile) {

        int rv = 0;
        if (argc < 2) {
                std::cerr << "Must supply iHS input file\n";
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
		} else if (strcmp(argv[argpos],"-chrlen") == 0) {
			chrfile.open(argv[argpos+1]);
			if (!chrfile) {
				std::cerr << "Unable to open chromosome length file " << argv[argpos+1] << "\n";
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
	std::cout << chr << "\t" << winstart << "\t" << winend << "\t";
	if (nsites > 0) {
		double pbig = nbig/nsites;
		std::cout << maxihs[1] << "\t" << static_cast<unsigned int>(maxihs[2]) << "\t" << pbig << "\t" << nsites << "\n";
	} else {
		std::cout << "NA\tNA\tNA\t0\n";
	}
}

void setLengthMap(std::ifstream &chrfile, std::map<std::string,unsigned int> &chrmap) {
	std::string chrstr;
	std::string chr;
	unsigned int pos;
	while(getline(chrfile, chrstr)) {
		std::stringstream ss(chrstr);
		ss >> chr >> pos;
		chrmap.insert( std::pair<std::string, unsigned int>(chr, pos));
	}
}

int calciHSWindows(std::ifstream &ihsstream, const unsigned int winsize, const double ihs_cutoff, std::ifstream &chrfile) {
	std::string sitestr;
	std::stringstream ss;
	std::vector<double> sitevec;
	sitevec.resize(6); // ihs norm should have 8 fields
	std::vector<double>::iterator it;
	std::string locus_id;
	std::string chr;
	unsigned int pos;
	unsigned int winstart = 1;
	unsigned int winend = winstart + (winsize-1);
	double maxihs [3]; // 0=|iHS|, 1=iHS, 2=position
	unsigned int nsites = 0;
	int nbig = 0;

	// initialize chromosome length map if using this
	unsigned int chrlen = 0;
 	std::map<std::string, unsigned int> lenmap;
	if (chrfile.is_open()) {
		setLengthMap(chrfile, lenmap);
		chrfile.close();
	}

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
		if (chr.empty()) {
			chr = newchr;
			chrlen = (!lenmap.empty() && lenmap.find(newchr) != lenmap.end()) ? lenmap[chr] : 0;
		}

		// check conditions for window printing
		if (!chr.empty() && newchr != chr) {
			printWindow(chr, winstart, winend, nbig, nsites, maxihs);
			while (winend < chrlen) {
				winstart = winend+1;
				winend = winstart + (winsize-1);
 				if (chrlen && winend > chrlen) winend = chrlen;
				printWindow(chr, winstart, winend, 0, 0, maxihs);
			}
			nsites = 0;
			if (!lenmap.empty() && lenmap.find(newchr) != lenmap.end()) {
				chrlen = lenmap[newchr];
			} else chrlen = 0;
			winstart = 1;
                        winend = winstart + (winsize-1);
			if (chrlen && winend > chrlen) winend = chrlen;
		} else if (pos >= winend) {
			printWindow(chr, winstart, winend, nbig, nsites, maxihs);
			winstart = winend+1;
			winend = winstart + (winsize-1);
			if (chrlen && winend > chrlen) winend = chrlen;
			nsites = 0;
			while (pos > winend) {
				printWindow(chr, winstart, winend, 0, 0, maxihs);
 				winstart = winend+1;
				winend = winstart + (winsize-1);
				if (chrlen && winend > chrlen) winend = chrlen;
			}
		}

		//update window
		double absihs = fabs(sitevec[4]);
		if (nsites == 0) {
			nbig = 0;
			updateMax(sitevec, pos, maxihs);
		}

		if (absihs > maxihs[0]) {
			updateMax(sitevec, pos, maxihs);
		}

		if (absihs > ihs_cutoff) {
			++nbig;
		}

		chr = newchr;
		++nsites;
	}

	// print last windows
	printWindow(chr, winstart, winend, nbig, nsites, maxihs);
	while (winend < chrlen) {
		winstart = winend+1;
		winend = winstart + (winsize-1);
 		if (winend > chrlen) winend = chrlen;
		printWindow(chr, winstart, winend, 0, 0, maxihs);
	}

	return 0;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::ifstream ihsfile;
	std::ifstream chrfile;
	unsigned int winsize = 100000;
	double ihs_cutoff = 2;

	if (!(rv=parseArgs(argc, argv, ihsfile, winsize, ihs_cutoff, chrfile))) {
		rv = calciHSWindows(ihsfile, winsize, ihs_cutoff, chrfile);
	}

	if (ihsfile.is_open()) ihsfile.close();

	return rv;
}
