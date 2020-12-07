/*
* xpehhWindows.cpp
* Extracts most extreme XPEHH score and proportion of XPEHH above/below a threshold in fixed windows
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

void info (unsigned int winsize) {
	std::cout << "\nUsage:\n"
	<< "xpehhWindow <selscan normalized XPEHH *.norm file> <cutoff> [options]\n"
	<< "\nInput file must have locus ID in format chr*_position\n"
	<< "cutoff (FLOAT): Calculate proportion of sites with EXPEHH less (if negative) or greater (if positive) than cutoff\n"
	<< "\nOptions:\n"
	<< "-winsize INT Window size (bp) [" << winsize << "]\n"
	<< "-chrlen FILE TSV-file with columns (1) chr (2) chromosome length (bp), and each row is a different chromosome\n"
	<< "\nOutput:\n"
	<< "(1) chromosome\n"
	<< "(2) window start\n"
	<< "(3) window stop\n"
	<< "(4) minimum (negative cutoff) or maximum (postive cutoff) XPEHH score\n"
	<< "(5) extreme XPEHH position\n"
	<< "(6) proportion XPEHH scores > or < cutoff\n"
	<< "(7) Number SNPs in window\n\n";

/*
* use -chrlen option to print windows up to the end of a chromsome even if they
* do not contain any SNPs (useful for maintaining relative distances in Manhatten plots)
*/
}

int parseArgs (int argc, char** argv, std::ifstream &xpehhfile, double &cutoff, unsigned int &windowsz, std::ifstream &chrfile) {

        int rv = 0;
        if (argc < 3) {
                std::cerr << "Must supply XPEHH file and cutoff value\n";
		info(windowsz);
                return 1;
        }

        xpehhfile.open(argv[1], std::ios_base::in | std::ios_base::binary);
        if (! xpehhfile) {
                std::cerr << "Unable to open XPEHH inpt file " << argv[1] << "\n";
                return -1;
        }

	cutoff = atof(argv[2]);
	if (cutoff == 0) {
		std::cerr << "WARNING: cutoff value of zero will calculate proportion of non-negative XPEHH scores\n";
	}

	int argpos = 3;
	while (argpos < argc) {
        	if (strcmp(argv[argpos], "-winsize") == 0) {
                	windowsz = atoi(argv[argpos+1]);
                	if (windowsz <= 0) {
                        	std::cerr << "Window size must be a positive integer\n";
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

void updateOutlier (const std::vector<double> &siteinfo, const unsigned int pos, double* statinfo) {
	statinfo[0] = siteinfo[6];
	statinfo[1] = siteinfo[5];
	statinfo[2] = pos;
}

void printWindow(std::string &chr, const unsigned int &winstart, const unsigned int &winend, double nbig,
	const unsigned int &nsites, double* outlierstat) {
	std::cout << chr << "\t" << winstart << "\t" << winend << "\t";
	if (nsites > 0) {
		double pbig = nbig/nsites;
		std::cout << outlierstat[0] << "\t" << static_cast<unsigned int>(outlierstat[2]) << "\t" << pbig << "\t" << nsites << "\n";
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

int calcXpehhWindows(std::ifstream &xpehhstream, const double cutoff, const unsigned int winsize, std::ifstream &chrfile) {
	std::string sitestr;
	std::stringstream ss;
	std::vector<double> sitevec;
	sitevec.resize(8); // expehh norm should have 10 fields
	std::vector<double>::iterator it;
	std::string locus_id;
	std::string chr;
	unsigned int pos;
	unsigned int winstart = 1;
	unsigned int winend = winstart + (winsize-1);
	double outlierstat [3]; // 0=normalized xpehh, 1=unormalized xpehh, 2=position
	unsigned int nsites = 0;
	int nbig = 0;

	// initialize chromosome length map if using this
	unsigned int chrlen = 0;
 	std::map<std::string, unsigned int> lenmap;
	if (chrfile.is_open()) {
		setLengthMap(chrfile, lenmap);
		chrfile.close();
	}

	// skip header
	getline(xpehhstream, sitestr);
	if (xpehhstream.eof()) {
		std::cerr << "Input XPEHH file had zero sites\n";
		return 0;
	}

	// process the rest of XPEHH file
	while(getline(xpehhstream, sitestr)) {
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
			printWindow(chr, winstart, winend, nbig, nsites, outlierstat);
			while (winend < chrlen) {
				winstart = winend+1;
				winend = winstart + (winsize-1);
				if (chrlen && winend > chrlen) winend = chrlen;
				printWindow(chr, winstart, winend, 0, 0, outlierstat);
			}
			winstart = 1;
			winend = winstart + (winsize-1);
			nsites = 0;
			if (!lenmap.empty() && lenmap.find(newchr) != lenmap.end()) {
				chrlen = lenmap[newchr];
			} else chrlen = 0;
		} else if (pos >= winend) {
			printWindow(chr, winstart, winend, nbig, nsites, outlierstat);
			winstart = winend+1;
			winend = winstart + (winsize-1);
			nsites = 0;
			while (pos > winend) {
				printWindow(chr, winstart, winend, 0, 0, outlierstat);
				winstart = winend+1;
				winend = winstart + (winsize-1);
				if (chrlen && winend > chrlen) winend = chrlen;
			}
		}

		//update window
		double normstat = sitevec[6];
		if (nsites == 0) {
			nbig = 0;
			updateOutlier(sitevec, pos, outlierstat);
		}

		if (cutoff < 0) {
			if (normstat < outlierstat[0]) updateOutlier(sitevec, pos, outlierstat);
			if (normstat < cutoff) ++nbig;
		} else {
			if (normstat > outlierstat[0]) updateOutlier(sitevec, pos, outlierstat);
			if (normstat > cutoff) ++nbig;
		}

		chr = newchr;
		++nsites;
	}

	// print last windows
	printWindow(chr, winstart, winend, nbig, nsites, outlierstat);
	while (winend < chrlen) {
		winstart = winend+1;
		winend = winstart + (winsize-1);
		if (winend > chrlen) winend = chrlen;
		printWindow(chr, winstart, winend, 0, 0, outlierstat);
	}

	return 0;
}

int main (int argc, char** argv) {
	int rv = 0;
	std::ifstream xpehhfile;
	std::ifstream chrfile;
	unsigned int winsize = 100000;
	double cutoff;

	if (!(rv=parseArgs(argc, argv, xpehhfile, cutoff, winsize, chrfile))) {
		rv = calcXpehhWindows(xpehhfile, cutoff, winsize, chrfile);
	}

	if (xpehhfile.is_open()) xpehhfile.close();

	return rv;
}
