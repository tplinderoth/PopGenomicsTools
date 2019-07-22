/*
 * dxyWindow.cpp - calculates sliding dxy windows from angsd doMaf output
 *
 * Notes:
 * - The ANGSD maf output for both populations should be generated using a fixed reference
 *   allele (-doMajorMinor 4 -ref ), such that the frequency for the same allele is used for
 *   both populations
 *
 * TODO:
 * 1) Smarter window calculations such that sequence gaps are handled better.
 *    Right now the script just fills up windows of 'winsize' sites, assuming
 *    they are all adjust (which is not true usually), and without regard to
 *    the start or stop positions.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

typedef std::pair<int, double> sitedata;

struct Mafsite {
	std::string chr;
	unsigned int position;
	char major;
	char minor;
	char ref;
	double freq;
	int nind;
};

void helpInfo (unsigned int winsize, unsigned int stepsize, int minind) {
	int w1 = 12;
	int w2 = 8;

	std::cout << "\ndxyWindow [options] [pop1 maf file] [pop2 maf file]\n"
	<< "\nOptions:\n"
	<< std::setw(w1) << std::left << "-winsize" << std::setw(w2) << std::left << "INT" << "Number of sites in window (0 for global calculation) [" << winsize << "]\n"
	<< std::setw(w1) << std::left << "-stepsize" << std::setw(w2) << std::left << "INT" << "Number of sites to progress window [" << stepsize << "]\n"
	<< std::setw(w1) << std::left << "-minind" << std:: setw(w2) << std::left << "INT" << "Minimum number of individuals in each population with data [" << minind << "]\n"
	<< "\nNotes:\n"
	<< "-winsize 1 -stepsize 1 calculates per site dxy\n"
	<< "Both input MAF files need to have the same chromosomes in the same order\n"
	<< "Assumes SNPs are biallelic across populations\n"
	<< "\nOutput:\n"
	<< "(1) chromosome\n"
	<< "(2) Window start\n"
	<< "(3) Window end\n"
	<< "(4) dxy\n"
	<< "(5) number sites in window\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &pop1is, std::ifstream &pop2is, unsigned int &winsize, unsigned int &stepsize,
		int &minind, int* infmt) {
	int rv = 0;

	if (argc < 3) {
		helpInfo(winsize, stepsize, minind);
		return 1;
	}

	// parse input maf files
	unsigned char magic [2] = {0};

	// open pop1 maf file and check whether it's gzipped
	const char* inmaf1 = argv[argc-2];
	pop1is.open(inmaf1, std::ios_base::in | std::ios_base::binary);
	if (! pop1is) {
		std::cerr << "Unable to open Pop1 MAF file: " << inmaf1 << "\n";
		return -1;
	}
	pop1is.read(reinterpret_cast<char*>(magic), sizeof(magic));
	infmt[0] = (magic[0] == 0x1f && magic[1] == 0x8b) ? 1 : 0;
	pop1is.seekg(0, std::ios_base::beg);

	// open pop2 maf file and check whether it's gzipped
	const char* inmaf2 = argv[argc-1];
	pop2is.open(inmaf2, std::ios_base::in | std::ios_base::binary);
	if (! pop2is) {
		std::cerr << "Unable to open Pop2 MAF file: " << inmaf2 << "\n";
		return -1;
	}
	pop2is.read(reinterpret_cast<char*>(magic), sizeof(magic));
	infmt[1] = (magic[0] == 0x1f && magic[1] == 0x8b) ? 1 : 0;
	pop2is.seekg(0, std::ios_base::beg);

	int argpos = 1;
	while (argpos < argc-2) {
		if (strcmp(argv[argpos], "-winsize") == 0) {
			winsize = atoi(argv[argpos+1]);
		} else if (strcmp(argv[argpos], "-stepsize") == 0) {
			stepsize = atoi(argv[argpos+1]);
		} else if (strcmp(argv[argpos], "-minind") == 0) {
			minind = atoi(argv[argpos+1]);
			if (minind <= 0) {
				std::cerr << "-minind must be at least 1\n";
				return -1;
			}
		}

		argpos += 2;
	}

	if (winsize > 0 && stepsize < 1) {
		std::cerr << "Must specify a -stepsize > 0 when -winsize is > 0\n";
		return -1;
	}

	return rv;
}

void tokenizeStr (std::string &str, Mafsite* mafinfo) {
	static std::stringstream ss;
	ss.clear();
	ss.str(str);

	ss >> mafinfo->chr;
	ss >> mafinfo->position;
	ss >> mafinfo->major;
	ss >> mafinfo->minor;
	ss >> mafinfo->ref;
	ss >> mafinfo->freq;
	ss >> mafinfo->nind;
}

std::vector<sitedata>::iterator calcWindow (std::vector<sitedata> &dxywin, std::string* chr, const unsigned int &winsize, const unsigned int &step, unsigned int* nsites) {

	double dxy = 0;
	unsigned int neffective = 0;
	unsigned int i = 0;

	for(i=0; i<*nsites; ++i) {
		if (dxywin[i].second != -9) {
			dxy += dxywin[i].second;
			++neffective;
		}
	}

	// print window info

	if (neffective > 0) {
		std::cout << *chr << "\t" << dxywin[0].first << "\t" << dxywin[*nsites-1].first << "\t" << dxy << "\t" << neffective << "\n";
	} else {
		std::cout << *chr << "\t" << dxywin[0].first << "\t" << dxywin[*nsites-1].first << "\tNA\t" << neffective << "\n";
	}

	// prepare window for new values
	static unsigned int nnew;
	if (*nsites == winsize) {
		// same chromosome
		nnew = winsize - step;
		for (i=0; i<nnew; ++i) {
			dxywin[i].first = dxywin[step+i].first;
			dxywin[i].second = dxywin[step+i].second;
		}
	} else {
		// new chromosome
		nnew = 0;
	}

	*nsites = nnew;
	return (dxywin.begin() + nnew);
}

int maf2dxy (std::ifstream &pop1is, std::ifstream &pop2is, int* infmt, unsigned int winsize, unsigned int stepsize, int minind) {
	// set up zipped file reading if necessary for two maf input files
	std::streambuf *inbuf1 = NULL;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> zipbuf1;
	zipbuf1.push(boost::iostreams::gzip_decompressor());
	if (infmt[0] == 0) {
		inbuf1 = pop1is.rdbuf();
	}
	else if (infmt[0] == 1) {
		zipbuf1.push(pop1is);
		inbuf1 = &zipbuf1;
	}
	std::istream p1is(inbuf1);

	std::streambuf *inbuf2 = NULL;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> zipbuf2;
	zipbuf2.push(boost::iostreams::gzip_decompressor());
	if (infmt[1] == 0) {
		inbuf2 = pop2is.rdbuf();
	}
	else if (infmt[1] == 1) {
		zipbuf2.push(pop2is);
		inbuf2 = &zipbuf2;
	}
	std::istream p2is(inbuf2);

	// skip maf headers and load in first position of both files into vectors

	Mafsite maf1site;
	std::string maf1line;
	getline(p1is, maf1line); // skip header
	getline(p1is, maf1line); // get first position with data
	tokenizeStr(maf1line, &maf1site);

	Mafsite maf2site;
	std::string maf2line;
	getline(p2is, maf2line); // skip header
	getline(p2is, maf2line); // first position with data
	tokenizeStr(maf2line, &maf2site);

	std::string chr = maf1site.chr, prevchr = maf1site.chr;
	if (maf2site.chr != chr) {
		std::cerr << "Chromosomes in MAF files differ\n";
		return -1;
	}

	// process maf files
	unsigned int nsites = 0; // number of sites parsed from maf files for current window
	std::vector<sitedata> dxywin;
	dxywin.resize(winsize);
	std::vector<sitedata>::iterator winiter = dxywin.begin();

	// variables for global dxy
	double dxy_global = 0;
	unsigned int neffective_global = 0;
	unsigned int start_global = 0;
	unsigned int end_global = 0;

	while (!maf1line.empty()) {
		// parse maf files such that their positions are synched, assumes both maf files have the same chromosomes in the same order
		if (maf1site.position != maf2site.position || maf1site.chr != maf2site.chr) {
			if ((maf1site.chr == maf2site.chr && maf1site.position < maf2site.position) || (maf1site.chr != maf2site.chr && maf2site.chr != chr)) {
				// need to catch maf1 up to maf2
				while(maf1site.position != maf2site.position) {
					if (!getline(p1is, maf1line)) break;
					tokenizeStr(maf1line, &maf1site);
				}
				if (maf1site.position != maf2site.position) break; // implies end of maf1 file was reached before catching maf2
			} else {
				// need to catch maf2 up to maf1
				while (maf2site.position < maf1site.position) {
					if (!getline(p2is, maf2line)) break;
					tokenizeStr(maf2line, &maf2site);
				}
				if (maf1site.position != maf2site.position) break; // implies end of maf2 file was reached before catching maf1
			}
		}
		chr = maf1site.chr;

		// check if window needs to be printed
		if (winsize > 0) {
			if (chr != prevchr && nsites > 0) {
				// calculate dxy window from previous chromosome
				winiter = calcWindow(dxywin, &prevchr, winsize, stepsize, &nsites);
			} else if (nsites == winsize) {
				winiter = calcWindow(dxywin, &chr, winsize, stepsize, &nsites);
			}
		}

		// update global dxy
		if (!start_global) start_global = maf1site.position;
		end_global = maf1site.position;
		double dxy = (maf1site.nind >= minind && maf2site.nind >= minind) ? maf1site.freq*(1.0-maf2site.freq) + maf2site.freq*(1.0-maf1site.freq) : -9;
		if (dxy != -9) {
			dxy_global += dxy;
			++neffective_global;
		}

		// update window
		if (winsize > 0) {
			winiter->first = maf1site.position;
			winiter->second = dxy;
			++winiter;
			++nsites;
		}

		// parse new lines from MAF files
		prevchr = chr;

		if (!getline(p1is, maf1line)) break;
		tokenizeStr(maf1line, &maf1site);

		if (!getline(p2is, maf2line)) break;
		tokenizeStr(maf2line, &maf2site);
	}

	// do last window calculations
	if (nsites > (winsize-stepsize) && nsites <= winsize) {
		winiter = calcWindow(dxywin, &chr, winsize, stepsize, &nsites);
	}

	// print global information
	if (winsize == 0) {
		std::cout << chr << "\t" << start_global << "\t" << end_global << "\t" << dxy_global << "\t" << neffective_global << "\n";
	} else {
		std::cerr << chr << "\t" << start_global << "\t" << end_global << "\t" << dxy_global << "\t" << neffective_global << "\n";
	}

	return 0;
}

int main (int argc, char** argv) {
	/*
	 * -winsize 0 calculates global per site dxy
	 * -winsize 1 -stepsize 1 calculates per site dxy
	 */

	int rv = 0;

	std::ifstream pop1is; // input file stream for pop1
	std::ifstream pop2is; // input file stream for pop2
	unsigned int winsize = 0; // window size
	unsigned int stepsize = 0; // step size
	int minind = 1; // minimum number of individuals with data in each population to calculate dxy for site
	int infmt [2]; // 0 = maf file unzipped, 1 = maf file gzipped

	if ((rv = parseArgs(argc, argv, pop1is, pop2is, winsize, stepsize, minind, infmt))) {
		if (rv == 1) rv = 0;
	} else {
		rv = maf2dxy(pop1is, pop2is, infmt, winsize, stepsize, minind);
	}

	if (pop1is.is_open()) pop1is.close();
	if (pop2is.is_open()) pop2is.close();

	return rv;
}
