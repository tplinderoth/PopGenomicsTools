/*
 * dxyWindow.cpp - calculates sliding dxy windows from angsd doMaf output
 *
 * Notes:
 * - The ANGSD maf output for both populations should be generated using a fixed reference
 *   allele (-doMajorMinor 4 -ref ), such that the frequency for the same allele is used for
 *   both populations
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
#include <map>

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

void helpInfo (unsigned int winsize, unsigned int stepsize, int minind, int fixedsite) {
	int w1 = 12;
	int w2 = 8;

	std::cout << "\ndxyWindow [options] <pop1 maf file> <pop2 maf file>\n"
	<< "\nOptions:\n"
	<< std::setw(w1) << std::left << "-winsize" << std::setw(w2) << std::left << "INT" << "Window size in base pairs (0 for global calculation) [" << winsize << "]\n"
	<< std::setw(w1) << std::left << "-stepsize" << std::setw(w2) << std::left << "INT" << "Number of base pairs to progress window [" << stepsize << "]\n"
	<< std::setw(w1) << std::left << "-minind" << std::setw(w2) << std::left << "INT" << "Minimum number of individuals in each population with data [" << minind << "]\n"
	<< std::setw(w1) << std::left << "-fixedsite" << std::setw(w2) << std::left << "INT" << "(1) Use fixed number of sites from MAF input for each window (window sizes may vary) or (0) constant window size [" << fixedsite << "]\n"
	<< std::setw(w1) << std::left << "-sizefile" << std::setw(w2) << std::left << "FILE" << "Two-column TSV file with each row having (1) chromsome name (2) chromosome size in base pairs\n" 
	<< "\nNotes:\n"
	<< "* -winsize 1 -stepsize 1 calculates per site dxy\n"
	<< "* -sizefile is REQUIRED(!) with -fixedsite 0 (the default)\n"
	<< "* Both input MAF files need to have the same chromosomes in the same order\n"
	<< "* Assumes SNPs are biallelic across populations\n"
	<< "* For global Dxy calculations only columns 4, 5, and 6 below are printed\n"
	<< "* Input MAF files can contain all sites (including monomorphic sites) or just variable sites\n"
	<< "* -fixedsite 1 -winsize 500 would for example ensure that all windows contain 500 SNPs\n"
	<< "\nOutput:\n"
	<< "(1) chromosome\n"
	<< "(2) Window start\n"
	<< "(3) Window end\n"
	<< "(4) dxy\n"
	<< "(5) number sites in MAF input that were analyzed\n"
	<< "(6) number of sites in MAF input that were skipped due to too few individuals\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &pop1is, std::ifstream &pop2is, unsigned int &winsize, unsigned int &stepsize,
		int &minind, int* infmt, std::ifstream &sizefile, int &fixedsite) {
	int rv = 0;

	if (argc < 3) {
		helpInfo(winsize, stepsize, minind, fixedsite);
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
		} else if (strcmp(argv[argpos], "-sizefile") == 0) {
			const char* sizefname = argv[argpos+1];
			sizefile.open(sizefname);
			if (!sizefile) {
				std::cerr << "Unable to open sizefile: " << sizefname << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos], "-fixedsite") == 0) {
			fixedsite = atoi(argv[argpos+1]);		
		} else {
			std::cerr << "Unknown command: " << argv[argpos] << "\n";
			return -1;
		}

		argpos += 2;
	}

	if (winsize > 0 && stepsize < 1) {
		std::cerr << "Must specify a -stepsize > 0 when -winsize is > 0\n";
		return -1;
	}

	if (!fixedsite && !sizefile) {
		std::cerr << "Must supply size file unless -fixedsite 1\n";
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

int parseSizes (std::ifstream &sizefstream, std::map<std::string, unsigned int>* sizemap) {
	std::string sizeline;
	while (getline(sizefstream, sizeline)) {
		std::stringstream ss(sizeline);
		std::string chrname;
		unsigned int chrlen = 0;
		ss >> chrname >> chrlen;
		if (!chrname.empty() && chrlen > 0) {
			sizemap->insert( std::pair<std::string, unsigned int>(chrname, chrlen));
		} else {
			std::cerr << "Unable to correctly parse chromosome size file\n";
			return -1;		
		}
	}
	return 0;
}

std::vector<sitedata>::iterator calcWindow (std::vector<sitedata> &dxywin, std::string* chr, const unsigned int &winsize, const unsigned int &step, unsigned int* nsites) {

	double dxy = 0;
	unsigned int neffective = 0;
	unsigned int i = 0;
	int nskip = 0;

	for(i=0; i<*nsites; ++i) {
		if (dxywin[i].second >= 0) {
			dxy += dxywin[i].second;
			++neffective;
		} else if (dxywin[i].second == -9) {
			++nskip;		
		}
	}

	// print window info
	std::cout << *chr << "\t" << dxywin[0].first << "\t" << dxywin[*nsites-1].first << "\t" << dxy << "\t" << neffective << "\t" << nskip << "\n";

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

/*
std::vector<sitedata>::iterator calcWindow (std::vector<sitedata> &dxywin, std::string* chr, const unsigned int &winsize, const unsigned int &step, unsigned int* nsites) {
	// old function before adding fixed_sites functionality

	double dxy = 0;
	unsigned int neffective = 0;
	unsigned int i = 0;

	for(i=0; i<*nsites; ++i) {
		if (dxywin[i].second == -9) {
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
*/

int maf2dxy (std::ifstream &pop1is, std::ifstream &pop2is, int* infmt, unsigned int winsize, unsigned int stepsize, int minind, int fixed_site, std::map<std::string, unsigned int> &chrsize) {
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
	unsigned int positer = 1;
	unsigned int lastpos;

	// variables for global dxy
	double dxy_global = 0;
	unsigned int neffective_global = 0;
	unsigned int nskip = 0;

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

		if (winsize > 0 && chr != prevchr) {
		// deal with a new chromosome
			if (!fixed_site) {
			// print windows up to the end of the chromosome
				if (chrsize.find(prevchr) != chrsize.end()) {
					lastpos = chrsize[prevchr];
				} else {
					std::cerr << "Unable to determine size for " << prevchr << "\n";
					return -1;
				}

				while (positer <= lastpos) {
					if (nsites == winsize) winiter = calcWindow(dxywin, &prevchr, winsize, stepsize, &nsites);
					winiter->first = positer;
					winiter->second = -7;
					++winiter;
					++nsites;
					++positer;
				}
				if (nsites > (winsize-stepsize)) {
					winiter = calcWindow(dxywin, &prevchr, winsize, stepsize, &nsites);
				}
		
			} else {
				if (nsites > 0) winiter = calcWindow(dxywin, &prevchr, winsize, stepsize, &nsites);
			}
			positer = 1; // start at first position on new chromosome
		}

		if (winsize > 0 && !fixed_site) {
		// catch window iterator up with the current site with data to process
			while (positer < maf1site.position) {
				if (nsites == winsize) winiter = calcWindow(dxywin, &chr, winsize, stepsize, &nsites);
				winiter->first = positer;
				winiter->second = -7;
				++winiter;
				++nsites;
				++positer;
			}
		}

		// check if window needs to be printed
		if (winsize > 0 && nsites == winsize) {
			winiter = calcWindow(dxywin, &chr, winsize, stepsize, &nsites);		
		}

		// update global dxy
		double dxy = (maf1site.nind >= minind && maf2site.nind >= minind) ? maf1site.freq*(1.0-maf2site.freq) + maf2site.freq*(1.0-maf1site.freq) : -9;
		if (dxy != -9) {
			dxy_global += dxy;
			++neffective_global;
		} else ++nskip;

		// update window
		if (winsize > 0) {
			winiter->first = maf1site.position;
			winiter->second = dxy;
			++winiter;
			++nsites;
			++positer;
		}

		// parse new lines from MAF files
		prevchr = chr;

		if (!getline(p1is, maf1line)) break;
		tokenizeStr(maf1line, &maf1site);

		if (!getline(p2is, maf2line)) break;
		tokenizeStr(maf2line, &maf2site);
	}

	// do last window calculations
	if (!fixed_site) {
		if (chrsize.find(chr) != chrsize.end()) {
			lastpos = chrsize[chr];
		} else {
			std::cerr << "Unable to determine size for " << chr << "\n";
			return -1;
		}

		while (positer <= lastpos) {
			if (nsites == winsize) winiter = calcWindow(dxywin, &chr, winsize, stepsize, &nsites);
			winiter->first = positer;
			winiter->second = -7;
			++winiter;
			++nsites;
			++positer;
		}
	}
	if (nsites > (winsize-stepsize) && nsites <= winsize) {
		winiter = calcWindow(dxywin, &chr, winsize, stepsize, &nsites);
	}

	// print global information
	if (winsize == 0) {
		std::cout << dxy_global << "\t" << neffective_global << "\t" << nskip << "\n";
	} else {
		std::cerr << dxy_global << "\t" << neffective_global << "\t" << nskip << "\n";
	}

	return 0;
}

/*
std::vector<sitedata>::iterator calcWindow (std::vector<sitedata> &dxywin, std::string* chr, const unsigned int &winsize, const unsigned int &step, unsigned int* nsites) {
	// backup of old function before implementing the fixedsites option

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

std::vector<sitedata>::iterator calcWindow_fixedsize (std::vector<sitedata> &dxywin, std::string* chr, const unsigned int &winsize, const unsigned int &step, 
unsigned int &winstart, unsigned int &winend, unsigned int* nsites) {

	double dxy = 0;
	unsigned int neffective = 0;
	unsigned int i = 0;
	int nskip = 0;

	for(i=0; i<*nsites; ++i) {
		if (dxywin[i].second > 0) {
			dxy += dxywin[i].second;
			++neffective;
		} else if (dxywin[9].second == -9) {
			++nskip;		
		}
	}

	// print window info
	std::cout << *chr << "\t" << winstart << "\t" << winend << "\t" << dxy << "\t" << neffective << "\t" << nskip << "\n";
	
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
*/

int main (int argc, char** argv) {
	/*
	 * -winsize 0 calculates global per site dxy
	 * -winsize 1 -stepsize 1 calculates per site dxy
	 */

	int rv = 0;

	std::ifstream pop1is; // input file stream for pop1
	std::ifstream pop2is; // input file stream for pop2
	std::ifstream chrfile; // 2-column file of (1) chromosome name (2) length in base pairs
	unsigned int winsize = 0; // window size
	unsigned int stepsize = 0; // step size
	int minind = 1; // minimum number of individuals with data in each population to calculate dxy for site
	int infmt [2]; // 0 = maf file unzipped, 1 = maf file gzipped
	int fixedsite = 0;
	std::map<std::string, unsigned int> chrsize;

	if ((rv = parseArgs(argc, argv, pop1is, pop2is, winsize, stepsize, minind, infmt, chrfile, fixedsite))) {
		if (rv == 1) rv = 0;
	} else {
		if (!fixedsite) {
			rv = parseSizes(chrfile, &chrsize);
		}
		if (!rv) rv = maf2dxy(pop1is, pop2is, infmt, winsize, stepsize, minind, fixedsite, chrsize);
	}

	if (pop1is.is_open()) pop1is.close();
	if (pop2is.is_open()) pop2is.close();

	return rv;
}
