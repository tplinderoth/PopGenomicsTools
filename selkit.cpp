/*
*
* selkit.cpp
*
* Requires bcftools to be installed and in PATH
*
*/

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include "selkit.h"

bool fexists(const char* filename) {
	std::ifstream ifile(filename);
	return ifile;
}

void hkaInfo () {
	int w1 = 12;
	std::cerr << "\nselkit hka [arguments]\n"
	<< "\nArguments:\n"
	<< std::setw(w1) << std::left << "-vcf" << "VCF file to analyze (must be bgzipped and indexed)\n"
	<< std::setw(w1) << std::left << "-popfile" << "TSV file with rows having VCF sample name followed by a 0 (outgroup) or 1 (ingroup) population identifier\n"
	<< std::setw(w1) << std::left << "-rf" << "TSV file with regions in format CHR POS END to calculate HKA statistic for\n"
	<< std::setw(w1) << std::left << "-r" << "Region supplied as a string in format 'chr:from-to' to calculate HKA statistic for\n"
	<< std::setw(w1) << std::left << "-out" << "Name of output file\n"
	<< std::setw(w1) << std::left << "-include" << "TSV file with regions in format CHR POS END to use for expectation parameter estimation [default: all sites]\n"
	<< std::setw(w1) << std::left << "-exclude" << "TSV file with regions in format CHR POS END to exclude for expectation calculation\n"
	<< std::setw(w1) << std::left << "-ps" << "Probability of a neutral site being polymorphic within the ingroup\n"
	<< std::setw(w1) << std::left << "-pd" << "Probability that a neutral site is fixed between the ingroup and outgroup\n"
	<< std::setw(w1) << std::left << "-passonly" << "Use only sites with PASS in VCF FILTER field\n"
	<< "\nHKA statistic calculated by comparing all ingroup individuals to all outgroup individuals.\n"
	<< "If supplying precalculated probabilities -ps and -pd must sum to 1.\n"
	<< "\n";
}

std::string formatRegion(std::string s) {
	std::string sfmt("");
	size_t l = s.length();
	size_t c = s.find(':');
	size_t d = s.find_last_of('-');
	if (c == std::string::npos || d == std::string::npos) {
		std::cerr << "ERROR: Invalid region format " << s << "\n";
		return "";
	}

	sfmt += s.substr(0,c) + "\t";
	sfmt += s.substr(c+1,(d-c-1)) + "\t";
	sfmt += s.substr(d+1,(l-d-1));

	return sfmt;
}

int hkaArgs (int argc, char** argv, std::string &vcf, std::string &popfile, std::string &rf, std::string &rstr, std::string &out, std::string &exfile, std::string &incfile,
	int &passonly, double &ps, double &pd, int &preprob) {
	if (argc < 3 || strcmp(argv[2],"-help") == 0) {
		hkaInfo();
		return 1;
	}

	int c = 2;
	while (c < argc) {
		if (strcmp(argv[c],"-vcf") == 0) {
			vcf = argv[c+1];
			if (!fexists(vcf.c_str())) {
				std::cerr << "Unable to locate VCF file " << vcf << "\n";
				return -1;
			}
		} else if (strcmp(argv[c],"-popfile") == 0) {
			popfile = argv[c+1];
			if (!fexists(popfile.c_str())) {
				std::cerr << "Unable to locate popfile " << popfile << "\n";
				return -1;
			}
		} else if (strcmp(argv[c],"-rf") == 0) {
			rf = argv[c+1];
			if (!fexists(rf.c_str())) {
				std::cerr << "Unable to locate file of regions to calculate statistics for: " << rf << "\n";
				return -1;
			}
		} else if (strcmp(argv[c],"-r") == 0) {
			rstr = formatRegion(argv[c+1]);
			if (rstr.empty()) return -1;
		} else if (strcmp(argv[c],"-out") == 0) {
			out = argv[c+1];
		} else if (strcmp(argv[c],"-include") == 0) {
			incfile = argv[c+1];
			if (!fexists(incfile.c_str())) {
				std::cerr << "Unable to locate file of regions to include for estimating parameters: " << incfile << "\n";
				return -1;
			}
		} else if (strcmp(argv[c],"-exclude") == 0) {
			exfile = argv[c+1];
			if (!fexists(exfile.c_str())) {
				std::cerr << "Unable to locate file of regions to exclude: " << exfile << "\n";
			}
		} else if (strcmp(argv[c],"-ps") == 0) {
			ps = atof(argv[c+1]);
			preprob = 1;
		} else if (strcmp(argv[c],"-pd") == 0) {
			pd = atof(argv[c+1]);
			preprob = 1;
		} else if (strcmp(argv[c],"-passonly") == 0) {
			passonly = 1;
			c--;
		} else {
			std::cerr << "Unknown argument: " << argv[c] << "\n";
			return -1;
		}

		c += 2;
	}

	return 0;
}

int bcftools_check () {
	FILE *fp = NULL;
	char version[40];
	fp = popen("bcftools -v","r");
	char* strptr = fgets(version, 9, fp);
	if (fp == NULL || strptr == NULL || strcmp(version,"bcftools") != 0) {
                std::cerr << "selkit hka requires bcftools to be installed and in user's PATH\n";
                return -1;
        }
        pclose(fp);
	return 0;
}

int indexPops (const std::string &vcf, const std::string &popfile, std::vector<int>* popmap) {
	std::ifstream popstream;
	popstream.open(popfile.c_str());
	if (!popstream) {
		std::cerr << "ERROR: Unable to open popfile " << popfile << "\n";
		return -1;
	}

	// parse popfile
	std::map<std::string,int> groupid;
	std::string id;
	int group;
	while (popstream >> id) {
		if(!(popstream >> group)) {
			std::cerr << "ERROR: Popfile is truncated\n";
			return -1;
		}
		groupid.insert(std::pair<std::string,int>(id,group));
	}

	// assign population membership to individual's vcf index
	int buffsize = 5000;
	char buf[buffsize];
	std::stringstream ss;
	std::string cmd("bcftools view -h " + vcf + " | tail -n 1");
	FILE *fp = NULL;
	fp = popen(cmd.c_str(), "r");
	if (fp == NULL) {
		std::cerr << "ERROR: Unable to read VCF " << vcf << "\n";
		return -1;
	}
	while(fgets(buf, buffsize, fp) != NULL) {
		ss << buf;
	}
	pclose(fp);

	std::string tok;
	int i = 0;
	while (ss >> tok) {
		if (i>8) {
			int popval = groupid.find(tok) != groupid.end() ? groupid[tok] : -9;
			popmap->push_back(popval);
		}
		++i;
	}

	return 0;
}

int pop2idx_sub(const std::vector<int> &full, std::vector<int>* subset, int id) {
	if (!subset->empty()) subset->clear();
	subset->reserve(full.size());
	for (size_t i=0; i<full.size(); ++i) {
		if (full[i] == id) subset->push_back(i+9);
	}
	return subset->size();
}

int countAlt(const std::string &geno) {
	// count number of alternate alleles
	int count = 0;

	// check for correct formatting and missing genotypes
	if (geno == ".") {
		return -9;
	} else if (geno.size() >= 3 && (geno[1] == '/' || geno[1] == '|')) {
		if (geno[0] == '.' || geno[1] == '.') return -9;
	} else {
		return -99;
	}

	// count alt alleles
	for (int i = 0; i<3; i+=2) {
		switch (geno[i]) {
			case '0' :
				break;
			case '1' :
				++count;
				break;
			default :
				return -99;
		}
	}

	return count;
}

unsigned long* countVarPatterns (const std::string &cmd, const std::vector<int> &popmap, int passonly, unsigned long *counts) {
	// count number of segregating sites within each population
	// and the number of fixed differences

	for (int i=0; i<3; ++i) counts[i] = 0;

	int buffsize = 1000;
	char buf [buffsize];
	std::string vcfstr;
	int nind = popmap.size()-1;
	FILE *fp = NULL;
	fp = popen(cmd.c_str(),"r");
	double site_counts [2][2];
	int nsites = 0; // counts number of sites processed for region

	while(fgets(buf, buffsize, fp) != NULL) {
		vcfstr += buf;
		if (vcfstr[vcfstr.length()-1] == '\n') {
			// Completed reading in site
			std::stringstream ss(vcfstr);
			std::string tok;
			for (int i = 0; i<9; ++i) ss >> tok; // skip nongenotype info
			int c = 0;
			for (int i=0; i<2; ++i) {
				for (int j=0; j<2; ++j) {
					site_counts[i][j] = 0;
				}
			}
			while(ss >> tok) {
				if (popmap[c] == 0 || popmap[c] == 1) {
					int altcount = countAlt(tok);
					if (altcount >= 0) {
						site_counts[popmap[c]][0] += altcount;
						site_counts[popmap[c]][1] += 2;
					} else if (altcount == -99 ){
						std::cerr << "ERROR: Bad genotype format: " << tok << "\n";
						return NULL;
					}
				}
				++c;
			}
			if (c < nind) {
				std::cerr << c << "\t" << "ERROR: VCF is truncated\n";
				return NULL;
			}

			// count patterns if site has data for both populations
			if (site_counts[0][1] > 0 && site_counts[1][1] > 0) {
				double p0altf = site_counts[0][0]/site_counts[0][1];
				double p1altf = site_counts[1][0]/site_counts[1][1];
				if (p0altf > 0 && p0altf < 1) ++counts[0]; // pop0 variable site
				if (p1altf > 0 && p1altf < 1) ++counts[1]; // pop1 variable site
				if ((p0altf == 0 || p0altf == 1) && (p1altf == 0 || p1altf == 1) && (p0altf != p1altf)) ++counts[2]; // fixed difference
				++nsites;
			}
			vcfstr.clear();
		}
	}

	if (nsites == 0) std::cerr << "WARNING: '" << cmd << "' returned zero sites\n";

	return counts;
}

int expectedParams(const std::string &vcf, const std::string &incfile, const std::string &exfile, const std::vector<int> &popmap, int passonly, double* ps, double* pd) {
	// make bcftools call
	std::string cmd("bcftools view -H");
	if (!incfile.empty()) cmd += (" -R " + incfile);
	if (!exfile.empty()) cmd += (" -T ^" + exfile);
	if (passonly) cmd += " -f PASS";
	cmd += (" " + vcf);

	unsigned long counts[3];
	if (countVarPatterns(cmd, popmap, passonly, counts) == NULL) {
		return -1;
	}

	unsigned long total_counts = 0;
	for (int i = 1; i<3; i++) total_counts += counts[i];
	if (total_counts == 0) {
		std::cerr << "ERROR: Cannot estimate variability parameters, zero informative sites\n";
		return -1;
	}

	*ps = static_cast<double>(counts[1])/total_counts;
	*pd = static_cast<double>(counts[2])/total_counts;

	return 0;
}

double hkaGOF (double s, double d, double es, double ed) {
	/*
	* calculates HKA statistic
	* X2 = (s - es)^2/es + (d - ed)^2/ed
	* s = observed number sites segregating within ingroup
	* es = expected number sites segregating within ingroup
	* d = observed number of fixed difference sites between ingroup and outgroup
	* ed = expected number of fixed difference sites between ingroup and outgroup
	*
	* making the observed counts doubles instead of ints in case they are calculated from likelihoods
	*/

	double x2;

	if (es <= 0 || ed <= 0) {
		std::cerr << "\nERROR:Invalid expected counts for HKA goodness-of-fit calculation\nes: " << es << "\ned: " << ed << "\n\n";
		x2 = -999;
		throw(-1);
	}

	x2 = pow((s - es),2.0)/es + pow((d - ed),2.0)/ed;
	return x2;
}

double* testRegion (const std::string &vcf, const std::string &region, const std::vector<int> &popmap, int passonly, double ps, double pd, double* stats) {

	// stats should be a double array of length 5 that stores [x2, s, es, d, ed]

	unsigned long counts [3];
	unsigned int a, b;
	std::stringstream ss(region);
	std::string chr, start, end;

	while (ss >> chr >> start >> end) {
		// calculate test region length
		std::stringstream start_stream(start);
		start_stream >> a;
		std::stringstream end_stream(end);
		end_stream >> b;
		if (a > b) {
			std::cerr << "ERROR: Region start greater than region end\n";
			return NULL;
		}

		// make bcftools call command
		std::string cmd("bcftools view -H -r " + chr + ":" + start + "-" + end);
		if (passonly) cmd += " -f PASS";
		cmd +=  " " + vcf;

		// count polymorphic and fixed sites
		if (countVarPatterns(cmd, popmap, passonly, counts) == NULL) {
                	return NULL;
		}
		stats[1] = counts[1]; // s
		stats[3] = counts[2]; // d

		unsigned int effective_length = 0;
		for (int i=1; i<3; ++i) effective_length += counts[i];
		std::cerr << "effective length: " << effective_length << "\n";

		// calculate expected counts
		stats[2] = ps * effective_length; // es
		stats[4] = pd * effective_length; // ed

		// calculate hka goodness-of-fit stat
		try {
			stats[0] = hkaGOF(stats[1], stats[3], stats[2], stats[4]);
		}
		catch (int e) {
			return NULL;
		}

		if (ss.rdbuf()->in_avail()) {
			chr.clear();
			start.clear();
			end.clear();
		}
	}

	// make sure region was processed correctly
	if (chr.empty() || start.empty() || end.empty()) {
		std::cerr << "ERROR: Invalid test region\n";
		return NULL;
	}

	return stats;
}

int hka (int argc, char** argv) {
	/*
	* main function for calculating HKA stat
	* Uses genotypes right now but can be easily adjusted to take likelihoods - implement this later
	*/

	int rv = 0;

	// input
	std::string vcf; // name of VCF to analyze
	std::string popfile; // pop sample
	std::string regionfile; // file listing regions to calculate HKA statistic for
	std::string region; // single string 'chr:to-from' region to calculate HKA statistic for
	std::string incfile; // file of regions to use for parameter estimation (versus using all sites)
	std::string exfile; // file of regions to exclude when calculating probabilities for expectations
	std::string outfile; // output file
	int passonly = 0; // if 1 analyze only sites with PASS in FILTER vcf field
	double ps = 0; // probability neutral site is variable within ingroup
	double pd = 0; // probability that neutral site is fixed between ingroup and outgroup
	int preprob = 0; // if 1, probabilities for expectation calculation have been supplied

	// parse arguments
	if ((rv = hkaArgs(argc, argv, vcf, popfile, regionfile, region, outfile, exfile, incfile, passonly, ps, pd, preprob)) < 0) {
		return -1;
	} else if (rv) {
		return 0;
	}

	// check that bcftools is installed and in PATH
	if ((rv=bcftools_check())) {
		return rv;
	}

	// open output stream
	std::ofstream out_stream;
	if (!outfile.empty()) {
		out_stream.open(outfile.c_str());
		if (!out_stream) {
			std::cerr << "ERROR: Unable to open output file " << outfile << "\n";
			return -1;
		}
	} else {
		std::cerr << "ERROR: No output file name, must supply -out\n";
		return -1;
	}

	// open input stream of regions to test
	std::ifstream rf_stream;
	if (!regionfile.empty()) {
		rf_stream.open(regionfile.c_str());
		if (!rf_stream) {
			std::cerr << "ERROR: Unable to open file of regions to test " << regionfile << "\n";
			return -1;
		}
	}


	// initialize groups to compare
	std::vector<int> popmap; // maps sample index in vcf to pop membership
	if((rv=indexPops(vcf, popfile, &popmap))) {
		return rv;
	}

	std::vector<int> p0idx; // field indices in VCF of group 0 individuals
	if (!pop2idx_sub(popmap, &p0idx, 0)) {
		std::cerr << "ERROR: Population 0 individuals not in VCF\n";
		return -1;
	}
	std::vector<int> p1idx; // field indices in VCF of group 1 individuals
	if (!pop2idx_sub(popmap, &p1idx, 1)) {
		std::cerr << "ERROR: Population 1 individuals not in VCF\n";
		return -1;
	}

	// Estimate parameters for obtaining expectations from the genome
	if (!preprob) {
		std::cerr << "Estimating parameters\n";
		if((rv = expectedParams(vcf, incfile, exfile, popmap, passonly, &ps, &pd))) {
			return rv;
		}
		std::cerr << "\nPARAMETER ESTIMATES\n" << "ps: " << ps << "\npd: " << pd << "\n\n";
	}

	if ( fabs(1.0 - ps - pd) > COMPRECISION ) {
		std::cerr << "ERROR: ps and pd do not sum to 1\n";
		return -1;
	}

	// Calculate observed values and calculate statistic
	std::cerr << "Counting observed polymorphisms and calculating HKA statistic\n";
	const int nstats = 5;
	double stats[nstats];
	int df = 1; // degrees of freedom for GOF test
	double p,pval;
	out_stream << "chr\tstart\tend\tX2\tpval\tobsS\texpS\tobsD\texpD\n";

	while (!region.empty() || getline(rf_stream, region)) {
		std::cerr << region << "\n";
		testRegion (vcf, region, popmap, passonly, ps, pd, stats);
		if (stats != NULL) {
			// calculate X2 stat p-value
			if((p = statdist::pchisq(stats[0], df)) == -999) {
				return -1;
			}
			pval = 1.0-p;
			// print results
			out_stream << region << "\t" << stats[0] << "\t" << pval;
			for (int i=1; i<nstats; i++) {
				out_stream << "\t" << stats[i];
			}
			out_stream << "\n";
		} else {
			std::cerr << "Failed to calculate HKA statistic\n";
			return -1;
		}
		region.clear();
	}

	return rv;
}

void mainInfo () {
	int w1 = 8;

	std::cerr << "selkit: A toolbox for popgen selection analyses\n"
	<< "\nselkit [command] [command arguments]\n"
	<< "\nCommands:\n"
	<< std::setw(w1) << std::left << "hka" << "HKA test\n"
	<< "\n";
}

int main (int argc, char** argv) {
	int rv = 0;
	if (argc < 2 || strcmp(argv[1],"-help") == 0) {
		mainInfo();
	} else if (strcmp(argv[1],"hka") == 0) {
		rv = hka(argc, argv);
	} else {
		std::cerr << "ERROR: Unknown command '" << argv[1] << "'\n";
		rv = -1;
	}

	return rv;
}

