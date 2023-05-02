/*
* relateStats.cpp
*
* Compile: g++ -O3 -o relateStats relateStats.cpp
*
*/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <math.h>
#include <algorithm>
#include "Matrix.h"

// type defintions
typedef std::unordered_map<std::string, unsigned int> vecmap;
typedef std::vector<std::pair<std::string, double>> sdpvec;

// define classes and class functions

// INDIV CLASS
class indiv
{
public:
	// functions
	indiv(std::string);
	std::string& setid(std::string);
	std::string id () const;
	void setPopValue (int value);
	int popValue () const;
	char sex; // M=male, F=female, '*'=missing
	int cohort; // year class of individual, -999 = missing
	bool group; // 0=exclude individual from population analysis, 1=include individual in pop analysis
	std::vector<std::string> offspring; // vector of offspring IDs
private:
	// members
	std::string _id; // individual ID
	int _popval; // population ID;
};

indiv::indiv (std::string id = "*")
        : sex('N'),
          cohort(-999),
          group(1),
	 _id(id),
	 _popval(0)
{}

std::string& indiv::setid (std::string id = "") {
	if (!id.empty()) _id = id;
	return _id;
}

std::string indiv::id() const {return _id;}

void indiv::setPopValue(int value) {_popval = value;}

int indiv::popValue() const {return _popval;}

// ARRAY CLASS
template <class T>
class Array {
public:
// public functions
T& operator[] (size_t i);
T operator[] (size_t i) const;
void setSize(size_t size);
size_t size() const;
Array ();
Array (const Array& oldarr);
~Array ();

private:
// private member variables
	T* data;
	size_t sz;
	float _lastidx;
	int _flag; // 0: okay, >0: bad
};

template<typename T> T& Array<T>::operator[] (size_t i) {
	if (isnan(_lastidx) || i > (size_t)_lastidx) {
		std::cerr << "Subscript " << i << " out of range \n";
		_flag = 1;
	}
	return data[i];
}

template<typename T> T Array<T>::operator[] (size_t i) const {
	if (isnan(_lastidx) || i > (size_t)_lastidx) {
		std::cerr << "Subscript " << i << " out of range\n";
		_flag = 1;
	}
	return data[i];
}

template<typename T> void Array<T>::setSize(size_t size) {
	sz = size;
	_lastidx = size-1;
	data = new T[size];
	for(unsigned int long i = 0; i < size; ++i) {
		data[i] = 0;
	}
}

template<typename T> size_t Array<T>::size() const {
	return sz;
}

template<typename T> Array<T>::Array ()
	: data(0),
	  sz(0),
	  _lastidx(NAN),
	  _flag(0)
{ }

template<typename T> Array<T>::Array (const Array& oldarr)
           : sz(oldarr.sz)
{
	data = new T[sz];
	for( size_t i = 0; i < sz; ++i)
        	data[i]= oldarr.data[i];
}

template<typename T> Array<T>::~Array () {
	delete [] data;
	data = 0;
}


// define functions

void helpinfo(const double c, const double min_r) {
	int w1 = 16;
	int w2 = 8;
	int w3 = 2;

	std::cout << "\nrelateStats <input>\n"
	<< "\nPossible Inputs:\n"
	<< std::setw(w1) << std::left << "--pedstat" << std::setw(w2) << std::left << "<INT> Pedigree-based statistics\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "1: Expected genetic contribution from Hunter etal 2019 (requires --anc)\n\n"
	<< std::setw(w1) << std::left << "--skewstat" << std::setw(2) << std::left << "<INT> Genetic skew statistics\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "1: Mosaic FSJ skew statistic\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "2: Matrix proportion\n\n"
	<< std::setw(w1) << std::left << "--out" << std::setw(w2) << std::left << "<STRING> Output name prefix\n"
	<< std::setw(w1) << std::left << "--ped" << std::setw(w2) << std::left << "<FILE> ped format file\n"
        << std::setw(w1) << std::left << "--rmat" << std::setw(w2) << std::left << "<FILE> Relatedness matrix\n"
	<< std::setw(w1) << std::left << "--pop" << std::setw(w2) << std::left << "<FILE> List of individual IDs to keep for analyses\n"
	<< std::setw(w1) << std::left << "--anc" << std::setw(w2) << std::left << "<FILE> List of ancestor IDs\n"
	<< std::setw(w1) << std::left << "--cohort" << std::setw(w2) << std::left << "<FILE> List of individual IDs to restrict genetic representation analyses to\n"
	<< std::setw(w1) << std::left << "--mincohort" << std::setw(w2) << std::left << "<INT> Exclude indviduals with cohort value below INT\n"
	<< std::setw(w1) << std::left << "--maxcohort" << std::setw(w2) << std::left << "<INT> Exclude individuals with cohort value above INT\n"
	<< std::setw(w1) << std::left << "--draw" << std::setw(w2) << std::left << "Output direct descendent pedigrees\n"
	<< std::setw(w1) << std::left << "--background_r" << std::setw(w2) << std::left << "<FLOAT> Background relatedness for Mosaic stat [" << c << "]\n"
	<< std::setw(w1) << std::left << "--min_r" << std::setw(w2) << std::left << "<FLOAT> Consider r values < FLOAT 0 for matrix proportion stat [" << min_r << "]\n"

	<< "\nPedigree statistics:\n"
	<< "--pedstat --out --ped --rmat --anc [--pop] [--cohort] [--mincohort] [--maxcohort] [--draw]\n"
	<< "\nSkew statistics:\n"
	<< "--skewstat --out --rmat --anc [--cohort] [--background_r] [--min_r]"

	<< "\n\nNotes:\n"
	<<"* Assumes first row of relatdness matrix contains individual IDs\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &ped_is, std::ifstream &rmat_is, std::string &outprefix, std::vector <int> &pedstat, std::vector <int> &skewstat, std::ifstream &anc_is,
	std::ifstream &pop_is, int &mincohort, int &maxcohort, int &draw, std::ifstream &cohort_is, double &background_r, double &min_r) {
	int rv = 0;
	int argpos = 1;
	if (argc < 2 || strcmp(argv[argpos], "-h") == 0 || strcmp(argv[argpos],"--help") == 0) {
		helpinfo(background_r, min_r);
		return 1;
	}

	while (argpos < argc) {
		if (strcmp(argv[argpos],"--ped") == 0) {
			const char* pedname = argv[argpos+1];
			ped_is.open(pedname, std::ios_base::in);
			if (!ped_is) {
				std::cerr << "Unable to open ped file: " << pedname << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--rmat") == 0) {
			const char* rmat_name = argv[argpos+1];
			rmat_is.open(rmat_name, std::ios_base::in);
			if (!rmat_is) {
				std::cerr << "Unable to open relatedness matrix file: " << rmat_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--out") == 0) {
			outprefix = argv[argpos+1];
		} else if (strcmp(argv[argpos],"--pedstat") == 0) {
			int s = atoi(argv[argpos+1]);
			if (s > 1) { // change if adding more pedigree stat options
				std::cerr << "Invalid --pedStat option " << s << "\n";
				return -1;
			}
			pedstat.push_back(s);
		} else if (strcmp(argv[argpos],"--skewstat") == 0) {
			int s = atoi(argv[argpos+1]);
			if (s == 1 || s == 2) { // change if adding more skew stat options
			} else {
				std::cerr << "Invalid --skewstat option " << s << "\n";
				return -1;
			}
			skewstat.push_back(s);
		} else if (strcmp(argv[argpos],"--anc") == 0) {
			const char* ancf_name = argv[argpos+1];
			anc_is.open(ancf_name, std::ios_base::in);
			if (!anc_is) {
				std::cerr << "Unable to open file of ancestor IDs: " << ancf_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--pop") == 0) {
			const char* popf_name = argv[argpos+1];
			pop_is.open(popf_name, std::ios_base::in);
			if (!pop_is) {
				std::cerr << "Unable to open file of focal population IDs: " << popf_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos], "--cohort") == 0) {
			const char* cohortf_name = argv[argpos+1];
			cohort_is.open(cohortf_name, std::ios_base::in);
			if (!cohort_is) {
				std::cerr << "Unable to open cohort ID file " << cohortf_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--mincohort") == 0) {
			mincohort = atoi(argv[argpos+1]);
		} else if (strcmp(argv[argpos],"--maxcohort") == 0) {
			maxcohort = atoi(argv[argpos+1]);
		} else if (strcmp(argv[argpos], "--draw") == 0) {
			draw = 1;
			--argpos;
		} else if (strcmp(argv[argpos],"--background_r") == 0) {
			background_r = std::stod(argv[argpos+1]);
			if (background_r < 0 || background_r > 1) {
				std::cerr << "background_r out of range (0,1)\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--min_r") == 0) {
			min_r = std::stod(argv[argpos+1]);
			if (min_r < 0 || min_r > 1) {
				std::cerr << "min_r out of range (0,1)\n";
				return -1;
			}
		} else {
			std::cerr << "Unknown argument: " << argv[argpos] << "\n";
			return -1;
		}
		argpos += 2;
	}

	return rv;
}

size_t getIDs (std::ifstream &is, std::vector<std::string>* ids) {
	std::string line;
	while(getline(is, line)) {
		ids->push_back(line);
	}
	ids->shrink_to_fit();
	return (ids->size());
}

int readPed(std::ifstream &ped_is, std::vector<indiv>* ped, std::unordered_map<std::string, unsigned int>* pedidx) {
	std::string line, tok;

	// figure out optional columns. Set fields are [0] = focal individual ID, [1] = parent 1, [2] = parent 2.
	getline(ped_is,line);
	std::stringstream ss(line);
	float fieldidx [2] = {NAN, NAN}; // [0] = sex field, [1] = cohort field

	int i = 0;
	while(getline(ss, tok, '\t')) {
		std::transform(tok.begin(), tok.end(), tok.begin(), ::toupper);
		if (tok == "SEX") {
			fieldidx[0] = i; // sex
		} else if (tok == "COHORT") {
			fieldidx[1] = i; // cohort
		}
		++i;
	}

	for (i=0; i<2; i++) {
		if (!isnan(fieldidx[i]) && fieldidx[i] < 3) {
			std::cerr << "Invalid field order detected in ped input\n";
			return 1;
		}
	}

	// store individuals
	unsigned int indidx = 0;
	while(getline(ped_is, line)) {
		ss.str(std::string());
		ss.clear();
		ss.str(line);
		i = 0;
		std::string focalid;
		while(getline(ss, tok, '\t')) {
			//std::cout << tok << "\n"; // debug
			if (i == 0) {
				focalid = tok;
				if (pedidx->find(tok) == pedidx->end()) {
				// new individual entry
					ped->push_back(indiv(tok));
					pedidx->insert(std::make_pair(tok,indidx));
					++indidx;
				}
			} else if ((i == 1 || i == 2) && tok != "*") {
				if (pedidx->find(tok) == pedidx->end()) {
					// parent not logged so need to add them
					ped->push_back(indiv(tok));
					pedidx->insert(std::make_pair(tok,indidx));
					++indidx;
				}
				// add focal individual to list of offspring from parent
				(*ped)[(*pedidx)[tok]].offspring.push_back(focalid);
			} else if (!isnan(fieldidx[0]) && i == (int)fieldidx[0]) {
				std::transform(tok.begin(), tok.end(), tok.begin(), ::toupper);
				char sexchar = '*';
				if (tok == "MALE" || tok == "M") {
					sexchar = 'M';
				} else if (tok == "FEMALE" || tok == "F") {
					sexchar = 'F';
				}
				(*ped)[(*pedidx)[focalid]].sex = sexchar;
			} else if (!isnan(fieldidx[1]) && i == (int)fieldidx[1]) {
				for (char const &c : tok) {
					if (!std::isdigit(c)) {
						tok = "*";
						break;
					}
				}
				if (tok != "*") {
					(*ped)[(*pedidx)[focalid]].cohort = stoi(tok);
				}
			}
			++i;
		}
	}

	ped->shrink_to_fit();

	return 0;
}

int readrMat (std::ifstream &rmat_is, std::unordered_map<std::string, unsigned int>* matidx, std::vector<std::string>* matids, Matrix<double> &rmat) {
	std::string line, tok;

	// parse IDs from matrix (this should be row 1)
	getline(rmat_is, line);
	std::stringstream ss(line);
	unsigned int i = 0, j;
	while (ss >> tok) {
		matidx->insert(std::make_pair(tok,i));
		matids->push_back(tok);
		++i;
	}
	matids -> shrink_to_fit();

	std::cerr << "Relatedness matrix header length: " << matids->size() << "\n"; // debug
	//std::cerr << "Number matrix columns: " << matidx->size() << "\n"; // debug

	// parse matrix r values
	rmat.allocate(matidx->size(), matidx->size(), NAN);
	unsigned int lastidx = matidx->size()-1;
	i = 0;
	while(getline(rmat_is, line)) {
		if (i > lastidx) {
			std::cerr << "Matrix row dimensions exceed header length\n";
			return 1;
		}
		ss.str(std::string());
                ss.clear();
                ss.str(line);
		j = 0;
		while (ss >> rmat[i][j]) {
			if (j > lastidx) {
				std::cerr << "Matrix row " << j << " entries exceed header length\n";
				return 1;
			}
			++j;
		}
		if (j != rmat.coln()) {
			std::cerr << "Number of matrix row " << j << " entries does not match header length\n";
			return 1;
		}
		++i;
	}
	if (i != rmat.rown()) {
		std::cerr << "Number of matrix rows does not match header length\n";
		return 1;
	}

	return 0;
}

int hunterStat (std::vector<indiv> &ped, std::unordered_map<std::string, unsigned int> &pedidx, Matrix<double> &rmat,
   std::unordered_map<std::string, unsigned int> &matidx, std::vector<std::string> &matids, std::string &outprefix, std::vector<std::string>* pop,
   std::vector<std::string> *anc, int mincohort, int maxcohort, int draw, std::vector<std::string>* cohort) {
	int rv  = 0;

	if (anc->size() < 1) {
		std::cerr << "Hunter stat requires passing ancestral individuals with --anc\n";
		return -1;
	}

	// open output streams
	std::ofstream outstream((outprefix + ".pedstat1").c_str(), std::ios_base::out);

	std::ofstream drawstream;
	if (draw) drawstream.open((outprefix + ".topo").c_str(), std::ios_base::out);

	// set focal population and count focal size
	int effective_pop_n = 0;
	int cohort_n = 0;
	for (std::vector<indiv>::iterator pediter = ped.begin(); pediter != ped.end(); ++pediter) {
		int v = 1;
		if (pop->size() > 0 && std::find(pop->begin(), pop->end(), pediter->id()) == pop->end()) v = 0;
		pediter->setPopValue(v);
		if (cohort->size() > 0 && std::find(cohort->begin(), cohort->end(), pediter->id()) == cohort->end()) pediter->group = 0;
		if (pediter->popValue() > 0) {
			if ((mincohort != -999 && pediter->cohort < mincohort) || (maxcohort != -999 && pediter->cohort > maxcohort)) continue;
			++effective_pop_n;
			if (pediter->group) ++cohort_n;
		}
	}
	std::cerr << "Focal population size: " << effective_pop_n << "\n";
	std::cerr << "Focal cohort size: " << cohort_n << "\n";

	// check ped file parsing (debug)
	/*
	std::cerr << "ped vector size: " << ped.size() << "\nped index vector size: " << pedidx.size() << "\n";
	for (std::vector<indiv>::iterator it = ped.begin(); it != ped.end(); it++) {
		std::cout << it->id() << " " << it->popValue() << ":";
		for (std::vector<std::string>::iterator child_iter = (*it).offspring.begin(); child_iter != (*it).offspring.end(); child_iter++) {
			std::cout << " " << *child_iter;
		}
		std::cout << "\n";
	}
	*/

	// check matrix parsing (debug)
	/*
	for (std::vector<std::string>::iterator it = matids.begin(); it != matids.end(); ++it) {
		std::cout << ((it == matids.begin())? "" : "\t" ) << *it;
	}
	std::cout << "\n";
	for (size_t i = 0; i < rmat.rown(); ++i) {
		for (size_t j = 0; j < rmat.coln(); ++j) {
			std::cout << ((j == 0) ? "" : "\t") << rmat[i][j];
		}
		std::cout << "\n";
	}
	*/

	// Check for pedigree individuals not present in the relatedness matrix
	for (vecmap::iterator it = pedidx.begin(); it != pedidx.end(); ++it) {
		if (ped[it->second].popValue() != 0 && matidx.find(it->first) == matidx.end()) {
			std::cerr << "warning: " << it->first << " missing from relatedness matrix\n";
		}
	}

	// process pedigree
	std::vector<std::pair<std::string, double>> ind_stats;
	ind_stats.reserve(anc->size());
	double pop_total = 0;
	std::unordered_map<std::string, int> seen;
	seen.reserve(rmat.rown());
	for (std::string const &ancid : *anc) {
		std::vector<std::vector<indiv*>> lineage;
		unsigned int gen = 0;
		lineage.resize(ped.size()+1); // this is maximum possible generations + 1
		lineage[gen].push_back(&ped[pedidx[ancid]]);
		size_t noffspring = lineage[0][0]->offspring.size();
		while (noffspring > 0) {
			size_t offgen = gen+1;
			noffspring = 0;
			for (std::vector<indiv*>::iterator ind_iter = lineage[gen].begin(); ind_iter != lineage[gen].end(); ++ind_iter) {
				if ((*ind_iter)->popValue() == 1) { // screen of population inclusion
					if ((*ind_iter)->cohort != -999 && ((mincohort != -999 && (*ind_iter)->cohort < mincohort) || (maxcohort != -999 && (*ind_iter)->cohort > maxcohort))) {
						continue; // skip if not in set year/cohort range
					}
					for (std::vector<std::string>::iterator offspring_iter = (*ind_iter)->offspring.begin(); offspring_iter != (*ind_iter)->offspring.end(); ++offspring_iter) {
						lineage[offgen].push_back(&ped[pedidx[*offspring_iter]]);
						noffspring++;
					}
				}
			}
			++gen;
		}

		// calculate statistics
		seen.clear();
		double ind_total = 0;
		gen = 0;
		std::stringstream ss;
		//std::cout << "## LINEAGE " << ancid << " ##\n"; // debug
		while (lineage[gen].size() > 0) {
			for (std::vector<indiv*>::iterator ind_iter = lineage[gen].begin(); ind_iter != lineage[gen].end(); ++ind_iter) {
				ss.str(std::string());
				ss.clear();
				ss << matidx[ancid];
				std::string elem(ss.str());
				ss.str(std::string());
				ss.clear();
				ss << matidx[(*ind_iter)->id()];
				elem += "_" + ss.str();
				//std::cout << elem << "\n"; // debug
				if (seen.find(elem) == seen.end()) {
					// relationship has not been recorded yet, so do so
					if ((*ind_iter)->group) {
						ind_total += rmat[matidx[ancid]][matidx[(*ind_iter)->id()]];
					}
					seen.insert({elem,1});
				}
			}
			++gen;
		}
		ind_stats.push_back(std::make_pair(ancid, ind_total));
		pop_total += ind_total;

		// draw pedigrees
		if (draw) {
			gen = 0;
			drawstream << "## LINEAGE " << ancid << " ##\n";
			while (lineage[gen+1].size() > 0) {
				std::string genstr("");
				if (lineage[gen+1].size() > 0) drawstream << gen+1 << ": ";
				for (std::vector<indiv*>::iterator ind_iter = lineage[gen].begin(); ind_iter != lineage[gen].end(); ++ind_iter) {
					if ((*ind_iter)->offspring.size() == 0) continue;
						genstr += "(" + (*ind_iter)->id() + ")";
					for (std::vector<std::string>::iterator offspring_iter = (*ind_iter)->offspring.begin(); offspring_iter != (*ind_iter)->offspring.end(); ++offspring_iter) {
						genstr += " " + *offspring_iter;
					}
					genstr += " ## ";
				}
				drawstream << genstr.substr(0,genstr.length()-4) << "\n";
				++gen;
			}
		}
	}

	// output statistics
	if (pop_total == 0) std::cerr << "warning: Total ancestral contribution is zero\n";
	outstream << "ID\tN_GENOME_COPIES\tGENOME_PROPORTION\n";
	for (sdpvec::iterator it = ind_stats.begin(); it != ind_stats.end(); ++it) {
		outstream << it->first << "\t" << it->second << "\t" << it->second/pop_total << "\n";
	}

	// close output stream
	outstream.close();
	if (drawstream.is_open()) drawstream.close();

	return rv;
}

void relateProb (double* p_relate, double* avg_r, std::vector<std::string>* relatives, double c, size_t cohort_n, const std::vector<std::string>* anc, unsigned int* popidx,
   Matrix<double> &rmat, std::unordered_map<std::string, unsigned int> &matidx, const std::vector<std::string> &matids) {
	// probability that a given ancestor will be related (above level c) to a randomly sampled cohort individual
	unsigned int i = 0;
	for (std::vector<std::string>::const_iterator anciter = anc->begin(); anciter != anc->end(); ++anciter) {
		size_t n_related = 0;
		relatives[i].reserve(cohort_n);
		double r_sum = 0;
		for (unsigned int j=0; j<cohort_n; ++j) {
			double r = rmat[matidx[*anciter]][popidx[j]];
			if (r >= c) {
				++n_related;
				r_sum += r;
				relatives[i].push_back(matids[popidx[j]]);
			}
		}
		p_relate[i] = static_cast<double>(n_related)/cohort_n;
		avg_r[i] = n_related > 0 ? r_sum/n_related : 0;
		++i;
	}
}

void rankProb (double* p_rank, const std::vector<std::string>* anc, size_t cohort_n, double c, unsigned int* popidx, Matrix<double> &rmat, std::unordered_map<std::string,
   unsigned int> &matidx) {
	// probability that a given ancestor will be more related to a randomly sampled cohort individual than any other
	// potential ancestors.

	unsigned int ancidx [anc->size()];

	// find row indices in relatedness matrix of ancestors
	unsigned int i = 0, j = 0;
	for (std::vector<std::string>::const_iterator anciter = anc->begin(); anciter != anc->end(); ++anciter) {
		ancidx[i] = matidx[*anciter];
		++i;
	}

	// for each ancestor count the number of other ancestors with lower relatedness to each cohort individual
	unsigned int d = anc->size()-1;
	for (i=0; i<anc->size(); ++i) {
		double psum = 0;
		size_t nrelatives = 0;
		for (j=0; j<cohort_n; ++j) {
			double r = rmat[ancidx[i]][popidx[j]];
			size_t nlower = 0;
			if (r >= c) {
				for (unsigned int k=0; k<anc->size(); ++k) {
					if (k != i && rmat[ancidx[k]][popidx[j]] < r) ++nlower;
				}
				++nrelatives;
			}
			psum += static_cast<double>(nlower)/d;
		}
		// average probability of being more related to a cohort individual than other ancestors
		p_rank[i] = nrelatives > 0 ? psum/nrelatives : 0;
	}
}

int mosaicStat (Matrix<double> &rmat, std::unordered_map<std::string, unsigned int> &matidx, const std::vector<std::string> &matids, const std::string &outprefix,
   const std::vector<std::string> *anc, const std::vector<std::string>* cohort, const double c) {
	/*
	* calculates skew statistic used in Mosaic FSJ study
	* c: background relatedness
	*/

	int rv = 0;

	if (anc->size() < 1) {
		std::cerr << "Mosaic stat requires passing ancestral individuals with --anc\n";
		return -1;
	}

	// open output streams
	std::ofstream stat_outstream((outprefix + ".skewstat1").c_str(), std::ios_base::out);
	std::ofstream relative_outstream((outprefix + ".relatives").c_str(), std::ios_base::out);

	// find derived population column indexes in relatedness matrix
	size_t cohort_n = 0, i = 0;
	unsigned int popidx [rmat.coln()];
	for (vecmap::const_iterator it = matidx.begin(); it != matidx.end(); ++it) {
		if (std::find(anc->begin(), anc->end(), it->first) == anc->end()) { // exclude ancestral individuals from derived group
			if (cohort->size() == 0 || (cohort->size() > 0 && std::find(cohort->begin(), cohort->end(), it->first) != cohort->end())) {
				// restrict derived pop to list of cohort IDs if provided
				popidx[i] = it->second;
				++i;
			}
		}
	}
	cohort_n = i;
	if (cohort->size() < 1) std::cerr << "Treating all non-ancestral individuals as focal cohort for Mosaic stat calculation\n";
	std::cerr << "focal cohort size: " << cohort_n << "\n" << "Number ancestors: " << anc->size() << "\n";
	std::cerr << "Background relatedness set to " << c << "\n";

	// check that focal IDs were parsed correctly (debug)
	/*
	for (i=0; i<cohort_n; ++i) {
		std::cout << matids[popidx[i]] << "\n";
	}
	*/

	// calculate statistic

	// probability that a given ancestor will be related (above level c) to a randomly sampled cohort individual;
	double p_relate [anc->size()];
	double avg_r [anc->size()];
	std::vector<std::string> relatives[anc->size()];
	relateProb(p_relate, avg_r, relatives, c, cohort_n, anc, popidx, rmat, matidx, matids);

	// probability that a given ancestor will be more related to a randomly sampled cohort individual than any other
	// potential ancestors.
	double p_rank [anc->size()];
	rankProb (p_rank, anc, cohort_n, c, popidx, rmat, matidx);

	// Write results
	stat_outstream << "ID\tSc\tSrank\trelate_prob\tavg_r\trank_prob\n";
	relative_outstream << "ID\trelatives\n";
	for (i=0; i<anc->size(); ++i) {
		double sc = p_relate[i]*avg_r[i];
		double srank = p_relate[i]*p_rank[i];
		stat_outstream << (*anc)[i] << "\t" << sc << "\t" << srank << "\t" << p_relate[i] << "\t" << avg_r[i] << "\t" << p_rank[i] << "\n";
		relative_outstream << (*anc)[i] << "\t";
		if (relatives[i].empty()) {
			relative_outstream << "*\n";
		} else {
			for (std::vector<std::string>::const_iterator it = relatives[i].begin(); it != relatives[i].end(); ++it) {
				relative_outstream << *it;
				if (it == relatives[i].end()-1)
					relative_outstream << "\n";
				else relative_outstream << ",";
			}
		}
	}

	if (stat_outstream.is_open()) stat_outstream.close();
	if (relative_outstream.is_open()) relative_outstream.close();

	return rv;
}

int findMatIndex (unsigned int* idx, vecmap &matidx, const std::vector<std::string>* ids) {
	unsigned int i = 0;
	for(std::vector<std::string>::const_iterator it = ids->begin(); it != ids->end(); ++it) {
		if (matidx.find(*it) == matidx.end()) {
			std::cerr << "error: " << *it << " not in r matrix\n";
			return -1;
		}
		idx[i] = matidx[*it];
		++i;
	}
	return i;
}

void matrixCounts (std::vector<double>* count_p, std::vector<double>* rp, const unsigned int* ancidx, const size_t anc_n, const unsigned int* cohortidx, const size_t cohort_n,
   Matrix<double> &rmat, const double min_r, size_t* count_total = NULL, double* r_total = NULL) {

	// for each ancestor count the number of other ancestors with lower relatedness to each cohort individual
	count_p->clear();
	rp->clear();
	unsigned int i = 0, j = 0;
	for (i=0; i<anc_n; ++i) {
		count_p->push_back(0);
		rp->push_back(0);
		for (j=0; j<cohort_n; ++j) {
			double r = rmat[ancidx[i]][cohortidx[j]];
			(*rp)[i] += r;
			if (r < min_r) r = 0;
			for (unsigned int k=0; k<anc_n; ++k) {
				if (k == i) continue;
				if (rmat[ancidx[k]][cohortidx[j]] < r) ++(*count_p)[i];
			}
		}
		//std::cout << (*count_p)[i] << "\t" << total_counts << "\t" << (*rp)[i] << "\t" << total_r << "\n";
	}

	// calculate normalizations
	unsigned int total_counts = (anc_n * cohort_n) - cohort_n; // this is all pairwise comparisons minus those between the focal ancestor and cohort
	double total_r = 0;
	for (i=0; i<anc_n; ++i) {
		for(j=0; j<cohort_n; ++j) {
			total_r += rmat[ancidx[i]][cohortidx[j]];
		}
	}
	if (count_total) *count_total = total_counts;
	if (r_total) *r_total = total_r;

	// calculate proportions
	i = 0;
	for (std::vector<double>::iterator it = count_p->begin(); it != count_p->end(); ++it) {
		*it /= total_counts;
		(*rp)[i] /= total_r;
		//std::cout << *it << "\t" << (*rp)[i] << "\n";
		++i;
	}
}

int matPstat (Matrix<double> &rmat, std::unordered_map<std::string, unsigned int> &matidx, const std::vector<std::string> &matids, const std::string &outprefix,
   const std::vector<std::string> *anc, const std::vector<std::string>* cohort, const double min_r) {

	int rv = 0;

	// open output streams
	std::ofstream outs((outprefix + ".skewstat2").c_str(), std::ios_base::out);

	// find group indices in relatedness matrix

	// * if neither ancestral or cohort IDs are provided calculate stat for all individuals in matrix using all pairwise r
	// * if ancestral IDs are not provided, but cohort IDs are, ancestral IDs set to all non-cohort individuals in matrix
	// * if only ancestral IDs are provided cohort individuals will be set to all non-ancestral individuals in matrix

	if (anc->empty() && cohort->empty()) {
		cohort = &matids;
		anc = &matids;
	}

	std::vector<std::string> compliment;
	if (!anc->empty() && cohort->empty()) {
		compliment.reserve(rmat.coln()-anc->size());
		for (std::vector<std::string>::const_iterator it = matids.begin(); it != matids.end(); ++it) {
			if (std::find(anc->begin(), anc->end(), *it) == anc->end()) {
				compliment.push_back(*it);
			}
		}
		cohort = &compliment;
	} else if (anc->empty() && !cohort->empty()) {
 		compliment.reserve(rmat.coln()-cohort->size());
		for (std::vector<std::string>::const_iterator it = matids.begin(); it != matids.end(); ++it) {
			if (std::find(cohort->begin(), cohort->end(), *it) == cohort->end()) {
				compliment.push_back(*it);
			}
		}
                anc = &compliment;
	}


	unsigned int ancidx [rmat.coln()];
	int anc_n = findMatIndex(ancidx, matidx, anc);

	unsigned int cohortidx [rmat.coln()];
	int cohort_n = findMatIndex(cohortidx, matidx, cohort);

	if (anc_n < 0 || cohort_n < 0) return -1;

	std::cerr << "ancestral n: " << anc_n << "\ncohort n: " << cohort_n << "\n";

/*
	std::cout << "ANCESTRAL IDs\n";
	for (unsigned int i=0; i<anc_n; ++i) {
		std::cout << matids[ancidx[i]] << "\n";
	}
	std::cout << "COHORT IDs\n";
	for (unsigned int i=0; i<cohort_n; ++i) {
		std::cout << matids[cohortidx[i]] << "\n";
	}
*/

	// calculate statistics
	std::vector<double> count_p;
	count_p.reserve(anc->size());
	std::vector<double> rp;
	rp.reserve(anc->size());

	size_t total_counts = 0;
	double r_total = 0;

	matrixCounts (&count_p, &rp, ancidx, anc_n, cohortidx, cohort_n, rmat, min_r, &total_counts, &r_total);
	std::cerr << "Total pairwise counts: " << total_counts << "\n" << "Total r: " << r_total << "\n";

	// write output
	outs << "ID\tScount\tSrsum\n";
	unsigned int i = 0;
	for (std::vector<double>::const_iterator it = count_p.begin(); it != count_p.end(); ++it) {
		outs << (*anc)[i] << "\t" << *it << "\t" << rp[i] << "\n";
		++i;
	}

	if (outs.is_open()) outs.close();

	return rv;

}

int main (int argc, char** argv) {
	int rv = 0;

	std::ifstream ped_is; // ped format input stream
	std::ifstream rmat_is; // relatedness matrix input stream
	std::string outprefix; // output prefix
	std::vector <int> pedstat; // type of pedigree statistic to calculate
	std::vector <int> skewstat; // type of skew statistic to calculate
	std::ifstream anc_is; // focal individual input stream
	std::ifstream pop_is; // focal population input stream
	std::ifstream cohort_is; // list of IDs to calculate representation among
	int maxcohort = -999, mincohort = -999;
	int draw = 0; // draw direct descendent pedigree if 1
	double background_r = 0.0; // background relatedness level used for Mosaic statistic
	double min_r = 0.0; // minimum r value, currently not used

	// parse arguments
	if ((rv = parseArgs(argc, argv, ped_is, rmat_is, outprefix, pedstat, skewstat, anc_is, pop_is, mincohort, maxcohort, draw, cohort_is, background_r, min_r))) {
		if (rv == 1) rv = 0;
		return rv;
	}

	// parse ID lists
	std::vector<std::string> anc; // vector of ancestral individual IDs
	if (anc_is.is_open() && !getIDs(anc_is, &anc)) {
		std::cerr << "error: Read zero ancestral IDs\n";
		return -1;
	}
	if (anc_is.is_open()) anc_is.close();

	std::vector<std::string> pop; // vector of individual IDs in focal population
	if (pop_is.is_open() && !getIDs(pop_is, &pop)) {
		std::cerr << "error: Read zero focal population IDs\n";
		return -1;
	}
	if (pop_is.is_open()) pop_is.close();

	std::vector<std::string> cohort; // vector of cohort IDs
	if (cohort_is.is_open() && !getIDs(cohort_is, &cohort)) {
		std::cerr << "error: Read zero cohort individual IDs\n";
		return -1;
	}

	// read in relatedness matrix
        std::unordered_map<std::string, unsigned int> matidx;
        std::vector<std::string> matids;
        Matrix<double> rmat;
        if (rmat_is.is_open()) {
		if ((rv = readrMat(rmat_is, &matidx, &matids, rmat))) return rv;
        	std::cerr << "Relatedness matrix dimensions: " << rmat.rown() << " x " << rmat.coln() << "\n";
		rmat_is.close();
	}

	// read in pedigree
        std::vector<indiv> ped;
        std::unordered_map<std::string, unsigned int> pedidx;
	if (ped_is.is_open()) {
		if ((rv = readPed(ped_is, &ped, &pedidx))) return rv;
		ped_is.close();
	}

	// calculate statistics
	std::vector<int>::iterator iter;

	for (iter = pedstat.begin(); iter != pedstat.end(); iter++) {
		if (*iter == 1) {
			std::cerr << "Calculating Hunter expected contribution\n";
			if (anc.size() < 1) {
				std::cerr << "Hunter stat requires passing ancestral individuals with --anc\n";
				return -1;
			}
			if ((rv = hunterStat(ped, pedidx, rmat, matidx, matids, outprefix, &pop, &anc, mincohort, maxcohort, draw, &cohort))) return rv;
		}
	}

	for (iter = skewstat.begin(); iter != skewstat.end(); ++iter) {
		if (*iter == 1) {
			std::cerr << "Calculating Mosaic skew statistic\n";
			if ((rv = mosaicStat(rmat, matidx, matids, outprefix, &anc, &cohort, background_r))) return rv;
		}
		if (*iter == 2) {
			std::cerr << "Calculating matrix proportion skew statistic\n";
			if ((rv = matPstat (rmat, matidx, matids, outprefix, &anc, &cohort, min_r))) return rv;
		}
	}

	return rv;
}
