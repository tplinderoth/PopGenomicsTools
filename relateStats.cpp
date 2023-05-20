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
	void setPop (const std::string &p);
	const std::string & pop ();

	// public members
	char sex; // M=male, F=female, '*'=missing
	float cohort; // year class of individual
	float cohort_last; // year (class) that individual is last alive (beyond this they are deceased)
	bool group; // 0=exclude individual from population analysis, 1=include individual in pop analysis
	std::vector<std::string> offspring; // vector of offspring IDs
	std::string parents [2]; // array of [parent1, parent2]
private:
	// private members
	std::string _id; // individual ID
	std::string _pop; // population ID;
};

indiv::indiv (std::string id = "*")
        : sex('N'),
          cohort(NAN),
	  cohort_last(NAN),
          group(1),
          parents({"*"}),
	  _id(id),
	  _pop("")

{}


std::string& indiv::setid (std::string id = "") {
	if (!id.empty()) _id = id;
	return _id;
}

std::string indiv::id() const {return _id;}

void indiv::setPop(const std::string &p) {_pop = p;}

const std::string & indiv::pop() {return _pop;}

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


struct timeSort {
	timeSort(bool descend = false) {this->rev = descend;}
	bool operator () (const std::pair<int, std::vector<std::string>> &p1, const std::pair<int, std::vector<std::string>> &p2) {
		if (rev) return (p1.first > p2.first);
		return (p1.first < p2.first);
	}

	bool rev;
};

// define functions

void helpinfo(const double c, const double min_r, const bool hunter_norm) {
	int w1 = 16;
	int w2 = 8;
	int w3 = 2;

	std::cout << "\nrelateStats <input>\n"
	<< "\nPossible Inputs:\n"
	<< std::setw(w1) << std::left << "--pedstat" << std::setw(w2) << std::left << "<INT> Pedigree-based statistics\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "1: Expected genetic contribution from Hunter etal 2019\n\n"
	<< std::setw(w1) << std::left << "--skewstat" << std::setw(2) << std::left << "<INT> Genetic skew statistics\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "1: Ranked relatedness among relatives\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "2: Relatedness matrix proportion\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "3: Average rank-weighted relatedness\n\n"
	<< std::setw(w1) << std::left << "--out" << std::setw(w2) << std::left << "<STRING> Output name prefix\n"
	<< std::setw(w1) << std::left << "--ped" << std::setw(w2) << std::left << "<FILE> ped format file\n"
        << std::setw(w1) << std::left << "--rmat" << std::setw(w2) << std::left << "<FILE> Relatedness matrix\n"
	<< std::setw(w1) << std::left << "--pop" << std::setw(w2) << std::left << "<FILE> Two-column, tab-delimited file specifying (1) individual ID & (2) population ID\n"
	<< std::setw(w1) << std::left << "--anc" << std::setw(w2) << std::left << "<FILE> List of ancestor IDs\n"
	<< std::setw(w1) << std::left << "--cohort" << std::setw(w2) << std::left << "<FILE> List of individual IDs to restrict genetic representation analyses to\n"
	<< std::setw(w1) << std::left << "--time2" << std::setw(w2) << std::left << "<INT> Sets descendant population to all extant individuals between ancestor's cohort and INT\n"
	<< std::setw(w1) << std::left << "--t2_only" << std::setw(w2) << std::left << "Restrict representation to individuals born during --time2\n"
	<< std::setw(w1) << std::left << "--max_norm" << std::setw(2) << std::left<< "<0|1> Normalize contributions with respect to maximum ancestral cohort contribution if 1 [" << hunter_norm << "]\n"
	<< std::setw(w1) << std::left << "--draw" << std::setw(w2) << std::left << "Output direct descendent pedigrees\n"
	<< std::setw(w1) << std::left << "--background_r" << std::setw(w2) << std::left << "<FLOAT> Background relatedness for skewstat 1 and 3 [" << c << "]\n"
	<< std::setw(w1) << std::left << "--min_r" << std::setw(w2) << std::left << "<FLOAT> Consider r values < FLOAT 0 for skewstat 2 [" << min_r << "]\n"

	<< "\nPedigree statistics:\n"
	<< "--pedstat --out --ped --rmat --anc [--time2] [--t2_only] [--cohort] [--max_norm] [--draw]\n"
	<< "\nSkew statistics:\n"
	<< "--skewstat --out --rmat --anc [--cohort] [--background_r] [--min_r]"

	<< "\n\nNotes:\n"
	<<"* Assumes first row of relatedness matrix contains individual IDs\n\n";
}

void countRemaining(const int &argc, const int &argpos, int n) {
	// n is the length of expected arguments
	static const std::string user_err("error: Unexpected number of user arguments\n");
	if ((argc-1)-(argpos+n) < 0) throw user_err;
}

int parseArgs (int argc, char** argv, std::ifstream &ped_is, std::ifstream &rmat_is, std::string &outprefix, std::vector <int> &pedstat, std::vector <int> &skewstat, std::ifstream &anc_is,
	std::ifstream &pop_is, bool &t2_only, int &draw, std::ifstream &cohort_is, double &background_r, double &min_r, float &time2, bool &max_norm) {
	int rv = 0;
	int argpos = 1;
	if (argc < 2 || strcmp(argv[argpos], "-h") == 0 || strcmp(argv[argpos],"--help") == 0) {
		helpinfo(background_r, min_r, max_norm);
		return 1;
	}

	while (argpos < argc) {
		if (strcmp(argv[argpos],"--ped") == 0) {
			countRemaining(argc, argpos, 1);
			const char* pedname = argv[argpos+1];
			ped_is.open(pedname, std::ios_base::in);
			if (!ped_is) {
				std::cerr << "Unable to open ped file: " << pedname << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--rmat") == 0) {
			countRemaining(argc, argpos, 1);
			const char* rmat_name = argv[argpos+1];
			rmat_is.open(rmat_name, std::ios_base::in);
			if (!rmat_is) {
				std::cerr << "Unable to open relatedness matrix file: " << rmat_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--out") == 0) {
			countRemaining(argc, argpos, 1);
			outprefix = argv[argpos+1];
		} else if (strcmp(argv[argpos],"--pedstat") == 0) {
			int s = atoi(argv[argpos+1]);
			if (s > 1) { // change if adding more pedigree stat options
				std::cerr << "Invalid --pedStat option " << s << "\n";
				return -1;
			}
			pedstat.push_back(s);
		} else if (strcmp(argv[argpos],"--skewstat") == 0) {
			countRemaining(argc, argpos, 1);
			int s = atoi(argv[argpos+1]);
			if (s == 1 || s == 2 || s == 3) { // change if adding more skew stat options
			} else {
				std::cerr << "Invalid --skewstat option " << s << "\n";
				return -1;
			}
			skewstat.push_back(s);
		} else if (strcmp(argv[argpos],"--anc") == 0) {
			countRemaining(argc, argpos, 1);
			const char* ancf_name = argv[argpos+1];
			anc_is.open(ancf_name, std::ios_base::in);
			if (!anc_is) {
				std::cerr << "Unable to open file of ancestor IDs: " << ancf_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--pop") == 0) {
			countRemaining(argc, argpos, 1);
			const char* popf_name = argv[argpos+1];
			pop_is.open(popf_name, std::ios_base::in);
			if (!pop_is) {
				std::cerr << "Unable to open population file: " << popf_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos], "--cohort") == 0) {
			countRemaining(argc, argpos, 1);
			const char* cohortf_name = argv[argpos+1];
			cohort_is.open(cohortf_name, std::ios_base::in);
			if (!cohort_is) {
				std::cerr << "Unable to open cohort ID file " << cohortf_name << "\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--t2_only") == 0) {
			countRemaining(argc, argpos, 0);
			t2_only = 1;
			--argpos;
		} else if (strcmp(argv[argpos], "--draw") == 0) {
			countRemaining(argc, argpos, 0);
			draw = 1;
			--argpos;
		} else if (strcmp(argv[argpos],"--background_r") == 0) {
			countRemaining(argc, argpos, 1);
			background_r = std::stod(argv[argpos+1]);
			if (background_r < 0 || background_r > 1) {
				std::cerr << "background_r out of range (0,1)\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--min_r") == 0) {
			countRemaining(argc, argpos, 1);
			min_r = std::stod(argv[argpos+1]);
			if (min_r < 0 || min_r > 1) {
				std::cerr << "min_r out of range (0,1)\n";
				return -1;
			}
		} else if (strcmp(argv[argpos],"--time2") == 0) {
			countRemaining(argc, argpos, 1);
			time2 = atoi(argv[argpos+1]);
		} else if (strcmp(argv[argpos],"--max_norm") == 0) {
			countRemaining(argc, argpos, 1);
			int v = atoi(argv[argpos+1]);
			if (v == 0 || v == 1) {
				max_norm = v;
			} else {
				std::cerr << "--max_norm only accepts values 0 or 1\n";
				return  -1;
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

template <typename T> void idMap(std::ifstream &is, std::unordered_map<std::string, T> &m) {
	std::string line;
	std::string key;
	T val;

	while(getline(is,line)) {
		std::stringstream ss(line);
		ss >> key;
		ss >> val;
		m.insert({key,val});
	}
}

bool isInt(const std::string &str) {
	int i = 0;
	for (const char &c : str) {
		if (!std::isdigit(c)) {
			if (i == 0 && c == '-') continue; // this allows for negative integers
			return false;
		}
	++i;
	}
	return true;
}

int readPed(std::ifstream &ped_is, std::vector<indiv>* ped, std::unordered_map<std::string, unsigned int>* pedidx, std::vector<std::string>* fields) {
	std::string line, tok;

	// figure out optional columns. Set fields are [0] = focal individual ID, [1] = parent 1, [2] = parent 2.
	getline(ped_is,line);
	std::stringstream ss(line);
	int nfields = 3;
	float fieldidx [nfields] = {NAN, NAN, NAN}; // [0] = sex field, [1] = cohort field, [2] = last cohort that individual is alive field
	int i = 0;
	while(getline(ss, tok, '\t')) {
		std::transform(tok.begin(), tok.end(), tok.begin(), ::toupper);
		if (tok == "SEX") {
			fieldidx[0] = i; // sex
			fields->push_back("sex");
		} else if (tok == "COHORT") {
			fieldidx[1] = i; // cohort
			fields->push_back("cohort");
		} else if (tok == "COHORT_LAST") {
			fieldidx[2] = i; // last alive cohort
			fields->push_back("cohort_last");
		}
		++i;
	}

	for (i=0; i<nfields; i++) {
		if (!isnan(fieldidx[i]) && fieldidx[i] < 3) {
			// first 3 fields should always be ID, parent1, parent2
			std::cerr << "Invalid field order detected in ped input\n";
			return 1;
		}
		// assumes first three fields are ID, PARENT 1, and PARENT 2
		if (i==0) fields->push_back("id");
		if (i==1) fields->push_back("p1");
		if (i==2) fields->push_back("p2");
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
				(*ped)[(*pedidx)[focalid]].parents[i-1] = tok;
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
				if (!isInt(tok)) tok = "*";
				if (tok != "*") (*ped)[(*pedidx)[focalid]].cohort = stoi(tok);
			} else if (!isnan(fieldidx[2]) && i == (int)fieldidx[2]) {
				if (!isInt(tok)) tok = "*";
				if (tok != "*") (*ped)[(*pedidx)[focalid]].cohort_last = stoi(tok);
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

int findMatIndex (unsigned int* idx, vecmap &matidx, const std::vector<std::string>* ids) {
	unsigned int i = 0;
	for(std::vector<std::string>::const_iterator it = ids->begin(); it != ids->end(); ++it) {
		if (matidx.find(*it) == matidx.end()) {
			std::string err("error :" + *it + " not in r matrix\n");
			std::cerr << err << "\n";
			throw err;
		}
		idx[i] = matidx[*it];
		++i;
	}
	return i;
}

void groupByTime(std::vector<std::pair<int,std::vector<std::string>>> &cohorts, const std::vector<indiv> &ped, bool descend = false) {
	size_t sz = cohorts.size();
	cohorts.resize(sz + ped.size());
	std::unordered_map<int, unsigned int> vecidx;
	unsigned int c = 0;
	for (const indiv &i : ped) {
		if (!isnan(i.cohort)) {
			if (vecidx.find(i.cohort) == vecidx.end()) {
				vecidx[i.cohort] = c;
				cohorts[vecidx[i.cohort]].first = i.cohort;
				//std::cout << cohorts[vecidx[i.cohort]].first << "\n";
				size_t prevsize = cohorts[vecidx[i.cohort]].second.size();
				cohorts[vecidx[i.cohort]].second.reserve(prevsize + ped.size());
				++c;
			}
			cohorts[vecidx[i.cohort]].second.push_back(i.id());
		}
	}
	cohorts.resize(c);

	for (unsigned int i = 0; i<cohorts.size(); ++i) {
		cohorts[i].second.shrink_to_fit();
	}

	// sort vector by time
	std::sort(cohorts.begin(), cohorts.end(), timeSort(descend));

/*
	// debug printing
	for (unsigned int i = 0; i<cohorts.size(); ++i) {
		std::cout << cohorts[i].first << "\n";
		for (std::vector<std::string>::const_iterator it = (cohorts[i].second).begin(); it != (cohorts[i].second).end(); ++it) {
			std::cout << *it;
			if (it == (cohorts[i].second).end()-1) {
				std::cout << "\n";
			} else {
				std::cout << "\t";
			}
		}
	}

*/
}

double sumRelatedness (const std::string &ancid, const std::vector<std::vector<indiv*>> &lineage, Matrix<double> &rmat, std::unordered_map<std::string, unsigned int> &matidx) {
	static std::vector<indiv*>::const_iterator ind_iter;
	static std::unordered_map<std::string, int> seen;
	seen.clear();
	seen.reserve(rmat.rown());
	double ind_total = 0;
	unsigned int gen = 0;
	std::stringstream ss;

	if (matidx.find(ancid) == matidx.end()) {
		std::string miserr_row("error: " + ancid + " not in relatedness matrix\n");
		std::cerr << miserr_row;
		throw miserr_row;
	}

	while (lineage[gen].size() > 0) {
		for (ind_iter = lineage[gen].begin(); ind_iter != lineage[gen].end(); ++ind_iter) {
			if (matidx.find(ancid) == matidx.end()) {
				std::string miserr_col("error: " + (*ind_iter)->id() + " not in relatedness matrix\n");
				std::cerr << miserr_col;
				throw miserr_col;
			}
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
				/*
				* relationship has not been recorded yet, so do so. Relantionships can get logged twice if an descendant (ind2) produces
				* offspring with their parent (ind1) because the offspring gets logged in lineage vector for both ind1 and ind2. This prevents
				* dobule counting.
				*/
				if ((*ind_iter)->group) {
					ind_total += rmat[matidx[ancid]][matidx[(*ind_iter)->id()]];
					//std::cout << "anc: " << ancid << "\t" << (*ind_iter)->id() << "\t" << (*ind_iter)->cohort << "\n"; // debug
				}
				seen.insert({elem,1});
			}
		}
		++gen;
	}
	return ind_total;
}

void collectDescendants (const std::string ancid, std::vector<indiv> &ped, std::unordered_map<std::string, unsigned int> &pedidx,
   std::vector<std::vector<indiv*>> &lineage, std::unordered_map<std::string, bool>* ex = NULL) {
	lineage.clear();
	lineage.resize(ped.size()+1); // this is maximum possible generations + 1
	unsigned int gen = 0;
	lineage[gen].push_back(&ped[pedidx[ancid]]);
	if (ex) (*ex)[ancid] = 1;
	size_t noffspring = lineage[0][0]->offspring.size();
		while (noffspring > 0) {
			size_t offgen = gen+1;
			noffspring = 0;
			for (std::vector<indiv*>::iterator ind_iter = lineage[gen].begin(); ind_iter != lineage[gen].end(); ++ind_iter) {
				for (std::vector<std::string>::iterator offspring_iter = (*ind_iter)->offspring.begin(); offspring_iter != (*ind_iter)->offspring.end(); ++offspring_iter) {
					if (ex == NULL || ex->find(*offspring_iter) == ex->end()) {
						lineage[offgen].push_back(&ped[pedidx[*offspring_iter]]);
						noffspring++;
					}
				}
			}
			++gen;
		}
}

int hunterStat (std::vector<indiv> &ped, std::unordered_map<std::string, unsigned int> &pedidx, Matrix<double> &rmat,
   std::unordered_map<std::string, unsigned int> &matidx, std::vector<std::string> &matids, std::string &outprefix, std::vector<std::string> *anc,
   const bool &t2_only, int draw, std::vector<std::string>* cohort, const std::vector<std::string> &fields, const float &time2, const bool max_norm) {
	int rv  = 0;

	if (anc->size() < 1) {
		std::string empty_anc_err("Hunter stat requires passing ancestral individuals with --anc\n");
		std::cerr << empty_anc_err;
		throw empty_anc_err;
	}

	if (t2_only && isnan(time2)) {
		std::string t2_dep_err("--t2_only requires --time2\n");
		std::cerr << t2_dep_err;
		throw t2_dep_err;
	}

	// open output streams
	std::ofstream outstream((outprefix + ".pedstat1").c_str(), std::ios_base::out);

	std::ofstream drawstream;
	if (draw) drawstream.open((outprefix + ".topo").c_str(), std::ios_base::out);

	// check for input fields
	//for (const std::string &f : fields) {
	//	std::cout << f << "\n"; // debug
	//}
	bool cohort_info = (std::find(fields.begin(), fields.end(), "cohort") != fields.end()) ? true : false;
	bool survive_info = (std::find(fields.begin(), fields.end(), "cohort_last") != fields.end()) ? true : false;

	std::vector<std::pair<int,std::vector<std::string>>> time_cohorts;
	std::unordered_map<int, unsigned int> timeidx;
	if (cohort_info) {
		groupByTime(time_cohorts, ped, true);
 		unsigned int k = 0;
		for (const std::pair<int,std::vector<std::string>> &t : time_cohorts) {
			timeidx[t.first] = k;
			++k;
		}
	}

	// check for necessary normalization input
	if (max_norm) {
		if (!cohort_info) {
			std::string cohort_err("error: Ancestral population normalization requires 'COHORT' in input ped\n");
			std::cerr << cohort_err;
			throw cohort_err;
		}
		if (!t2_only && !survive_info) {
			std::string cohort_last_err("error: Ancestral population normalization requires requires 'COHORT_LAST' in input ped\n");
			std::cerr << cohort_last_err;
			throw cohort_last_err;
		}
	}

	// mask individuals from statistic calculation
	std::vector<indiv>::iterator pediter;

	if (!isnan(time2)) {
		for (pediter = ped.begin(); pediter != ped.end(); ++pediter) {
			if (t2_only) {
				// mask individuals not born in time2 or with missing cohort info
				if (isnan(pediter->cohort) || pediter->cohort != time2) pediter->group = 0;
			} else {
				// mask individuals not alive at time of focal descendant cohort or with missing survival information
				if (isnan(pediter->cohort) || isnan(pediter->cohort_last) || pediter->cohort > time2  || pediter->cohort_last < time2) pediter->group = 0;
			}
		}
	}


	float maxcohort = NAN;
	float minanc = NAN;
	if (!cohort->empty()) {
		for (pediter = ped.begin(); pediter != ped.end(); ++pediter) {
			// mask individuals not in specified cohort (and record cohort of oldest individual)
			if (std::find(cohort->begin(), cohort->end(), pediter->id()) == cohort->end()) {
				pediter->group = 0;
			} else if (cohort_info) {
				if (!isnan(pediter->cohort) && (isnan(maxcohort) || pediter->cohort > maxcohort)) maxcohort = pediter->cohort;
			}
		}
	}

	// process pedigree

	std::unordered_map<std::string, double> ind_stats;
	ind_stats.reserve(anc->size());
	double pop_total = 0;
	std::vector<std::vector<indiv*>> lineage;
	std::vector<indiv*>::iterator ind_iter;
	std::vector<std::string>::iterator offspring_iter;

	// calculate expected genetic contribution for each focal ancestor
	std::unordered_map<float, int> seen_times;
	seen_times.reserve(time_cohorts.size());

	for (std::string const &ancid : *anc) {
		float anc_age = ped[pedidx[ancid]].cohort;
		if (!isnan(anc_age)) {
			seen_times[anc_age] = 1; // used to track the different year classes represented among focal ancestors
			if (isnan(minanc) || anc_age < minanc) minanc = anc_age;
		}

		// find all direct descedants of ancestor
		collectDescendants(ancid, ped, pedidx, lineage);

		// sum contributions of ancestor down its direct line of descent
		try {
			double ind_total = sumRelatedness(ancid, lineage, rmat, matidx);
			ind_stats[ancid] = ind_total;
			pop_total += ind_total;
		} catch (const std::string &err) {
			throw;
		}

		// draw lineage pedigree for ancestor
		if (draw) {
			unsigned int gen = 0;
			drawstream << "## LINEAGE " << ancid << " ##\n";
			while (lineage[gen+1].size() > 0) {
				std::string genstr("");
				if (lineage[gen+1].size() > 0) drawstream << gen+1 << ": ";
				for (ind_iter = lineage[gen].begin(); ind_iter != lineage[gen].end(); ++ind_iter) {
					if ((*ind_iter)->offspring.size() == 0) continue;
						genstr += "(" + (*ind_iter)->id() + ")";
					for (offspring_iter = (*ind_iter)->offspring.begin(); offspring_iter != (*ind_iter)->offspring.end(); ++offspring_iter) {
						genstr += " " + *offspring_iter;
					}
					genstr += " ## ";
				}
				drawstream << genstr.substr(0,genstr.length()-4) << "\n";
				++gen;
			}
		}
	}

	if (!isnan(maxcohort) && !isnan(minanc) && minanc > maxcohort) {
		std::cerr << "warning: At least some individuals in cohort were born before ancestors\n";
	}

	// calculate contributions by cohort time for Hunter stat normalization
	std::unordered_map<int, double> time_norm;

	if (max_norm) {
		// calculate maximum ancestral contribution by time period
		std::vector<std::pair<int,std::vector<std::string>>>::const_iterator past_iter;

		// go through time periods found for all focal ancestors
		for (std::unordered_map<float, int>::const_iterator time_iter = seen_times.begin(); time_iter != seen_times.end(); ++time_iter) {
			time_norm[time_iter->first] = 0;
			if (time_iter->first > time2) continue; // ancestors born after focal descendant cohort cannot contribute to it
			double* total_ptr = &time_norm[time_iter->first];
			std::unordered_map<std::string, bool> seen_anc;
			seen_anc.reserve(ped.size());

			// starting at the time of ancestor's birth work backwards through time (this focuses on individuals that are extant at the time of ancestor's birth)
			for (past_iter = time_cohorts.begin()+timeidx[time_iter->first]; past_iter != time_cohorts.end(); ++past_iter) {

				// go through all extant individuals for current time period
				for (const std::string &id : past_iter->second) {
					indiv* ancptr = &ped[pedidx[id]];
					if (isnan(ancptr->cohort_last) || ancptr->cohort_last < time_iter->first) continue; // died before ancestor birth or unknown survival

					// Check if the contribution of this focal period ancestor was already calculated (no need to do this twice)
					if (past_iter->first == time_iter->first && ancptr->cohort == time_iter->first) {
						*total_ptr += ind_stats[id];
						seen_anc[id] = 1;
						continue;
					}

					// collect all direct descendants of focal individual
					collectDescendants(id, ped, pedidx, lineage, &seen_anc);

					// sum contribution of focal individual down its direct line of descent
					try {
						*total_ptr += sumRelatedness(id, lineage, rmat, matidx);
					} catch (const std::string &err) {
						throw;
					}
					seen_anc[id] = 1;
				}
			}
		}

	}

	// output statistics
	if (pop_total == 0) std::cerr << "warning: Total focal ancestral contribution is zero\n";

	std::string header("ID\tN_GENOME_COPIES\tP_ANC_FOCAL");
	if (max_norm) {
		header += "\tP_ANC_MAX\t";
		std::vector<int> t2vec;
		t2vec.reserve(time_norm.size());
		for (const auto &kv : time_norm) {
			t2vec.push_back(kv.first);
		}
		std::cerr << "\n# Maximum population contribution to descedants by ancestral cohort\n";
		std::sort(t2vec.begin(),t2vec.end());
		for (const int &t2 : t2vec) {
			std::cerr << t2 << ": " << time_norm[t2] << "\n";
		}
		std::cerr << "\n";
	}
	outstream << header << "\n";
	for (std::string const &ancid : *anc) {
		double ngenome = ind_stats[ancid];
		double proportion_focal_anc = pop_total > 0 ? ngenome/pop_total : NAN;
		outstream << ancid << "\t" << ngenome << "\t" << proportion_focal_anc;
		if (max_norm) {
			indiv* indptr = &ped[pedidx[ancid]];
			double proportion_anc_max = (!isnan(indptr->cohort) && time_norm[indptr->cohort] > 0) ? ngenome/time_norm[indptr->cohort] : NAN;
			outstream << "\t" << proportion_anc_max;
		}
		outstream << "\n";
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

bool pairDescendSort(const std::pair<std::string, double> &p1, const std::pair<std::string, double> &p2) {
	return (p1.second > p2.second);
}

void sortedMatrix (std::vector<sdpvec>& rank_vec, const std::vector<std::string>* anc, const std::vector<std::string>* cohort, Matrix<double> &rmat, vecmap& matidx) {
	if (!rank_vec.empty()) rank_vec.clear();
	rank_vec.resize(cohort->size());

	// ancestor and cohort indices in relatedness matrix
	unsigned int anc_idx [anc->size()];
	unsigned int cohort_idx [cohort->size()];
	try {
		findMatIndex(anc_idx, matidx, anc);
		findMatIndex(cohort_idx, matidx, cohort);
	} catch (const std::string &ex) {
		throw;
	}

	for (unsigned int i = 0; i<rank_vec.size(); ++i) {
		(rank_vec[i]).reserve(anc->size());
		unsigned int j = 0;
		for (const std::string& ancid : *anc) {
			rank_vec[i].push_back(std::make_pair(ancid, rmat[cohort_idx[i]][anc_idx[j]]));
			++j;
		}
		// sort ancestors based on ascending r
		std::sort(rank_vec[i].begin(), rank_vec[i].end(), pairDescendSort);

		// debug print
		//std::cout << (*cohort)[i] << "\n";
		//for (sdpvec::iterator it = rank_vec[i].begin(); it != rank_vec[i].end(); ++it) {
		//	std::cout << it->first << ":" << it->second;
		//	if (it == rank_vec[i].end()-1) {
		//		std::cout << "\n";
		//	} else {
		//		std::cout << " ";
		//	}
		//}
	}
}

int rankWtRelate (Matrix<double> &rmat, std::unordered_map<std::string, unsigned int> &matidx, const std::vector<std::string> &matids, const std::string &outprefix,
   const std::vector<std::string> *anc, const std::vector<std::string>* cohort, const double c) {
	// want to compute stat: 1/ncohort * sum_i_ncohort(1/rank_i*r_i)
	int rv = 0;

	if (anc->empty()) {
		std::cerr << "error: Mosaic stat 2 requires passing ancestral IDs with --anc\n";
		return -1;
	}

	// open output stream
	std::ofstream stat_outstream((outprefix + ".skewstat3").c_str(), std::ios_base::out);

	// if no cohort IDs, se them to compliment of ancestors in relatedness matrix
	const std::vector<std::string>* cohort_ptr = cohort;
	std::vector<std::string> cohort_local;
	if (cohort_ptr->empty()) {
		std::cerr << "Treating all non-ancestral individuals in relatedness matrix as the focal cohort\n";
		cohort_local.reserve(rmat.rown()-anc->size());
		for (const std::string& id : matids) {
			if (std::find(anc->begin(), anc->end(), id) == anc->end()) cohort_local.push_back(id);
		}
		cohort_ptr = &cohort_local;
	}
	std::cerr << "Focal cohort size: " << cohort_ptr->size() << "\n";

	// for each individual in focal cohort rank ancestors according to relatedness
	std::vector<sdpvec> anc_rank;
	try {
		sortedMatrix(anc_rank, anc, cohort_ptr, rmat, matidx);
	} catch (const std::string & ex) {
		return -1;
	}

	// index array position by ancestor name
	std::unordered_map<std::string, unsigned int> arridx;
	arridx.reserve(anc->size());
	for (unsigned int idx = 0; idx < anc->size(); ++idx) {
		arridx[(*anc)[idx]] = idx;
	}

	// calculate skew statistic
	double skew [anc->size()] = {0};
	unsigned int k = 0;
	for (const sdpvec &cohortind : anc_rank) {
		//std::cout << (*cohort)[k] << ", " << cohortind.size() << "\n"; // debug
		unsigned int i = 1;
		for (const std::pair<std::string, double> &ancind : cohortind) {
			// Find probability of sampling a random ancestor with lower relatedness to cohort individual than the focal ancestor
			double r = ancind.second < c ? 0 : ancind.second;
			size_t nlower = 0;
			for (sdpvec::const_iterator it = cohortind.begin()+i; it != cohortind.end(); ++it) {
				if (it->second < r) ++nlower;
			}
			//std::cout << ancind.first << "\t" << ancind.second << "\t" << nlower << "\n"; // debug
			// add prob_more_related x prob_IBD to ancestor array
			skew[arridx[ancind.first]] += nlower*r;
			++i;
		}
		//++k;
	}

	size_t lowerd = anc->size()-1;
	for (k = 0; k < anc->size(); ++k) {
		skew[k] /= (lowerd*cohort->size()); // converts counts of individuals with lower relatedness into probability and takes average over cohort
	}

	// print skew statistic
	stat_outstream << "ID\tSwtr\n";
	k = 0;
	for (const std::string &ancid : *anc) {
		stat_outstream << ancid << "\t" << skew[k] << "\n";
		++k;
	}

	if (stat_outstream.is_open()) stat_outstream.close();

	return rv;
}

int mosaicStat (Matrix<double> &rmat, std::unordered_map<std::string, unsigned int> &matidx, const std::vector<std::string> &matids, const std::string &outprefix,
   const std::vector<std::string> *anc, const std::vector<std::string>* cohort, const double c) {
	/*
	* calculates skew statistic used in Mosaic FSJ study
	* c: background relatedness
	*/

	int rv = 0;

	if (anc->size() < 1) {
		std::cerr << "error: Mosaic stat requires passing ancestral individuals with --anc\n";
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
	unsigned int cohortidx [rmat.coln()];
	int anc_n = 0, cohort_n = 0;
	try {
		anc_n = findMatIndex(ancidx, matidx, anc);
		cohort_n = findMatIndex(cohortidx, matidx, cohort);
	} catch (const std::string &ex) {
		return -1;
	}

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
	std::ifstream pop_is; // population ID input stream
	std::ifstream cohort_is; // list of IDs to calculate representation among
	float time2 = NAN; // descendant population time point
	bool t2_only = 0; // calculate representation only among individuals born during time2
	bool max_norm = 1; // if 1 normalize contributions using max ancestral cohort contribution as in Hunter etal 2019
	int draw = 0; // draw direct descendent pedigree if 1
	double background_r = 0.0; // background relatedness level used for Mosaic statistic
	double min_r = 0.0; // minimum r value, currently not used

	// parse arguments
	try {
		if ((rv = parseArgs(argc, argv, ped_is, rmat_is, outprefix, pedstat, skewstat, anc_is, pop_is, t2_only, draw, cohort_is, background_r, min_r, time2, max_norm))) {
			if (rv == 1) rv = 0;
			return rv;
		}
	} catch (const std::string &err) {
		std::cerr << err;
		return -1;
	}

	// parse ID lists
	std::vector<std::string> anc; // vector of ancestral individual IDs
	if (anc_is.is_open() && !getIDs(anc_is, &anc)) {
		std::cerr << "error: Read zero ancestral IDs\n";
		return -1;
	}
	if (anc_is.is_open()) anc_is.close();

	std::vector<std::string> cohort; // vector of cohort IDs
	if (cohort_is.is_open() && !getIDs(cohort_is, &cohort)) {
		std::cerr << "error: Read zero cohort individual IDs\n";
		return -1;
	}

	// parse ID maps (e.g. for population)
	std::unordered_map<std::string, std::string> pop; // map relating IDs of individuals to their population
	idMap<std::string>(pop_is, pop);
	if (pop_is.is_open()) pop_is.close();

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
	std::vector<std::string> fields; // stores pedigree fields
	fields.reserve(6);
        std::vector<indiv> ped;
        std::unordered_map<std::string, unsigned int> pedidx;
	if (ped_is.is_open()) {
		if ((rv = readPed(ped_is, &ped, &pedidx, &fields))) return rv;
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
			try {
				hunterStat(ped, pedidx, rmat, matidx, matids, outprefix, &anc, t2_only, draw, &cohort, fields, time2, max_norm);
			} catch (const std::string & err) {
				return -1;
			}
		}
	}

	for (iter = skewstat.begin(); iter != skewstat.end(); ++iter) {
		if (*iter == 1) {
			std::cerr << "Calculating Mosaic skew statistic\n";
			try {
				if ((rv = mosaicStat(rmat, matidx, matids, outprefix, &anc, &cohort, background_r))) return rv;
			} catch (const std::string & err) {
				return -1;
			}
		}
		if (*iter == 2) {
			std::cerr << "Calculating matrix proportion skew statistic\n";
			if ((rv = matPstat (rmat, matidx, matids, outprefix, &anc, &cohort, min_r))) return rv;
		}
		if (*iter == 3) {
			std::cerr << "Calculating average rank-weighted relatedness\n";
			if ((rv = rankWtRelate(rmat, matidx, matids, outprefix, &anc, &cohort, background_r))) return rv;
		}
	}

	return rv;
}
