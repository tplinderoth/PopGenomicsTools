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

// define classes and class functions

// INDIV CLASS
class indiv
{
public:
	// functions
	indiv(std::string);
	std::string& setid(std::string);
	std::string id () const;
	char sex; // M=male, F=female, '*'=missing
	int cohort; // year class of individual, -999 = missing
	std::vector<std::string> offspring; // vector of IDs of offspring
private:
	// members
	std::string _id; // individual ID
};

indiv::indiv (std::string id = "*")
        : sex('N'),
          cohort(-999),
	 _id(id)
{}

std::string& indiv::setid (std::string id = "") {
	if (!id.empty()) _id = id;
	return _id;
}

std::string indiv::id() const {return _id;}

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

void helpinfo() {
	int w1 = 14;
	int w2 = 8;
	int w3 = 2;

	std::cout << "\nrelateStats [options] <input>\n"
	<< "\nInput:\n"
	<< std::setw(w1) << std::left << "--ped" << std::setw(w2) << std::left << "<FILE> ped-format file\n"
	<< std::setw(w1) << std::left << "--rmat" << std::setw(w2) << std::left << "<FILE> Relatedness matrix\n"
	<< std::setw(w1) << std::left << "--out" << std::setw(w2) << std::left << "<STRING> Output name prefix\n"

	<< "\nOptions\n"
	<< std::setw(w1) << std::left << "\n--pedStat" << std::setw(w2) << std::left << "<INT> Pedigree-based statistics\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "1: Expected genetic contribution from Hunter etal 2019\n"
	<< std::setw(w1) << std::left << "\n--skewStat" << std::setw(2) << std::left << "<INT> Genetic skew statistics\n"
	<< std::setw(w3) << std::left << "" << std::setw(w2) << std::left << "1: Mosaic FSJ skew statistic (in development)\n"

	<< "\nNotes:\n"
	<<"* Assumes first row of relatdness matrix contains individual IDs\n\n";
}

int parseArgs (int argc, char** argv, std::ifstream &ped_is, std::ifstream &rmat_is, std::string &outprefix, std::vector <int> &pedstat, std::vector <int> &skewstat, std::ifstream &anc_is,
	std::ifstream &pop_is, int &mincohort, int &maxcohort) {
	int rv = 0;
	int argpos = 1;
	if (argc < 2 || strcmp(argv[argpos], "-h") == 0 || strcmp(argv[argpos],"--help") == 0) {
		helpinfo();
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
		} else if (strcmp(argv[argpos],"--pedStat") == 0) {
			int s = atoi(argv[argpos+1]);
			if (s > 1) { // change if adding more pedigree stat options
				std::cerr << "Invalid --pedStat option " << s << "\n";
				return -1;
			}
			pedstat.push_back(s);
		} else if (strcmp(argv[argpos],"--skewStat") == 0) {
			int s = atoi(argv[argpos+1]);
			if (s > 1) { // change if adding more skew stat options
				std::cerr << "Invalid --skewStat option " << s << "\n";
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
		} else if (strcmp(argv[argpos],"--mincohort") == 0) {
			mincohort = atoi(argv[argpos+1]);
		} else if (strcmp(argv[argpos],"--maxcohort") == 0) {
			maxcohort = atoi(argv[argpos+1]);
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

int hunterStat (std::ifstream &ped_is, std::ifstream &rmat_is, std::string &outprefix, std::vector<std::string>* pop,
   std::vector<std::string> *anc, int mincohort, int maxcohort) {
	int rv  = 0;

	// open output stream
	std::ofstream outstream((outprefix + ".pedstat1").c_str(), std::ios_base::out);

	// read in pedigree
	std::vector<indiv> ped;
	std::unordered_map<std::string, unsigned int> pedidx;
	if ((rv = readPed(ped_is, &ped, &pedidx))) return rv;

	// check ped file parsing (debug)
	/*
	std::cerr << "ped vector size: " << ped.size() << "\nped index vector size: " << pedidx.size() << "\n";
	for (std::vector<indiv>::iterator it = ped.begin(); it != ped.end(); it++) {
		std::cout << it->id() << ":";
		for (std::vector<std::string>::iterator child_iter = (*it).offspring.begin(); child_iter != (*it).offspring.end(); child_iter++) {
			std::cout << " " << *child_iter;
		}
		std::cout << "\n";
	}
	*/

	// read in relatedness matrix
	std::unordered_map<std::string, unsigned int> matidx;
	std::vector<std::string> matids;
	Matrix<double> rmat;
	if ((rv = readrMat(rmat_is, &matidx, &matids, rmat))) return rv;
	std::cerr << "Relatedness matrix dimensions: " << rmat.rown() << " x " << rmat.coln() << "\n";

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

	// For for pedigree individuals not present in the relatedness matrix
	for (vecmap::iterator it = pedidx.begin(); it != pedidx.end(); ++it) {
		if (matidx.find(it->first) == matidx.end()) {
			std::cerr << "warning: " << it->first << " missing from relatedness matrix\n";
		}
	}

	// calculate statistics

	// close output stream
	outstream.close();

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
	int maxcohort = -999, mincohort = -999;

	// parse arguments
	if ((rv = parseArgs(argc, argv, ped_is, rmat_is, outprefix, pedstat, skewstat, anc_is, pop_is, mincohort, maxcohort))) {
		if (rv == 1) rv = 0;
		return rv;
	}

	// parse ID lists
	std::vector<std::string> anc; // vector of ancestral individual IDs
	if (anc_is && !getIDs(anc_is, &anc)) {
		std::cerr << "error: Read zero ancestral IDs\n";
		return -1;
	}
	std::vector<std::string> pop; // vector of individual IDs in focal population
	if (pop_is && !getIDs(pop_is, &pop)) {
		std::cerr << "error: Read zero focal population IDs\n";
		return -1;
	}

	// calculate statistics
	std::vector<int>::iterator iter;

	for (iter = pedstat.begin(); iter != pedstat.end(); iter++) {
		if (*iter == 1) {
			std::cerr << "Calculating Hunter expected contribution\n";
			if ((rv = hunterStat(ped_is, rmat_is, outprefix, &pop, &anc, mincohort, maxcohort))) return rv;
		}
	}

	// close input streams
	ped_is.close();
	rmat_is.close();

	return rv;
}
