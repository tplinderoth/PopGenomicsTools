/*
* ancAllele.cpp
*
* Use parsimony to determine ancestral allele from genome alignments
* Uses reference fasta and genome variants in paftools call format
*
* Compile: g++ -O3 -o ancAllele ancAllele.cpp
*
*/

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

struct alignidx {
	std::string id;
	std::string file;
	std::map<std::string,unsigned int> offset;
};

void info (const std::string &refname) {
	int w=12;
	std::cerr << "\nInfer ancestral alleles for reference regions covered by one contig per aligned query assembly\n\n"
	<< "ancAllele [options]\n\n"
	<< std::setw(w) << std::left << "-fasta" << "Indexed fasta file (required)\n"
	<< std::setw(w) << std::left << "-varlist" << "File with columns (1) assembly/species name, (2) paftools variant file (required)\n"
	<< std::setw(w) << std::left << "-chrlist" << "File with each row specifying a chromosome/scaffold to include\n"
	<< std::setw(w) << std::left << "-refname" << "Name of reference assembly written to output [" << refname << "]\n"
	<< "\n";
}

int parseArgs (int argc, char** argv, std::string &fasta, std::string &varlist, std::string &chrlist, int &maxdepth, std::string &refname) {
	int n = 1;
	int rv = 0;

	if (argc < 5 || (argc == 2 && strcmp(argv[1], "-help") == 0)) {
		info(refname);
		return 1;
	}

	std::ifstream infile;

	while (n < argc) {
		if (strcmp(argv[n], "-fasta") == 0) {
			fasta = argv[n+1];
			infile.open(fasta.c_str());
			if (!infile) {
				std::cerr << "Could not open FASTA file " << fasta << "\n";
				return -1;
			}
			infile.close();
		} else if (strcmp(argv[n], "-varlist") == 0) {
			varlist = argv[n+1];
			infile.open(varlist.c_str());
			if (!infile) {
				std::cerr << "Could not open variant file list " << varlist << "\n";
				return -1;
			}
			infile.close();
		} else if (strcmp(argv[n],"-chrlist") == 0) {
			chrlist = argv[n+1];
			infile.open(chrlist.c_str());
			if (!infile) {
				std::cerr << "Could not open list of chromosomes to keep " << chrlist << "\n";
				return -1;
			}
			infile.close();
		} else if (strcmp(argv[n], "-maxdepth") == 0) {
			maxdepth = atoi(argv[n]);
			if (maxdepth < 1) {
				std::cerr << "-maxdepth must be >= 1\n";
				return -1;
			}
		} else if (strcmp(argv[n], "-refname") == 0) {
			refname = argv[n+1];
		} else {
			std::cerr << "Uknown argument " << argv[n] << "\n";
			return -1;
		}

		n += 2;
	}

	return rv;
}

int parseFasta (const std::string* fasta, std::ifstream &fasta_stream, std::map<std::string,unsigned int*> &faidx_map) {

	fasta_stream.open(fasta->c_str());
	if (!fasta_stream) {
		std::cerr << "Could not open FASTA " << *fasta << "\n";
		return -1;
	}

	std::string faidxname(*fasta);
        faidxname += ".fai";
	std::ifstream faidx_stream(faidxname.c_str());
	if (!faidx_stream) {
		std::cerr << "Could not open FASTA index file " << faidxname << "\n";
		return -1;
	}

	std::string faidxline;

	while(getline(faidx_stream, faidxline)) {
		std::stringstream ss(faidxline);
		std::string id;
		ss >> id;
		faidx_map[id] = new unsigned int[4];
		for (int i=0; i<4; ++i) {
			ss >> faidx_map[id][i];
		}
	}

	return 0;
}

int buildVarIdx (const std::string* varlist, std::vector<alignidx> &idx) {

	std::ifstream liststream(varlist->c_str());
	if (!liststream) {
		std::cerr << "Could not open alignment variant file list " << *varlist << "\n";
		return -1;
	}

	std::string listline;

	while(getline(liststream, listline)) {
		std::stringstream ss(listline);
		alignidx assembly;
		ss >> assembly.id;
		ss >> assembly.file;
		idx.push_back(assembly);
	}

	for (std::vector<alignidx>::iterator it = idx.begin(); it != idx.end(); ++it) {
		const char* fname = (it->file).c_str();
		std::cerr << "indexing " << fname << "\n";
		std::ifstream varfstream(fname);
		if (!varfstream) {
			std::cerr << "Could not open variant file " << fname << "\n";
			return -1;
		}
		std::string varline, tok, chr;
		unsigned int pos = varfstream.tellg();
		while (getline(varfstream, varline)) {
			std::stringstream ss(varline);
			ss >> tok;
			if (tok == "R") {
				ss >> tok;
				if (chr.empty() || tok != chr) {
					chr = tok;
					(it->offset)[chr] = pos;
				}
			}
			pos = varfstream.tellg();
		}
		varfstream.close();
	}

	return 0;
}

int parseChrList (const std::string* fasta, const std::string* chrlist, std::vector<std::string> &id) {

	int c = 0;
	std::string chrname;

	if (chrlist->empty()) {

		std::string faidxname(*fasta);
        	faidxname += ".fai";
        	std::ifstream faidx_stream(faidxname.c_str());
		if (!faidx_stream) {
                	std::cerr << "Could not open FASTA index file " << faidxname << "\n";
                	return -1;
		}

		std::string faidxline;

		while(getline(faidx_stream, faidxline)) {
			std::stringstream ss(faidxline);
			ss >> chrname;
			id.push_back(chrname);
			++c;
		}

		std::cerr << "No scaffold subset file -> processing all " << c << " scaffolds in reference\n";

	} else {

		std::ifstream scafifs(chrlist->c_str());

		if (!scafifs) {
			std::cerr << "Could not open file listing chromosomes to process " << *chrlist << "\n";
			return -1;
		}

		while(getline(scafifs, chrname)) {
			id.push_back(chrname);
			++c;
		}

		std::cerr << "Processing " << c << " scaffolds from " << *chrlist << "\n";
	}

	return 0;
}

unsigned int maxseq (const std::map<std::string,unsigned int*> &faidx_map) {

	unsigned int max = 0;

	for (std::map<std::string,unsigned int*>::const_iterator it = faidx_map.begin(); it != faidx_map.end(); ++it) {
		if ((it->second)[0] > max)
			max = (it->second)[0];
	}

	return max;
}

void allocateSeq (std::vector<std::string*> &seq, unsigned int seqlen, int nassembly) {
	seq.resize(seqlen);
	for (std::vector<std::string*>::iterator seqiter = seq.begin(); seqiter != seq.end(); ++seqiter) {
		*seqiter = new std::string[nassembly];
	}
}

std::string getAA(const std::string* alleles, int arrsize) {
	// identifies ancestral allele as the most frequent among assemblies

	std::map<std::string,int> counts;

	for (int i=0; i<arrsize; ++i) {
		if (alleles[i] != "N") {
			++counts[alleles[i]];
		}
	}

	if (counts.empty()) return "N";

	int nseen = 0; // number of alleles seen with max frequency
	int maxcount = 0; // frequency of most common allele
	std::string aa;
	for (std::map<std::string,int>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
		if (it->second >= maxcount) {
			nseen = it->second == maxcount ? nseen+1:1;
			maxcount = it->second;
			aa = it->first;
		}
	}

	if (maxcount < 2 || nseen > 1) aa = "N";

	return aa;
}

int doAnc (const std::string &fasta, const std::string &varlist, const std::string &chrlist, const std::string &refname) {

	// infers ancestral alleles for reference regions covered by one query contig per assembly

	int rv = 0;

	// parse fasta and fasta index
	if (fasta.empty()) {
		std::cerr << "Fasta file not found and is required\n";
		return -1;
	}
	std::cerr << "reading in fasta index\n";
	std::ifstream fasta_stream;
	std::map<std::string,unsigned int*> faidx_map;
	if ((rv=parseFasta(&fasta, fasta_stream, faidx_map))) {
		return rv;
	}

	// parse list of alignment variant files
	if (varlist.empty()) {
		std::cerr << "Alignment variant file list not found and is required\n";
		return -1;
	}
	std::cerr << "building variant file index\n";
	std::vector<alignidx> varidx; // stores byte start position of each chromosome for each assembly
	if ((rv=buildVarIdx(&varlist, varidx))) {
		return rv;
	}

	// parse list of chromosomes to process
	std::vector<std::string> chrid;
	if ((rv=parseChrList(&fasta, &chrlist, chrid))) {
		return -1;
	}

	// print header
	std::string outheader("chr\tpos");
	outheader += "\t" + refname;
	for (std::vector<alignidx>::iterator assemblyiter = varidx.begin(); assemblyiter != varidx.end(); ++assemblyiter) {
		outheader +=  "\t" + assemblyiter->id;
	}
	outheader += "\tAA";
	std::cout << outheader << "\n";

	// parse reference fasta and aligned assemblies and find ancestral alleles
	std::cerr << "Parsing alignments and finding parsimonious alleles\n";

	int nassembly = varidx.size() + 1;
	unsigned int maxlen = maxseq(faidx_map);
	std::vector<std::string*> seq;
	allocateSeq(seq, maxlen, nassembly);

	for (std::vector<std::string>::iterator chriter = chrid.begin(); chriter != chrid.end(); ++chriter) {
		std::cerr << *chriter << ": " << refname;

		// parse reference alleles
		unsigned int pos = 0;
		fasta_stream.seekg(faidx_map[*chriter][1]);

		char base;
		while (fasta_stream.get(base) && pos < faidx_map[*chriter][0]) {
			if (base != '\n') {
				seq[pos][0] = toupper(base);
				++pos;
			}
		}

		// parse aligned assembly alleles
		int asmn = 0;
		for (std::vector<alignidx>::iterator assemblyiter = varidx.begin(); assemblyiter != varidx.end(); ++assemblyiter) {
			std::cerr << " " << assemblyiter->id;

			++asmn;
			std::ifstream align_stream((assemblyiter->file).c_str());
			if (!align_stream) {
				std::cerr << "Could not open alignment variant file " << assemblyiter->file << "\n";
				return -1;
			}

			align_stream.seekg((assemblyiter->offset)[*chriter]);

			pos = 0;
			std::string alnline, vtype, chr, ref, alt;
			int depth;
			unsigned int startpos, endpos, regstart = 0, regend = 0;
			double mapq;

			while(getline(align_stream, alnline)) {
				std::stringstream alnss(alnline);
				alnss >> vtype >> chr >> startpos >> endpos;

				if (vtype == "R") {
				// new alignment region covered by one query contig
					if (chr == *chriter) {
						pos = regstart = startpos;
						regend = endpos;
					} else {
						// end of chromosome
						break;
					}
				} else {
				// postions within aligned region
					alnss >> depth >> mapq >> ref >> alt;

					if (pos >= regstart && pos < regend) {
					// ensures that parsed positions are within region covered by one query contig

						while (pos < startpos) {
							seq[pos][asmn] = seq[pos][0];
							++pos;
						}

						if (ref == "-") {
						// insertion to reference
							int insertsz = alt.size();
							seq[pos-1][asmn].reserve(insertsz+1);
							for (int i=0; i<insertsz; ++i) {
								seq[pos-1][asmn].push_back(toupper(alt[i]));
							}
						} else if (alt == "-") {
						// deletion from reference
							seq[pos][asmn] = seq[pos][asmn];
							while (pos < endpos) {
								seq[pos][asmn] = "*";
								++pos;
							}
						} else {
						// SNP
							if (ref.size() != 1 || alt.size() != 1) {
								std::cerr << "Problematic entry at " << chr << " " << pos+1 << " in " << assemblyiter->file << "\n";
								return -1;
							}
							seq[pos][asmn] = toupper(alt[0]);
							++pos;
						}

					}
				}
			}

			// add alleles from last aligned region
			while (pos < regend) {
				seq[pos][asmn] = seq[pos][0];
				++pos;
			}

			align_stream.close();

		}

		std::cerr << "\n";

		// find parsimonious ancestral allele
		for (unsigned int site = 0; site < faidx_map[*chriter][0]; ++site) {
			std::cout << *chriter << "\t" << site+1;

			for (int assem = 0; assem < nassembly; ++assem) {
				if (seq[site][assem].empty()) seq[site][assem] = "N";
				std::cout << "\t" << seq[site][assem];
			}

			std::string aa = getAA(seq[site], nassembly);
			std::cout << "\t" << aa;

			std::cout << "\n";

			// clear alleles in sequence structure for next chromosome
			for (int assem = 0; assem < nassembly; ++assem) {
				seq[site][assem].clear();
			}
		}

	}

	// clean up allocated memory

	for (std::vector<std::string*>::iterator seqiter = seq.begin(); seqiter != seq.end(); ++seqiter) {
		delete [] *seqiter;
	}

	for (std::map<std::string,unsigned int*>::iterator mapiter = faidx_map.begin(); mapiter != faidx_map.end(); ++mapiter) {
		delete [] mapiter->second;
	}

	return 0;
}

int main(int argc, char** argv) {

	int rv = 0;

	std::string fasta;
	std::string varlist;
	std::string chrlist;
	int maxdepth = 1; // not using this but leaving it for possible future use
	std::string refname("ref");

	if ((rv=parseArgs(argc, argv, fasta, varlist, chrlist, maxdepth, refname))) {
		return rv;
	}

	if ((rv = doAnc(fasta, varlist, chrlist, refname))) {
		std::cerr << "doAnc() error --> exiting\n";
	}

	return rv;
}
