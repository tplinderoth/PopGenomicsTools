/*
* selkit.h
*
* Toolkit for popgen analyses
* Requires bcftools (http://samtools.github.io/bcftools/bcftools.html) to be installed and in user's PATH
*/

# ifndef _SELKIT_H_
# define _SELKIT_H_

#include <vector>
#include "StatDist.h"

# define COMPRECISION 1e-7 // max precision for comparing floats

// General functions
int main (int argc, char** argv);
void mainInfo ();
int bcftools_check ();
bool fexists(const char* filename);

// HKA test functions
int hka (int argc, char** argv);
double* testRegion (const std::string &vcf, const std::string &region, const std::vector<int> &popmap, int passonly, double ps0, double ps1, double pd01, double* stats);
double hkaGOF (double s0, double s1, double d01, double es0, double es1, double ed01);
int expectedParams(const std::string &vcf, const std::string &exfile, const std::vector<int> &popmap, int passonly, double* ps0, double* ps1, double* pd01);
unsigned long* countVarPatterns (const std::string &cmd, const std::vector<int> &popmap, int passonly, unsigned long *counts);
int countAlt(const std::string &geno);
int pop2idx_sub(const std::vector<int> &full, std::vector<int>* subset, int id);
int indexPops (const std::string &vcf, const std::string &popfile, std::vector<int>* popmap);
int hkaArgs (int argc, char** argv, std::string &vcf, std::string &popfile, std::string &rf, std::string &out, std::string &exfile, int &passonly, double &ps0, double &ps1, double &pd01, int &preprob);
void hkaInfo ();

# endif /* SELKIT_H_ */
