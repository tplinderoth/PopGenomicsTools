PopGenomicsTools
================

### A collection of small programs for population genetic analyses

Contact: Tyler Linderoth, lindero1@msu.edu
__________________________________________________________________

To compile all programs `cd` into the PopGenomicsTools directory and then type
```
make
```

### relateStats

Calculate expected genetic contribution given pedigree input or other genetic contribution statistics based on relatedness matrices.

#### Compile
`cd` into the PopGenomicsTools directory and then compile using
```
g++ -Wall -O3 -o relateStats relateStats.cpp
```

#### Runing relateStats

Help can be printed by running relateStats without arguments, `./relateStats`, or using `./relateStat --help`

Which will print
```
relateStats <input>

Possible Inputs:
--pedstat       <INT> Pedigree-based statistics
  1: Expected genetic contribution from Hunter etal 2019

--skewstat      <INT> Genetic skew statistics
  1: Ranked relatedness among relatives
  2: Relatedness matrix proportion
  3: Average rank-weighted relatedness

--out           <STRING> Output name prefix
--ped           <FILE> ped format file
--rmat          <FILE> Relatedness matrix
--pop           <FILE> Two-column, tab-delimited file specifying (1) individual ID & (2) population ID
--anc           <FILE> List of ancestor IDs
--cohort        <FILE> List of individual IDs to restrict genetic representation analyses to
--time2         <INT> Sets descendant population to all extant individuals between ancestor's cohort and INT
--t2_only       Restrict representation to individuals born during --time2
--max_norm      <0|1> Normalize contributions with respect to maximum ancestral cohort contribution if 1 [1]
--draw          Output direct descendent pedigrees
--background_r  <FLOAT> Background relatedness for skewstat 1 and 3 [0]
--min_r         <FLOAT> Consider r values < FLOAT 0 for skewstat 2 [0]

Pedigree statistics:
--pedstat --out --ped --rmat --anc [--time2] [--t2_only] [--cohort] [--max_norm] [--draw]

Skew statistics:
--skewstat --out --rmat --anc [--cohort] [--background_r] [--min_r]

Notes:
* Assumes first row of relatedness matrix contains individual IDs
```
#### Inputs

There are two main types of inputs, a **--ped** pedigree file and/or a **--rmat** relatedness matrix.

**--ped** : This is a tab-delimited ped format file with required columns (1) individual ID, (2) parent 1 ID, (3) parent 2 ID. A header is assumed and required. The names 
for the first three columns can be anything. Additional columns (in any order) can contain 'sex' ('M','male','F','female', and '*' for missing), an integer-valued 'cohort' 
(e.g. 0, 1, 2, 2003, 2004, 2021), and integer-valued 'cohort_last' (format same as cohort and indicates the last time point at which an individual is alive). Additional columns 
must be named 'sex', 'cohort', and 'cohort_last' if used. A few line example:

```
ID      SIRE_ID DAM_ID  SEX     COHORT COHORT_LAST
SR-BK   LK-WS   PK-BS   FEMALE  2004   2009
SR-PK   LK-WS   PK-BS   FEMALE  2004   2016
SR-OK   LK-WS   PK-BS   *       2004   2004
SR-LK   BK-GS   BK-BS   *       2004   *
SR-FK   KRS-    GK-BS   *       2007   2008
S-WK    *       *       MALE    *      2005
OK-FS   *       *       FEMALE  2008   2010
```

**--rmat** : This is a whitespace-delimited, square matrix of pairwise relatedness values. This file must contain a header (as the first row) with the IDs for the 
individuals represented by the rows/colums. Diagonal entries must be a numeric value, typically 1. An example of a small relatedness matrix:
```
APZ-F     PGZ-F     ROZ-F     ORZ-F     OOZ-F
1         0.079762  2.5e-05   0.019089  0.031298
0.079762  1         0.000388  0.037417  0.174358
2.5e-05   0.000388  1         4e-06     1e-06
0.019089  0.037417  4e-06     1         0.047881
0.031298  0.174358  1e-06     0.047881  1

``` 

**--pop** : Two-column, tab-delimimted file with columns specifying (1) individual ID, and (2) population name. Each row is an individual. Example:
```
WSA-K	Site1
K-RSW	Site18
SR-PK	Site13
WSF-K	Site13
K-GSW	Site13
SRL-K	Texaco
SAG-K	Texaco
```

**--anc**, **--cohort** : Each of these arguments are used to restrict analyses to subsets of individuals and take as input a one-column file of individuals IDs. Each row is one individual. Example:
```
SR-PK
SR-FK
S-WK
OK-FS
WSF-K
WSA-K
WK-SB
K-GSW
K-RSW
```
Both **--pedstat** and **--skewstat** calculate the genetic representation of each ancestor in the **--anc** list among the 'focal cohort' (usually a descendant group
of individuals) specified in the file passed to **--cohort**.

**--pedstat** : This takes an INT argument and performs analyses based on an input pedigree. A description of analyses (INT arguments) follows.

1 : Calculate the genetic contribution of each individual in **--anc** to the focal cohort as in 
[Hunter etal 2019](https://academic.oup.com/jhered/article/110/4/433/5525396). The focal cohort is specified with **--cohort** or **--time2**. 
If focal cohort IDs are not supplied assumes all non-ancestral individuals in relatedness matrix are in the focal cohort.
**Requirements: --rmat, --anc, --out**.

**--skewstat** : This takes an INT argument and performs analyses based on a relatedness matrix. A description of analyses (INT arguments) follows.

1 : Calculates genetic contribution based on ranked relatedness among cohort relatives. Specifically, for each ancestral individual in **--anc** this quantifies their 
genetic representation in the focal cohort based on the proportion of individuals that they are related to and their rank among all potential ancestors to these individuals. 
**Requirements: --rmat, --anc, --out**.

2 : Calculates genetic representation for an ancestor based on the ratio of their relatedness with the focal cohort to the total relatedness between all ancestors and the 
focal cohort. Also calculates the proportion of pairwise {ancestor, focal cohort individual} comparisons for which the relatedness involving 
the focal ancestor is higher. If neither **--anc** or **--cohort** are supplied statistics are calculated based on all pairwise comparisons in the relatedness matrix. If 
only **--anc** is supplied, assumes all non-ancestral individuals in relatedness matrix are focal cohort individuals. If only **--cohort** is supplied, assumed all other individuals in 
relatedness matrix are ancestors. **Requirements: --rmat, --out**.

3 : Calculates genetic contribution based on rank-weighted relatedness to focal cohort individuals. Specifically, for each ancestral individual in **--anc** this is 
the probability that the ancestor will have the highest ancestral relatedness to a random focal cohort individual and that two alleles drawn from this pair will be 
identical by descent (IBD). **Requirments: --rmat, --anc, --out**.

**--out** : output file name prefix.

Note that multiple analyses can be run in a single call, e.g. `./relateStats --skewstat 1 --skewstat 2 --pedstat 1 ...` is valid.

#### More running options

**--time2** : Treat all individuals alive in the population at the time/cohort value passed this argument as the focal cohort. Requires 'cohort' to be present 
in ped input.

**--t2_only** : Treat only individuals born at **--time2** as the focal cohort.

**--max_norm** : Normalize contribution of each ancestor based on the maximum contribution that any individual alive in the population at the time of the ancestor's birth 
could make to the focal cohort. This is the normalization introduced in [Hunter etal 2019](https://academic.oup.com/jhered/article/110/4/433/5525396). This option applies when 
--pedstat is used and is on by default (1) but can disabled with '0'. Note that if used these normalized values are printed in addition to values normalized with respect to the 
total contribution among the individuals in **--anc**.

**--background_r** : Pairs of individuals with relatedness above this value are considered relatives.

**--min_r** : Sets relatedness values below this level to zero for --skewstat 2 **S<sub>count</sub>** calculation. This can help reduce noise from very low relatedness values.

**--draw** : Specifying this dumps a '.topo' file containing a crude representation of each ancestral lineage pedigree from ancestor down through descedants. This is useful for visualizing 
lines of descent.

#### Output

The following describes the output columns for each analysis.

**--pedstat 1**

.pedstat1<br>
**(1) ID** : Ancestor ID.<br>
**(2) N_GENOME_COPIES** : Expected number of ancestor genome copies in focal cohort.<br>
**(3) P_ANC_FOCAL** : Ancestor's proportion of the total expected genomic copies from all ancestors in **--anc** among the focal cohort.<br>
**(4) P_ANC_MAX** : Ancestor's proportion of the total expected genomic copies in the focal cohort out of the max contribution possible for any individual alive in the population 
when the ancestor was born/came into existence (e.g. through migration). This is the Hunter *et al.* (2019) normalized value. Note that here 'alive in the population' really means present 
in the ped input.<br>

.topo<br>
Pedigree representation of each ancestral lineage.

**--skewstat 1**

.skewstat1<br>
**(1) ID** : Ancestor ID.<br>
**(2) S<sub>c</sub>** : Joint probability that the ancestor is a relative of a random focal cohort individual and two alleles from this pair are IBD.<br>
**(3) S<sub>rank</sub>** : Joint probility that the ancestor is a relative of a random focal cohort individual and their relatedness is higher than for any other ancestor.<br>
**(4) relate_prob** : Probability that the ancestor and a focal cohort individual are related above level **--background_r**.<br>
**(5) avg_r** : Probability that two alleles drawn from the ancestor and a random focal cohort relative are IBD.<br>
**(6) rank_prob** : Probability that the relatedness between the ancestor and a random focal cohort relative is higher than for any other potential ancestor.<br>

.relatives<br>
**(1) ID** : Ancestor ID.<br>
**(2) relatives** : Comma-delimited list of IDs of focal cohort individuals with relatedness to ancestor above **--background_r**. '*' indicates no relatives.<br>

**--skewstat 2**

.skewstat2<br>
**(1) ID** : Ancestor ID.<br>
**(2) S<sub>count</sub>** : Proportion of all pairwise {ancestor, focal cohort individual} comparisons for which the relatedness involving the ancestor is higher.<br>
**(3) S<sub>rsub</sub>** : Proportion of the total relatedness between ancestors and focal cohort indivivduals from this ancestor.<br>

**--skewstat 3**

.skewstat3<br>
**(1) ID** : Ancestor ID.<br>
**(2) S<sub>wtr</sub>** : Joint probability that the ancestor has the highest ancestral relatedness to a random focal cohort individual and alleles drawn from this pair are IBD.
 
__________

### betaAFOutlier.R

Allele frequency difference p-values based on a null distribution generated through simulations under the Balding-Nichols model.

#### Running betaAFOutlier.R

Help

```
betaAFOutlier.R <pop1n,pop2n> <fst> <allele frequency file> <out prefix> [options]

pop1n,pop2n: Comma-separated diploid sample sizes for populations 1 and 2.
fst: Genome-wide average FST value between groups being compared or '9' to estimate FST from input.
allele frequency file: TSV file with columns (1) chromosome, (2) position, (3) pop1 allele frequency (4) pop2 allele frequency. Assumes header.
out prefix: Output file name prefix to which '.dist', '.ancf', '.sigtest', and '.pdf' are appended

Optional Input
--ancmethod <0|1>: Ancestral allele frequency method where 0 (default) uses empirical ancestral frequency prior, 1 uses expected frequency prior
--ancbin <FLOAT in range [0,1]>: Ancestral allele frequency prior bin width [default: 0 (no binning)]
--plotqq: Generate qq-plot of p-values
--seed <INT>: Set a specific seed
```

__________

### selkit

A toolbox for selection analysis from population genetic data.
Available tests:
* Hudson–Kreitman–Aguadé (HKA) test

#### Running selkit

**HKA** help obtained with `./selkit hka`

```
selkit hka [arguments]

Arguments:
-vcf        VCF file to analyze (must be bgzipped and indexed)
-popfile    TSV file with rows having VCF sample name followed by a 0 (outgroup) or 1 (ingroup) population identifier
-rf         TSV file with regions in format CHR POS END to calculate HKA statistic for
-r          Region supplied as a string in format 'chr:from-to' to calculate HKA statistic for
-out        Name of output file
-include    TSV file with regions in format CHR POS END to use for expectation parameter estimation [default: all sites]
-exclude    TSV file with regions in format CHR POS END to exclude for expectation calculation
-ps         Probability of a neutral site being polymorphic within the ingroup
-pd         Probability that a neutral site is fixed between the ingroup and outgroup
-passonly   Use only sites with PASS in VCF FILTER field
-ratesonly  Only calculate the parameters used for HKA expectations

HKA statistic calculated by comparing all ingroup individuals to all outgroup individuals.
If supplying precalculated probabilities -ps and -pd must sum to 1.
```
__________

### dxyWindow

Calculates Dxy statistic in sliding windows across the genome.

#### Running dxyWindow

Help

```
dxyWindow [options] <pop1 maf file> <pop2 maf file>

Options:
-winsize      INT     Window size in base pairs (0 for global calculation) [0]
-stepsize     INT     Number of base pairs to progress window [0]
-minind       INT     Minimum number of individuals in each population with data [1]
-fixedsite    INT     (1) Use fixed number of sites from MAF input for each window (window sizes may vary) or (0) constant window size [0]
-sizefile     FILE    Two-column TSV file with each row having (1) chromsome name (2) chromosome size in base pairs
-skip_missing INT     Do not print windows with zero effective sites if INT=1 [0]

Notes:
* -winsize 1 -stepsize 1 calculates per site dxy
* -sizefile is REQUIRED(!) with -fixedsite 0 (the default)
* Both input MAF files need to have the same chromosomes in the same order
* Assumes SNPs are biallelic across populations
* For global Dxy calculations only columns 4, 5, and 6 below are printed
* Input MAF files can contain all sites (including monomorphic sites) or just variable sites
* -fixedsite 1 -winsize 500 would for example ensure that all windows contain 500 SNPs

Output:
(1) chromosome
(2) Window start
(3) Window end
(4) dxy
(5) number sites in MAF input that were analyzed
(6) number of sites in MAF input that were skipped due to too few individuals
```
__________

### fstWindow

Calculates F<sub>ST</sub> in sliding windows across the genome.

#### Running fstWindow

Help

```
Usage:
fstWindow [ANGSD fst variance component file] [window size (number sites)] [step size (number sites)]
default window size: 1
default step size: 1

Output:
(1) chromosome
(2) window start
(3) window end
(4) window midpoint position
(5) Fst
(6) Number sites in window
```
__________

### hetWindow

Calculates heterozygosity in sliding windows across the genome.

#### Running hetWindow

Help

```
Usage:
hetWindow [genotypes file] [window size (number sites)] [step size (number sites)]
default window size: 1
default step size: 1

Output:
(1) chromosome
(2) window start
(3) window end
(4) window midpoint position
(5) heterozygosity
(6) Number sites in window
```
__________

### ihsWindow

iHS statistic analysis in sliding windows across the genome for selection analyses.

#### Running ihsWindow

```
Usage:
ihsWindow [selscan normalized iHS *.norm file] [options]

Assumes iHS locus ID in format chr*_position

Options:
-winsize INT Window size (bp) [100000]
-cutoff FLOAT Determine fraction of sites with |iHS| > cutoff [2]
-chrlen FILE TSV-file with columns (1) chr (2) chromosome length (bp), and each row is a different chromosome

Output:
(1) chromosome
(2) window start
(3) window stop
(4) most extreme iHS score
(5) extreme iHS position
(6) proportion |iHS| > cutoff
(7) Number SNPs in window
```
__________

### xpehhWindow

Sliding window analysis for selection using the cross-population EHH statistic.

#### Running xpehhWindow

Help

```
Usage:
xpehhWindow <selscan normalized XPEHH *.norm file> <cutoff> [options]

Input file must have locus ID in format chr*_position
cutoff (FLOAT): Calculate proportion of sites with EXPEHH less (if negative) or greater (if positive) than cutoff

Options:
-winsize INT Window size (bp) [100000]
-chrlen FILE TSV-file with columns (1) chr (2) chromosome length (bp), and each row is a different chromosome

Output:
(1) chromosome
(2) window start
(3) window stop
(4) minimum (negative cutoff) or maximum (postive cutoff) XPEHH score
(5) extreme XPEHH position
(6) proportion XPEHH scores > or < cutoff
(7) Number SNPs in window
```
__________

### Other code

* pafAlleles : Used for genomic assembly liftover
* ngsMisc.pl : Toolkit for processing and manipulating NGS data
