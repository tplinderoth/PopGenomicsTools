#!/usr/bin/perl

# hmm2froh.pl

# The calculation of Froh implemented in this script is consistent with how it was defined in the supplemental materials of the BCFtools paper by Narasimhan etal (2016):
# "To characterize the total autozygosity for each population, we counted the total number of callable bases in the autozygous sections and divided that by the total number of 
# callable bases in the genome. We define callable regions as those that lie within the 1000 genomes accessibility mask." This script uses this approach, i.e., it estimates 
# Froh as the fraction of sites in unmasked sections of autozygous intervals out the total number of unmasked sites in the genome.

# Note that Narasimhan etal (2016) also calculate Froh based on "...the number of autozygous sites detected divided by the total number of variant sites in the individual sample, Froh."
# This alternate method relies on SNPs being approximately evenly distributed throughout the genome to give accurate results. This is not the approach implemented in this script.

use warnings;
use strict;
use Getopt::Long;

my $states_file = undef;
my $mask_file = undef;
my $minq = 20;
my $minlength = 1;
my $exclude_file = undef;

die(qq/
hmm2froh.pl <inputs>

Required input:
--states      FILE    Ouput from BCFtools ROH specifying autozygous or non-autozygous site states
--mask        FILE    TSV file of unmasked sites \(callable sites\) in CHR\\tSTART\\tEND format

Optional input:
--minq        FLOAT   Discard autozygous and nonautozygous tracts >= minlength with average quality below FLOAT [$minq]
--minlength   INT     Minimum ROH length [$minlength]
--exclude     FILE    File listing scaffolds \(one per row\) to exclude from analyses.

Notes:
- The input file for --states is the output of BCFtools roh run with output type set to "--output-type sr".
- The mask file should contain regions specifying all callable variable and monomorphic sites in the genome.
- Assumes site states file and mask file are sorted, i.e., scaffolds and sites are in the same order in both files.
- Assumes that all SNPs in the states file are in callable regions of the genome.
- Returns the fraction of unmasked sites in ROH intervals out of the total number of unmasked sites in the genome.

Ouput:
[1] FROH: Proportion of accessible sites in the genome that are autozygous.
[2] ANALYZED_ROH: Number of analyzed ROH.
[3] LONGEST_ROH: Longest ROH that passed quality checks.
[4] TOTAL_ACCESSIBLE_SITES: Number of accessible sites \(--mask\) minus those ignored due to low quality \(--minq\)
[5] TOTAL_NUMBER_ROH: Total number of potential ROH of any length.
[6] QC_FAIL_ROH: Number of ROH above the minimum length cutoff discarded due to low quality.
[7] FAIL_ROH_SITES: Total number of sites in ROH discarded due to low quality.
[8] QC_FAIL_NONAUTO_REGIONS: Number of nonautozygous regions above the minimum length cutoff discarded due to low quality.
[9] FAIL_NONAUTO_SITES: Total number of sites in nonautozygous regions discarded due to low quality.
\n/) if (!@ARGV || scalar @ARGV < 4);

my $rv = GetOptions('states=s' => \$states_file, 'mask=s' => \$mask_file, 'minq=f' => \$minq, 'minlength=i'=> \$minlength, 'exclude=s' => \$exclude_file);
die("--> input error, exiting prematurely\n") unless ($rv);

die("Error: --states is required\n") if ! $states_file;
die("Error: --mask is required\n") if ! $mask_file;
open(my $hmmfh, '<', $states_file) or die("Error: unable to open file of site states, $states_file\n");
open(my $maskfh, '<', $mask_file) or die("Error: unable to open accessability mask file, $mask_file\n");

die("Error: --minq must be >= 0\n") if ($minq < 0);
die("Error: --minlength must be > 0\n") if ($minlength < 1);

# parse scaffolds to ignore
my %exscaff;
if ($exclude_file) {
	open(my $exfh, '<', $exclude_file) or die("Error: unable to open file of scaffolds to ignore, $exclude_file\n");
	while (<$exfh>) {
		chomp;
		$exscaff{$_} = 1;
	}
	close $exfh;
}

# process 

my $total_unmasked = 0; # total number of unmasked sites
my $n_autozygous = 0; # number of unmasked autozygous sites
my $ntracts_keep = 0; # number of retained roh tracts
my $ntracts_total = 0; # total number of roh tracts
my $nroh_qc_fail = 0; # number of roh tracts above minlength, but ignored due to quality cutoffs
my $roh_sites_fail = 0; # number of potential roh sites ignored due to quality cutoffs
my $nhwe_qc_fail = 0; # number of nonautozygous tracts amove minlength, but ignored due to quality cutoffs
my $hwe_sites_fail = 0; # number of nonautozygous sites ignored due to quality cutoffs
my @interval = (); # scaffold, start, end, state [0 = nonautozygous, 1 = autozygous]
my $longest = 0; # longest ROH region
my $qsum;
my $nsnps;
my $maskstr;

# calculate total unmasked genome size and store accessible regions
my @access;
while (<$maskfh>) {
	chomp;
	my @tok = split(/\t/, $_);
	next if exists $exscaff{$tok[0]};
	$total_unmasked += $tok[2] - $tok[1] + 1;
	push @access, [@tok];
}

# initialize interval info
while (<$hmmfh>) {
	if ($_ =~ /^ROH/) {
		chomp;
		my @tok = split(/\t/, $_);
		next if exists $exscaff{$tok[0]};
		@interval = @tok[2,3,3,4];
		$qsum = $tok[5];
		$nsnps = 1;
		$ntracts_total = 1 if $interval[3] == 1;
		last;
	}
}

my $m = -1; # accessible region index
while (<$hmmfh>) {
	next unless $_ =~ /^ROH/;
	chomp;
	my @tok = split(/\t/, $_);
	next if exists $exscaff{$tok[0]};
	if (eof($hmmfh) || $tok[4] != $interval[3] || $tok[2] ne $interval[0]) {
		# hit either
		# (1) the end of the states file
		# (2) a different state
		# (3) a new scaffold

		# process the last interval
		my $reglen = $interval[2] - $interval[1] + 1;
		my $avgq = $qsum/$nsnps;
		if ($reglen >= $minlength) {
			# start on mask region overlapping the previous interval to ensure that mask regions spanning two intervals are accounted for
			$m = $m - 1;
			# catch up the accessible regions to interval in order to figure out overlaps between interval and unmasked regions
			while ($m < $#access && (($access[$m]->[0] ne $interval[0]) || ($access[$m]->[0] eq $interval[0] && $access[$m]->[2] < $interval[1]))) {
				# read through accessible regions of the genome until the end of the accessible region is at or beyond the start of the interval
				$m++;
			}
			# Calculate region overlap and progress the masked regions until no more overlaps are found with the current interval
			my $ovl;
			do {
				$ovl = 0;
				if ($access[$m]->[0] eq $interval[0]) {
					if ($interval[1] >= $access[$m]->[1] && $interval[1] <= $access[$m]->[2] && $interval[2] > $access[$m]->[2]) {
					           #-------------# interval
						#-------# mask
						$ovl = $access[$m]->[2] - $interval[1] + 1;
					} elsif ($interval[1] < $access[$m]->[1] && $interval[2] >= $access[$m]->[1] && $interval[2] <= $access[$m]->[2]) {
						#-----------# interval
						         #-----# mask
						$ovl = $interval[2] - $access[$m]->[1] + 1;
					} elsif ($interval[1] < $access[$m]->[1] && $interval[2] > $access[$m]->[2]) {
						#----------------# interval
						     #------# mask
						$ovl = $access[$m]->[2] - $access[$m]->[1] + 1;
					} elsif ($interval[1] >= $access[$m]->[1] && $interval[2] <= $access[$m]->[2]) {
						      #------# interval
						#-------------------# mask
						$ovl = $interval[2] - $interval[1] + 1;
					}
				}
				die("Error: unexpected overlap, $ovl, comparing regions '@interval' and '@{$access[$m]}'") if $ovl < 0;
				die("Error: accessability mask does not contain all sites in --states file") if ($m == $#access && !$ovl);
				# update the number of autozygous sites in accessible ROH regions
				if ($ovl) {
					if ($avgq >= $minq) {
						$n_autozygous += $ovl if $interval[3] == 1; # update number of accessable autozygous sites
					} else {
						$total_unmasked -= $ovl; # subtract ambiguous accessable regions from the total accessable genome length since they are ignored
						$roh_sites_fail += $ovl if ($interval[3] == 1);
						$hwe_sites_fail += $ovl if ($interval[3] == 0);
					}
					# advance mask
					$m++;
				}
			} until (!$ovl || $m > $#access);
		}

		# update tallies
		if ($interval[3] == 1) {
			# autozygous interval
			if ($reglen >= $minlength) {
				if ($avgq >= $minq) {
					$ntracts_keep++;
					$longest = $reglen if $reglen > $longest;
				} else {
					$nroh_qc_fail++;
				}
			}
		} else {
			# non-autozygous interval
			$nhwe_qc_fail++ if ($reglen >= $minlength && $avgq < $minq);
		}

		# set info for new interval
		@interval = @tok[2,3,3,4];
		$qsum = $tok[5];
		$nsnps = 1;
		$ntracts_total++ if ($interval[3] == 1);
	} elsif ($tok[2] eq $interval[0] && $tok[4] == $interval[3]) {
		# extend current interval
		$interval[2] = $tok[3];
		$qsum += $tok[5];
		$nsnps++;
	}
}
close $hmmfh;

# calculate autozyogosity
print STDERR "Warning: there were no accessible sites in file provided with --mask\n" if $total_unmasked < 1;
my $froh = $total_unmasked > 0 ? $n_autozygous/$total_unmasked : -9;

# print output
print STDOUT "FROH\tANALYZED_ROH\tLONGEST_ROH\tAUTOZYGOUS_SITES\tTOTAL_QCPASS_ACCESSIBLE_SITES\tTOTAL_NUMBER_ROH\tQC_FAIL_ROH\tFAIL_ROH_SITES\tQC_FAIL_NONAUTO_REGIONS\tFAIL_NONAUTO_SITES\n";
print STDOUT "$froh\t$ntracts_keep\t$longest\t$n_autozygous\t$total_unmasked\t$ntracts_total\t$nroh_qc_fail\t$roh_sites_fail\t$nhwe_qc_fail\t$hwe_sites_fail\n";

exit;
