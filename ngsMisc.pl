#!/usr/bin/perl

# ngsMisc.pl

use warnings;
use strict;
use Getopt::Long;

my $version = '1.0.0';

die(qq/
ngsMisc.pl version $version
Toolkit for processing and manipulating NGS data

Usage:
ngsMisc.pl [command]

Commands:
avgDistMat   Collapse pairwise distance matrix into matrix of averages
\n/) if (!@ARGV || $ARGV[0] eq '--help' || $ARGV[0] eq 'h');

my $command = shift @ARGV;
if ($command eq 'avgDistMat') {
	avgDistMat();
} else {
	die("Unknown command $command\n");
}

exit;

sub avgDistMat {

die(qq/
ngsMisc.pl avgDistMat [ID file] [square PHYLIP distance matrix]

ID file:
A 2-column TSV file where column 1 is individual name in the matrix (first column of the matrix file) 
and column 2 is the group name.e.g.:

CopChr2	Copadichromis_chrysonotus
CopChr5	Copadichromis_chrysonotus
CopChr8	Copadichromis_chrysonotus
MayZeb2	Maylandia_zebra
MayZeb3	-9
MayZeb5	Maylandia_zebra

The output matrix is the average distance between the groups in column 2. Individuals can be
excluded from the averaging by using '-9' as the group label.
\n/) if (!@ARGV || scalar(@ARGV) < 2);

# read in the IDs
open(my $idfh, '<', $ARGV[0]) or die("Unable to locate ID file $ARGV[0]: $!\n");
my (%id, %group, @order);
while (<$idfh>) {
	chomp;
	if ($_ =~ /^(\S+)\s+(\S+)/) {
		next if $2 eq '-9';
		$id{$1} = $2;
		if (!exists $group{$2}) {
			@{$group{$2}} = ();
			push @order, $2;
		}
	}
}
close $idfh;

# read in matrix
open(my $mfh, '<', $ARGV[1]) or die("Unable to locate distance matrix file $ARGV[1]: $!\n");
my @dist;
my $idx = 0;
while (<$mfh>) {
	@{$dist[$idx]} = split(/\s+/, $_);
	next if (scalar(@{$dist[$idx]}) < 2);
	my $label = shift @{$dist[$idx]};
	push @{$group{$id{$label}}}, $idx;
	$idx++;
}
close $mfh;

# average distances
my @avgdist;
for (my $i=0; $i <= $#order; $i++) {
	for (my $j=0; $j <= $#order; $j++) {
		$avgdist[$i]->[$j] = 0;
	}
}

for (my $i=0; $i <= $#order; $i++) {
	my $x = $order[$i];
	for (my $j=0; $j <= $#order; $j++) {
		my $y = $order[$j];
		if ($x eq $y) {
			next;
		} elsif ($avgdist[$j]->[$i] != 0) {
			$avgdist[$i]->[$j] = $avgdist[$j]->[$i];
		} else {
			foreach my $xidx (@{$group{$x}}) {
				foreach my $yidx (@{$group{$y}}) {
					$avgdist[$i]->[$j] += $dist[$xidx]->[$yidx];
				}
			}
			$avgdist[$i]->[$j] = sprintf "%.10f", $avgdist[$i]->[$j]/(scalar @{$group{$x}} * scalar @{$group{$y}});
		}
	}
}

# print averge distance
$"="\t";
print STDOUT scalar(@order), "\n";
my $k = 0;
foreach (@order) {
	print STDOUT "$_\t@{$avgdist[$k]}\n";
	$k++;
}

}
