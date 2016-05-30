#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use Cwd 'abs_path';

if ($#ARGV != 2) {
	print STDERR "Usage: ./PSARAlign.pl <alignment file> <PSAR parameter file> <output directory>\n";
	print STDERR "Ex: ./PSARAlign.pl examples/test.fa examples/parameters.txt myout\n";
	exit;
}

my $aln_f = shift;
my $par_f = shift;
my $out_dir = shift;

system("mkdir -p $out_dir");
$aln_f = abs_path($aln_f);
$par_f = abs_path($par_f);
$out_dir = abs_path($out_dir);

chdir($out_dir);

# Run PSAR
`$Bin/psar $aln_f $par_f samples`;

# Check for the samples
my $numsamples = `ls samples/ | wc -l`;
chomp($numsamples);
if ($numsamples == 0) {
	`cp $aln_f revised.aln.fa`;
	exit;	
}

# Compute posterior probabilities of aligning two residues
`$Bin/compute_probs.pl samples > post_probs.txt`;

# Convert an alignment to sequences
open(O,">seq.fa");
open(F,"$aln_f");
my ($spc, $seq) = ("", "");
while(<F>) { 
	chomp;
	if (length($_) == 0 || $_ =~ /^#/) { next; }

	if ($_ =~ /^>(\S+)/) {
		if (length($spc) > 0) {
			$seq =~ s/\-//;
			print O ">$spc\n";
			print O "$seq\n";
		}
		
		$spc = $1;
		$seq = "";
	} else {
		$seq .= "$_";
	}
}
if (length($spc) > 0) {
	$seq =~ s/\-//;
	print O ">$spc\n";
	print O "$seq\n";
}
close(F);
close(O);

# Alignment by FSA
`sed -e 's/\-//g' $aln_f > seq.fa`;
`$Bin/fsa-1.15.6/src/main/fsa --noanchored --load-probs post_probs.txt seq.fa > revised.aln.fa`; 
