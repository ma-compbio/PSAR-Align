#!/usr/bin/perl -w

=head1 NAME

accuracy.pl

=head1 SYNOPSIS

Compute expected accuracy of an alignment under FSA's statistical model.

Written by Robert Bradley.

=cut

use strict;

use Stockholm;
use FSA::Model;

my $acc_threshold = 0.9;
my $cert_threshold = 5.0;
my $usage = "\nUsage: $0 -f <sequence file> [-p <.probs file>] [<MFA file(s)> <Stockholm file(s)>]

Calculates Acc, Sn, PPV, Certainty and Consistency for input alignments.

     [--annotate] annotates the input alignment and prints it as a Stockholm-format alignment
     [--accuracy-threshold <double>] threshold for calling a character correctly aligned (default is $acc_threshold)
     [--certainty-threshold <double>] threshold for calling a character certain (default is $cert_threshold)
     [--hard] use a hard cutoff for certainty calculation instead of a \"soft\" logarithmic transform
     [--model-certainty-dist] show distribution of \"slopes\" for model certainty; suppress all other output
     [--precision <double>] precision of output
     [--terse] terse output format

";

my @alignments;
my ($seqfile, $probsfile);
my $annotate = 0;
my $use_log = 1;
my $certainty_slopes = 0;
my $precision = 3;
my $terse = 0;

while (@ARGV) {
  my $arg = shift;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif (($arg eq "-f")) { $seqfile = shift; }
    elsif (($arg eq "-p")) { $probsfile = shift; }
    elsif (($arg eq "--annotate")) { $annotate = 1; }
    elsif (($arg eq "--accuracy-threshold")) { $acc_threshold = shift; }
    elsif (($arg eq "--certainty-threshold")) { $cert_threshold = shift; }
    elsif (($arg eq "--hard")) { $use_log = 0; }
    elsif (($arg eq "--model-certainty-dist")) { $certainty_slopes = 1; }
    elsif (($arg eq "--precision")) { $precision = shift; }
    elsif (($arg eq "--terse")) { $terse = 1; }
    else { die $arg, "\n", $usage; }
  }
  else { push @alignments, $arg; }
}

# test files ok
unless (defined $seqfile && -r $seqfile) { die "You must provide a readable FSA input sequence file (don't forget the -f switch!).\n"; }
if (!defined $probsfile) { $probsfile = $seqfile . ".probs"; }
if (! -r $probsfile) { die ".probs file '$probsfile' not readable."; }

# initialize FSA::Model object
my $model = FSA::Model->from_probsfile ($seqfile, $probsfile);

# loop over alignments
foreach my $align (@alignments) {

  my $stock = Stockholm->from_file ($align);
  $stock->assert_flush();
  $stock->drop_allgaps_columns();
  if (!$terse) {
    print "Looking at alignment '$align'...\n";
  }
  
  # if requested, just show the distribution of certainty slopes
  if ($certainty_slopes) {
    $model->sparse_matrices->show_model_certainty_slopes ($precision);
    exit 0;
  }

  # calculate accuracies
  $model->calc_accuracies ($stock, $acc_threshold, $cert_threshold, $use_log);

  # print the annotated alignment if requested
  if ($annotate) {
    $model->annotate_alignment ($stock, $precision);
    print $stock->to_string();
  }
  # else just show results
  $model->show_accuracies ($precision, $terse);

}
