#!/usr/bin/perl -w

=head1 NAME

cmpalign.pl

=head1 SYNOPSIS

Compare two alignments.

Written by Robert Bradley.

=cut

use strict;
use Stockholm;
use Stockholm::Database;

my $usage = "\nUsage: $0 <reference alignment> <comparison alignment>

Compares two alignments.
Aligments can be in Stockholm, multi-FASTA, ClustalW or MSF formats.

     [--precision <double>] precision of output
     [--terse] terse output format
     [--log] log progress
     [--lazy] allow sequence data in two alignments to differ
     [--rna] allow U and T to be interchangeable
     [--pairwise-reference] reference alignment is a Stockholm database of pairwise alignments
                             (for e.g. SABmark, which provides pairwise reference alignments)

\n";

my $lazy = 0;
my $rna = 0;
my $log = 0;
my ($reffile, $cmpfile);
my $precision = 3;
my $terse = 0;
my $pairwiseref = 0;

my @argv;
while (@ARGV) {
  my $arg = shift;
  if ($arg =~ /^-/) {
    if (($arg eq "-h") || ($arg eq "--help")) { print $usage; exit; }
    elsif ($arg eq "--precision") { $precision = shift; }
    elsif ($arg eq "--terse") { $terse = 1; }
    elsif ($arg eq "--log") { $log = 1; }
    elsif ($arg eq "--lazy") { $lazy = 1; }
    elsif ($arg eq "--rna") { $rna = 1; }
    elsif ($arg eq "--pairwise-reference") { $pairwiseref = 1; }
    else { die $usage; }
  } else {
    push @argv, $arg;
  }
}
if (@argv == 2) {
  ($reffile, $cmpfile) = @argv;
}
else {
  die $usage;
}

my $gapChars = '-._';

my $refdb;  # fake a Stockholm::Database
if ($pairwiseref) { $refdb = Stockholm::Database->from_file ($reffile); }
else { push @{$refdb}, Stockholm->from_file ($reffile); }
my $cmp = Stockholm->from_file ($cmpfile);
$cmp->assert_flush();
$cmp->drop_allgaps_columns();

# make sure that the 2 alignments contain the same sequence names and data
# convert seqs to uppercase as side effect
my %seqname;
foreach my $ref (@{$refdb}) {
  foreach my $seq (@{$ref->seqname}) {
    $seqname{$seq} = 1;

    $ref->assert_flush();
    $ref->drop_allgaps_columns();

    unless (defined $cmp->seqdata->{$seq}) { die "Sequence data for '$seq' not in comparison alignment; I'm bailing.\n"; }
    # make sure that they're the same sequence
    if (ungapped ($ref->seqdata->{$seq}, $rna) ne ungapped ($cmp->seqdata->{$seq}, $rna)) {
      if ($lazy) { warn "Sequence data for '$seq' doesn't match:\n",ungapped ($ref->seqdata->{$seq}, $rna),"\n",ungapped ($cmp->seqdata->{$seq}, $rna),"\n  (but I'm doing the comparison regardless)\n"; }
      else { die "Sequence data for '$seq' doesn't match:\n",ungapped ($ref->seqdata->{$seq}, $rna),"\n",ungapped ($cmp->seqdata->{$seq}, $rna),"\n"; }
    }

    # convert everything to uppercase for safety
    $ref->seqdata->{$seq} = uc ($ref->seqdata->{$seq});
    $cmp->seqdata->{$seq} = uc ($cmp->seqdata->{$seq});

  }
}

# keep seqnames sorted (otherwise nasty things happen!)
my @seqname_sorted = sort keys %{seqname};

# total number of characters (sum-of-pairs)
my $total_char = 0;
# total number of residue pairs (sum-of-pairs) in each alignment
my ($sn1, $sn2) = (0, 0); # (ref, cmp)

# (num_aligned_pairs, num_unaligned_1, num_unaligned_2)
my @overlap = (0,0,0);

# (num_aligned_pairs, num_unaligned_1, num_unaligned_2)
my @stats_ref = (0,0,0);
my @stats_cmp = (0,0,0);

# for each pair of sequences

# for logging
my $num_seqs = scalar @seqname_sorted;
my $num_seq_pairs = $num_seqs * ($num_seqs - 1) / 2;
my $cnt = 0;
for (my $s1 = 0; $s1 < @seqname_sorted; $s1++) {
  for (my $s2 = $s1+1; $s2 < @seqname_sorted; $s2++) {
    my ($seqname1, $seqname2) = ($seqname_sorted[$s1], $seqname_sorted[$s2]);

    # get the reference alignment containing the sequences
    # (necessary for the case where the reference alignment is a Stockholm::Database
    # of pairwise comparisons)
    my $ref = get_alignment_with ($refdb, $seqname1, $seqname2);

    # increment total number of characters
    $total_char += length (ungapped ($ref->seqdata->{$seqname1}, $rna)) + length (ungapped ($ref->seqdata->{$seqname2}, $rna));

    # for the reference alignment
    # $align_cmp->{$i} = $j means that the positions ($i,$j) are aligned
    my ($row0, $row1) = ($ref->seqdata->{$seqname1}, $ref->seqdata->{$seqname2});
    my $align_ref = aligned_pairwise ($row0, $row1);
    my @indel_ref = @{indel_single ($row0,$row1)};
    $stats_ref[0] += scalar keys %$align_ref;
    $stats_ref[1] += scalar keys %{$indel_ref[0]};
    $stats_ref[2] += scalar keys %{$indel_ref[1]};

    # for the other alignment
    ($row0, $row1) = ($cmp->seqdata->{$seqname1}, $cmp->seqdata->{$seqname2});
    my $align_cmp = aligned_pairwise ($row0, $row1);
    my @indel_cmp = @{indel_single ($row0,$row1)};
    $stats_cmp[0] += scalar keys %$align_cmp;
    $stats_cmp[1] += scalar keys %{$indel_cmp[0]};
    $stats_cmp[2] += scalar keys %{$indel_cmp[1]};

    # compare aligned characters
    foreach my $i (sort { $a <=> $b } keys %{$align_ref}) {
      if (defined (my $j = $align_ref->{$i})) {
	$sn1++; # reference alignment
	if (defined $align_cmp->{$i} && $j == $align_cmp->{$i}) { # b/c $i isn't necessarily aligned in cmp
	  ++$overlap[0];
	}
      }
    }

    # compare unaligned (gapped) characters
    foreach my $i (sort { $a <=> $b } keys %{$indel_ref[0]}) {
      if (defined $indel_cmp[0]->{$i}) { ++$overlap[1]; }
    }
    foreach my $i (sort { $a <=> $b } keys %{$indel_ref[1]}) {
      if (defined $indel_cmp[1]->{$i}) { ++$overlap[2]; }
    }

    # now get counts for comparison alignment (for e.g. specificity measure)
    foreach my $i (sort { $a <=> $b } keys %{$align_cmp}) {
      if (defined (my $j = $align_cmp->{$i})) {
	$sn2++;
      }
    }

    # log progress if desired
    if (($cnt++ % 100 == 0) && $log) {
      my $percent_done = int (100 * $cnt / $num_seq_pairs + 0.5);
      warn "Processed sequence pair '$seqname1' and '$seqname2'; ", $percent_done, "% (", $cnt, "/", $num_seq_pairs, ") complete.\n";
    }

  }
}

# calculate and print the results
my $acc  = (2*$overlap[0] + $overlap[1]+ $overlap[2]) / $total_char;
my $sn = ($sn1 > 0) ? $overlap[0] / $sn1 : 0;
my $ppv = ($sn2 > 0) ? $overlap[0] / $sn2 : 0;

my $summary = $stats_ref[0] . " " . $stats_ref[1] . " " . $stats_ref[2]
  . " " . $stats_cmp[0] . " " . $stats_cmp[1] . " " . $stats_cmp[2]
  . " " . $overlap[0] . " " . $overlap[1] . " " . $overlap[2];

# terse results only if requested
if ($terse) {
  print "$summary\n";
}

# else pretty results
else {
  print "Alignment '$reffile' has ", $stats_ref[0], " aligned character pairs, ", $stats_ref[1], " insertions, ", $stats_ref[2], " deletions.\n";
  print "Alignment '$cmpfile' has ", $stats_cmp[0], " aligned character pairs, ", $stats_cmp[1], " insertions, ", $stats_cmp[2], " deletions.\n";
  print "Overlap is ", $overlap[0], " aligned character pairs, ", $overlap[1], " insertions, ", $overlap[2], " deletions.\n";
  printf "Acc %.${precision}f\n", $acc;
  printf "Sn  %.${precision}f\n", $sn;
  printf "PPV %.${precision}f\n", $ppv;
  print "Summary: $summary\n";
}

# Remove gaps from a string.
# Convert to uppercase for safety.
# Converts T to U if requested.
sub ungapped {
  my ($row, $rna) = @_;
  if (!defined $rna) { $rna = 0; }

  my $new = uc ($row);
  if ($rna) { $new =~ tr/T/U/; }
  $new =~ s/[$gapChars]//g;

  return $new;
}

# Takes 2 strings (alignments rows) as input.
# Returns a hash reference where $align_cmp->{$i} = $j means that the positions ($i,$j) are aligned.
sub aligned_pairwise {
  my ($row0, $row1) = @_;

  my $align_cmp = {};
  my $cols = length ($row0);
  my ($i, $j) = (0, 0); # keep track of positions in seq0, seq1
  for (my $col = 0; $col < $cols; $col++) {
    my $c0 = substr ($row0, $col, 1);
    my $c1 = substr ($row1, $col, 1);
    if (!($c0 =~ /[$gapChars]/)) {
      if (!($c1 =~ /[$gapChars]/)) {
	$align_cmp->{$i} = $j;
      }
      $i++;
    }
    if (!($c1 =~ /[$gapChars]/)) {
      $j++;
    }
  }

  return $align_cmp;
}

# Takes 2 strings (alignment rows) as input.
# Returns a hash reference where $indel->[0]->{$i} = 1 means that position $i (alignment coords) is gapped in sequence 0.
sub indel_single {
  my ($row0, $row1) = @_;

  my $indel = [];
  my $cols = length ($row0);
  my ($i, $j) = (0, 0); # keep track of positions in seq0, seq1
  for (my $col = 0; $col < $cols; $col++) {
    my $c0 = substr ($row0, $col, 1);
    my $c1 = substr ($row1, $col, 1);
    if (!($c0 =~ /[$gapChars]/) && ($c1 =~ /[$gapChars]/)) { $indel->[0]->{$i} = 1; ++$i; } # ungapped in 0 (i)
    elsif (($c0 =~ /[$gapChars]/) && !($c1 =~ /[$gapChars]/)) { $indel->[1]->{$j} = 1; ++$j; } # ungapped in 1 (j)
    elsif (!($c0 =~ /[$gapChars]/) && !($c1 =~ /[$gapChars]/)) { ++$i; ++$j; }
  }

  return $indel;
}

# Find the alignment containing the passed sequences.
sub get_alignment_with {
  my ($db, $seq1, $seq2) = @_;

  foreach my $stk (@$db) {
    if (exists $stk->seqdata->{$seq1} && exists $stk->seqdata->{$seq2}) {
      return $stk;
    }
  }

  die "Couldn't find a reference alignment containing '$seq1' and '$seq2'.";

}
