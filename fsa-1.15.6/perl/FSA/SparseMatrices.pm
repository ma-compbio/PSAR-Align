#!/usr/bin/perl -w

=head1 NAME

FSA::SparseMatrices.pm

=head1 SYNOPSIS

Perl module encapsulating sparse alignment posterior probability matrices.

Written by Robert Bradley.

=head1 METHODS

=cut

package FSA::SparseMatrices;

use strict;
use Carp;

use vars '@ISA';

sub new {
  my ($class) = @_;

  my $self = {
	      'seq_keys' => {},    # sequence keys (names)
	      'seq_len' => {},     # sequence lengths
	      'match_probs' => {}, # seq1->seq2->pos1->pos2; transpose stored also
	      'gap_probs' => {},   # seq1->seq2->pos1
	      'sum_probs' => {},   # seq1->seq2->pos1
	      'max_probs' => {}    # seq1->seq2->pos1 is a 2-element hash (element 0 is max, 1 is altmax)
	     };

  bless $self, ref ($class) || $class;

  return $self;
}

# catch methods
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD; # $AUTOLOAD contains the fully qualified name of the original subroutine
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq 'DESTROY';

  # check for directives accessor, e.g. $self->directives_debug or $self->directives_('debug')
  if ($sub =~ /^directives_(\S*)$/i) {
    my $flag = lc ($1);
    $flag = shift @args unless (length $flag); # this catches the argument 'debug' in the second example usage given above
    if (!defined $self->{'directives'}->{$flag}) {
      $self->{'directives'}->{$flag} = "";     # if no such flag exists, create one with the value ""
    }                                          # we therefore have to test 'if ($self->directives_debug)' rather than 'if (defined $self->directives_debug)'
    return $self->{'directives'}->{$flag};     # the second will always be true because AUTOLOAD will create an empty flag for us
  }

  # check for ordinary accessors
  # This has the effect of automatically implementing getter and setter methods.
  # If there's an instance variable $name, then the getter method is $self->name
  # and the setter method is $self->name('newName')
  if (exists $self->{$sub}) {
    if (@args > 1) { croak "Usage: $sub() or $sub(newValue)"; }
    return
      @args                      # if @args > 0
      ? $self->{$sub} = $args[0] # { $self->{$sub} = $args[0]; return $args[0]; }
      : $self->{$sub};           # else { return $self->{$sub}; }
  }

  croak "Unsupported method: $sub\n";
}

# Load from a FSA .probs file.
sub from_probsfile {
  my ($class, $probsfile) = @_;

  my $self = $class->new;

  # parse post probs
  open PROBS,"<$probsfile" or croak "Couldn't open '$probsfile' for reading.\n";

  warn "Loading posterior alignment probabilities from '$probsfile'.\n";

  while (<PROBS>) {

    chomp;
    next if $_ eq '';

    # comment
    if (/^;/) { next; }

    # probability entry
    # (0, 4022) ~ (1, -1) => 0.773072
    # handle case of scientific notation
    elsif (/\((\d+),\s(-?\d+)\)\s~\s\((\d+),\s(-?\d+)\)\s=>\s(\d?\.?\d+)\*?[eE]?(-\d+)?/) {

      my ($seq1, $pos1) = ($1, $2);
      my ($seq2, $pos2) = ($3, $4);
      my $prob = $5;
      my $exponent = $6;

      if (defined $exponent) {
	$prob *= 10 ** $exponent;
      }

      # store
      $self->store_entry ($seq1, $seq2, $pos1, $pos2, $prob);

    }

    else {
      warn "Skipping input line '$_'.\n";
    }
    
  }

  close PROBS;

  return $self;
}

# Store a match or gap probability.
# Update all information about the sparse matrices accordingly.
sub store_entry {
  my ($self, $seq1, $seq2, $pos1, $pos2, $prob) = @_;

  # sequence keys
  $self->seq_keys->{$seq1} = 1;
  $self->seq_keys->{$seq2} = 1;

  # sequence lengths
  # remember 0-based indexing!
  if (!defined $self->seq_len->{$seq1}) { $self->seq_len->{$seq1} = $pos1; }
  else { $self->seq_len->{$seq1} = ($self->seq_len->{$seq1} < $pos1 + 1) ? $pos1 + 1 : $self->seq_len->{$seq1}; }
  if (!defined $self->seq_len->{$seq2}) { $self->seq_len->{$seq2} = $pos2; }
  else { $self->seq_len->{$seq2} = ($self->seq_len->{$seq2} < $pos2 + 1) ? $pos2 + 1 : $self->seq_len->{$seq2}; }

  # match probs
  # store transpose
  if (($pos1 >= 0) && ($pos2 >= 0)) {

    # match_probs
    $self->match_probs->{$seq1}->{$seq2}->{$pos1}->{$pos2} = $prob;
    $self->match_probs->{$seq2}->{$seq1}->{$pos2}->{$pos1} = $prob;

    # max_probs: seq1

    # if max_prob is undefined, then store new value
    if (!defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0]) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1} = [];
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] = $prob;
    }
    # if new max, then store new value and update alternate max
    elsif ($self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] < $prob) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] = $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0];
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] = $prob;
    }
    # if new alternate max because alternate max is undefined
    elsif (!defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1]) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] = $prob;
    }
    # if new alternate max
    elsif ($self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] < $prob) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] = $prob;
    }

    # max_probs: seq2

    # if max_prob is undefined, then store new value
    if (!defined $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0]) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2} = [];
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0] = $prob;
    }
    # if new max, then store new value and update alternate max
    elsif ($self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0] < $prob) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] = $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0];
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0] = $prob;
    }
    # if new alternate max because alternate max is undefined
    elsif (!defined $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1]) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] = $prob;
    }
    # if new alternate max
    elsif ($self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] < $prob) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] = $prob;
    }

    # sum prob: seq1
    if (!defined $self->sum_probs->{$seq1}->{$seq2}->{$pos1}) {
      $self->sum_probs->{$seq1}->{$seq2}->{$pos1} = $prob;
    } else {
      $self->sum_probs->{$seq1}->{$seq2}->{$pos1} += $prob;
    }

    # sum prob: seq2
    if (!defined $self->sum_probs->{$seq2}->{$seq1}->{$pos2}) {
      $self->sum_probs->{$seq2}->{$seq1}->{$pos2} = $prob;
    } else {
      $self->sum_probs->{$seq2}->{$seq1}->{$pos2} += $prob;
    }

  }

  # seq1 ~ -
  elsif (($pos1 >= 0) && ($pos2 == -1)) {

    # gap_probs
    $self->gap_probs->{$seq1}->{$seq2}->{$pos1} = $prob;

    # max_prob: seq1

    # if max_prob is undefined, then store new value
    if (!defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0]) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1} = [];
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] = $prob;
    }
    # if new max, then store new value and update alternate max
    elsif ($self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] < $prob) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] = $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0];
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] = $prob;
    }
    # if new alternate max because alternate max is undefined
    elsif (!defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1]) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] = $prob;
    }
    # if new alternate max
    elsif ($self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] < $prob) {
      $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1] = $prob;
    }

  }

  # - ~ seq2
  elsif (($pos1 == -1) && ($pos2 >= 0)) {

    # gap_probs
    $self->gap_probs->{$seq2}->{$seq1}->{$pos2} = $prob;

    # max_prob: seq2

    # if max_prob is undefined, then store new value
    if (!defined $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0]) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2} = [];
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0] = $prob;
    }
    # if new max, then store new value and update alternate max
    elsif ($self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0] < $prob) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] = $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0];
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[0] = $prob;
    }
    # if new alternate max because alternate max is undefined
    elsif (!defined $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1]) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] = $prob;
    }
    # if new alternate max
    elsif ($self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] < $prob) {
      $self->max_probs->{$seq2}->{$seq1}->{$pos2}->[1] = $prob;
    }

  } else {
    croak "Unreachable on input line:\n$_";
  }

}

# Write formatted .probs file to a string.
# In the format produced by FSA.
sub to_string {
  my ($self) = @_;

  my $s;
  my @seqs = sort { $a <=> $b } keys %{$self->seq_keys};
  for (my $i = 0; $i < @seqs; ++$i) {
    my $seq1 = $seqs[$i];

    for (my $j = $i + 1; $j < @seqs; ++$j) {
      my $seq2 = $seqs[$j];

      if (!$self->exists_matrix ($seq1, $seq2)) { next; }

      $s .= "
; Sparse posterior probability matrix for sequences $seq1 and $seq2
; Format is:
;   ($seq1, position_1) ~ ($seq2, position_2) => prob
; which means that ($seq1, position_1) is aligned to ($seq2, position_2) with probability prob.
;   ($seq1, position_1) ~ ($seq2, -1) => prob
; means that ($seq1, position_1) is aligned to a gap in $seq2 with probability prob.
; sequence is 0-based and position is 0-based
";

      $s .= "; match posteriors\n";
      foreach my $pos1 (sort { $a <=> $b } keys %{$self->match_probs->{$seq1}->{$seq2}}) {
	foreach my $pos2 (sort { $a <=> $b } keys %{$self->match_probs->{$seq1}->{$seq2}->{$pos1}}) {
	  $s .= "($seq1, $pos1) ~ ($seq2, $pos2) => " . $self->match_probs->{$seq1}->{$seq2}->{$pos1}->{$pos2} . "\n";
	}
      }
      $s .= "\n";

      $s .= "; gap posteriors\n";
      foreach my $pos1 (sort { $a <=> $b } keys %{$self->gap_probs->{$seq1}->{$seq2}}) {
	$s .= "($seq1, $pos1) ~ ($seq2, -1) => " . $self->gap_probs->{$seq1}->{$seq2}->{$pos1} . "\n";
      }
      $s .= "\n";
      foreach my $pos2 (sort { $a <=> $b } keys %{$self->gap_probs->{$seq2}->{$seq1}}) {
	$s .= "($seq1, -1) ~ ($seq2, $pos2) => " . $self->gap_probs->{$seq2}->{$seq1}->{$pos2} . "\n";
      }
      $s .= "\n";

    }
  }

  return $s;
}

# Do we have information for a particular sequence pair?
sub exists_matrix {
  my ($self, $seq1, $seq2) = @_;

  if (defined $self->match_probs->{$seq1} and
      defined $self->match_probs->{$seq1}->{$seq2}) {
    return 1;
  }

  return 0;
}

# Get match probability.
sub get_match_prob {
  my ($self, $seq1, $seq2, $pos1, $pos2) = @_;

  if (!$self->exists_matrix ($seq1, $seq2)) {
    croak "No probability matrix for sequences '$seq1' and '$seq2'.";
  }

  # if no entry, return 0
  if (!defined $self->match_probs->{$seq1}->{$seq2}->{$pos1} or
      !defined $self->match_probs->{$seq1}->{$seq2}->{$pos1}->{$pos2}) {
    return 0;
  }

  return assert_valid_prob ($self->match_probs->{$seq1}->{$seq2}->{$pos1}->{$pos2});
}

sub assert_valid_prob {
  my ($prob) = shift;

  if ($prob > 1.0001) { die "Probability of ", $prob, " exceeds bounds.\n"; }
  
  return ($prob > 1.0) ? 1.0 : $prob;
}

# Get gap probability.
sub get_gap_prob {
  my ($self, $seq1, $seq2, $pos1) = @_;

  if (!$self->exists_matrix ($seq1, $seq2)) {
    croak "No matrix for sequences '$seq1' and '$seq2';"
  }

  # there must exist an entry (FSA always outputs one)
  if (!defined $self->gap_probs->{$seq1}->{$seq2}->{$pos1}) {
    croak "No gap probability ($seq1, $pos1) ~ ($seq2, -1); I'm quitting.";
  }

  return $self->gap_probs->{$seq1}->{$seq2}->{$pos1};
}

# Get sum_pos2 P((seq1, pos1) ~ (seq2, pos2)).
# Cover the case of no entry for a match probability > 0.01.
# This coarse-graining of probabilities is important,
# since sum_prob is used as a denominator.
sub get_sum_prob {
  my ($self, $seq1, $seq2, $pos1) = @_;

  my $p = (defined $self->sum_probs->{$seq1}->{$seq2}->{$pos1}) ? $self->sum_probs->{$seq1}->{$seq2}->{$pos1} : 0.01;
  return assert_valid_prob ($p);
}

# Get max_pos2 P((seq1, pos1) ~ (seq2, pos2)),
# where pos2 can be a gap as well.
sub get_max_prob {
  my ($self, $seq1, $seq2, $pos1) = @_;

  # assert that max prob is defined
  # (because pos2 can be a gap, a max prob is always defined)
  if (!defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0]) {
    croak "A max prob is undefined; this should never happen!\n";
  }

  return assert_valid_prob ($self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0]);
}

# Get altmax_pos2 P((seq1, pos1) ~ (seq2, pos2)),
# where pos2 can be a gap as well.
# Note that we need to know $pos2 here (in contrast to get_max_prob),
# in order to determine whether the probability corresponding to
# the alignment ($seq1, $pos1) ~ ($seq2, $pos2)
# is the max, altmax or otherwise.
sub get_altmax_prob {
  my ($self, $seq1, $seq2, $pos1, $pos2) = @_;

  # assert that max prob is defined
  if (!defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0]) {
    croak "A max prob is undefined; this should never happen!\n";
  }
  # note that an altmax might not be defined,
  # ie if there is no match probability entry and so the max prob is
  # the gap prob

  my $p;
  if ($pos1 >= 0 && $pos2 >= 0) {
    $p = $self->get_match_prob ($seq1, $seq2, $pos1, $pos2);
  }
  elsif ($pos1 >= 0 && $pos2 == -1) {
    $p = $self->get_gap_prob ($seq1, $seq2, $pos1);
  }
  else {
    croak "No good; tried to get the altmax for a gap!\n";
  }

  # if $p is the max, then return the second-best alternative
  if (abs ($p - $self->get_max_prob ($seq1, $seq2, $pos1)) < 0.0001) {
    if (defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1]) {
      return $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1];
    } else {
      return 0.01;
    }
  }
  # else return the max
  else {
    return $self->get_max_prob ($seq1, $seq2, $pos1);
  }

}

# Compute the average certainty of the model.
# This corresponds to:
#  for each sequence seq1
#     for each sequence seq2
#        for each position pos1 in seq1
#           += P(pos1 ~ best)
#           += P(pos1 ~ next-best)
# Note that this calculation is only a measure of the width of the 
# probability distribution in the DP matrix, and as such is alignment-independent.
sub get_model_certainty {
  my ($self) = @_;

  my $max = 0;
  my $altmax = 0;

  foreach my $seq1 (keys %{$self->max_probs}) {
    for (my $pos1 = 0; $pos1 < $self->seq_len->{$seq1}; ++$pos1) {
      foreach my $seq2 (keys %{$self->max_probs->{$seq1}}) {
	# remember that altmax automatically gets initialized to minimum value
	# when max is initialized
	$max += $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0];
	if (defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1]) {
	  $altmax += $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1];
	} else {
	  $altmax += 0.01;
	}
      }
    }
  }

  return ($max / $altmax);
}

# As above, but show the entire distribution of slopes (for all characters).
sub show_model_certainty_slopes {
  my ($self, $precision) = @_;

  foreach my $seq1 (keys %{$self->max_probs}) {
    for (my $pos1 = 0; $pos1 < $self->seq_len->{$seq1}; ++$pos1) {
      foreach my $seq2 (keys %{$self->max_probs->{$seq1}}) {
	if (defined $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1]) {
	  printf "%.${precision}f\n", $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] / $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[1];
	} else {
	  printf "%.${precision}f\n", $self->max_probs->{$seq1}->{$seq2}->{$pos1}->[0] / 0.01;
	}
      }
    }
  }

}


1
