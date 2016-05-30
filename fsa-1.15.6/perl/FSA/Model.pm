#!/usr/bin/perl -w

=head1 NAME

FSA::Model.pm

=head1 SYNOPSIS

Perl module encapsulating FSA's statistical model.

Written by Robert Bradley.

=head1 METHODS

=cut

package FSA::Model;

use strict;
use Carp;

use FSA::SparseMatrices;
use Stockholm;

use vars '@ISA';

sub new {
  my ($class) = @_;

  my $self = {
	      'sparse_matrices' => "", # SparseMatrices object
	      'seqname_map' => {},     # map sequence names to their SparseMatrices indices
	      'acc' => '',             # Acc (total accuracy)
	      'sn' => '',              # Sn (total sensivity)
	      'ppv' => '',             # PPV (total positive predictive value)
	      'cert' => '',            # Certainty (total certainty)
	      'cons' => '',            # Consistency (total consistency)
	      'acc_annot' => "",       # Acc per-column annotation
	      'sn_annot' => "",        # Sn per-column annotation
	      'ppv_annot' => "",       # PPV per-column annotation
	      'cert_annot' => "",      # Certainty per-column annotation
	      'cons_annot' => "",      # Consistency per-column annotation
	      'model_cert' => "",      # Model certainty (approximate average slope of probability distributions)
	      'exp_correct' => 0,      # expected number of correctly-aligned characters
	      'exp_incorrect' => 0,    # expected number of correctly-aligned characters
	      'total_chars' => 0       # total number of characters
	     };

  bless $self, ref ($class) || $class;

  return $self;
}

# Catch methods.
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
      : $self->{$sub};            # else { return $self->{$sub}; }
  }

  croak "Unsupported method: $sub\n";
}

# Initialization routines.
sub _initialize {
  my ($class, $seqfile) = @_;

  if (!defined $seqfile) { croak "You must define a sequence file!\n"; }

  my $self = $class->new;

  # create map from sequence names to indices
  unless (Stockholm->detect_FASTA ($seqfile)) {
    croak "'$seqfile' isn't a valid FASTA file.\n";
  }
  my $fasta = Stockholm->from_file ($seqfile);
  for (my $seq = 0; $seq < @{$fasta->seqname}; ++$seq) {
    $self->seqname_map->{@{$fasta->seqname}[$seq]} = $seq;
  }

  return $self;
}

# Load model from a FSA .probs file.
sub from_probsfile {
  my ($class, $seqfile, $probsfile) = @_;

  my $self = $class->_initialize ($seqfile);

  # initialize sparse matrices
  $self->sparse_matrices (FSA::SparseMatrices->from_probsfile ($probsfile));

  return $self;
}

# Initialize 0-1 posterior probabilities using the passed reference alignment.
sub from_alignment {
  my ($class, $seqfile, $alignfile) = @_;

  my $self = $class->_initialize ($seqfile);
  $self->sparse_matrices (FSA::SparseMatrices->new);

  # initialize sparse matrices from the passed alignment
  my $stock = Stockholm->from_file ($alignfile);
  # sanity checks
  $stock->assert_flush();
  $stock->drop_allgaps_columns();
  for (my $col = 0; $col < $stock->columns; ++$col) {

    foreach my $seq1 (@{$stock->seqname}) {
      # convert to index
      my $i = $self->seqname_map->{$seq1};
      my $ii = ($stock->is_gap ($seq1, $col)) ? -1 : ($stock->map_align_to_seq_coords ($col, $col, $seq1, 0))[0];
      
      foreach my $seq2 (@{$stock->seqname}) {
	my $j = $self->seqname_map->{$seq2};
	my $jj = ($stock->is_gap ($seq2, $col)) ? -1 : ($stock->map_align_to_seq_coords ($col, $col, $seq2, 0))[0];

	# skip gap-gap entries
	if (($ii < 0) && ($jj < 0)) { next; }

	# store with unit probability
	$self->sparse_matrices->store_entry ($i, $j, $ii, $jj, 1.0);
      
      }
    }

  }

  return $self;
}


# Calculate alignment accuracy information.
# Use either a "soft" logarithmic transform
# or a hard threshold for the certainty calculation.
sub calc_accuracies {
  my ($self, $stock, $acc_threshold, $cert_threshold, $use_log) = @_;

  if (!defined $acc_threshold) { $acc_threshold = 0.9; }
  if (!defined $cert_threshold) { $cert_threshold = 5.0; }
  if (!defined $use_log) { $use_log = 1; }

  # initialize temporary variables to 0
  # (these will be used to populate the alignment-specific variables in $self)

  # accuracy (Acc)
  my $acc_num = 0;
  my $acc_denom = 0;
  my $acc_annot = "";

  # positive predictive value (PPV)
  my $ppv_num = 0;
  my $ppv_denom = 0;
  my $ppv_annot = "";  

  # sensitivity (Sn)
  my $sn_num = 0;
  my $sn_denom = 0;
  my $sn_annot = "";

  # certainty
  my $cert_num = 0;
  my $cert_denom = 0;
  my $cert_annot = "";

  # consistency
  my $cons_num = 0;
  my $cons_denom = 0;
  my $cons_annot = "";

  # expected numbers of correctly and incorrectly-aligned characters
  # (meeting our $acc_threshold)
  my $exp_correct = 0;
  my $exp_incorrect = 0;
  my $total_chars = 0;

  # first calculate Acc and PPV
  # loop over columns
  for (my $col = 0; $col < $stock->columns; ++$col) {

    # per-column values
    my $col_acc_num = 0;
    my $col_acc_denom = 0;
    
    my $col_ppv_num = 0;
    my $col_ppv_denom = 0;

    my $col_sn_num = 0;
    my $col_sn_denom = 0;
    
    my $col_cert_num = 0;
    my $col_cert_denom = 0;

    my $col_cons_num = 0;
    my $col_cons_denom = 0;

    # loop over seqs in column
    foreach my $seq1 (@{$stock->seqname}) {

      # convert to index
      my $i = $self->seqname_map->{$seq1};

      # estimate whether it's correctly-aligned
      my $char_acc_num = 0;
      my $char_acc_denom = 0;

      # loop over all seqs
      foreach my $seq2 (@{$stock->seqname}) {

	my $j = $self->seqname_map->{$seq2};

	# skip the diagonal
	next if ($i == $j);

	# check that matrix is available
	if (!$self->sparse_matrices->exists_matrix ($i, $j)) { next; }
	
	# skip gap-gap entries
	if ($stock->is_gap ($seq1, $col) && $stock->is_gap ($seq2, $col)) { next; }

	# if seq1 is gapped
	if ($stock->is_gap ($seq1, $col)) {

	  # get position in seq j
	  my $jj = ($stock->map_align_to_seq_coords ($col, $col, $seq2, 0))[0]; # final argument to map_align_to_seq_coords
                                                                                # specifies 0-based coordinates         
	  my $p = $self->sparse_matrices->get_gap_prob ($j, $i, $jj);

	  $col_acc_num += $p;
	  $col_acc_denom += 1;

	}

	# if seq2 is gapped
	elsif ($stock->is_gap ($seq2, $col)) {

	  # get position in seq i
	  my $ii = ($stock->map_align_to_seq_coords ($col, $col, $seq1, 0))[0];

	  my $p = $self->sparse_matrices->get_gap_prob ($i, $j, $ii);

	  $col_acc_num += $p;
	  $col_acc_denom += 1;

	  $char_acc_num +=  $p;
	  $char_acc_denom += 1;
	  
	  # numerator for sn doesn't increase
	  # but denominator does
	  $col_sn_denom += $self->sparse_matrices->get_sum_prob ($i, $j, $ii);

	  # both numerator and denominator for consistency increase
	  $col_cons_num += $p;
	  $col_cons_denom += $self->sparse_matrices->get_max_prob ($i, $j, $ii);

	  if ($use_log) {
	    $col_cert_num += $p;
	    $col_cert_denom += $self->sparse_matrices->get_altmax_prob ($i, $j, $ii, -1);
	  } else {
	    $col_cert_num += ($p / $self->sparse_matrices->get_altmax_prob ($i, $j, $ii, -1) >= $cert_threshold) ? 1 : 0;
	    $col_cert_denom += 1;
	  }

	}

	# if match
	else {

	  # get position in seqs i and j
	  my $ii = ($stock->map_align_to_seq_coords ($col, $col, $seq1, 0))[0];
	  my $jj = ($stock->map_align_to_seq_coords ($col, $col, $seq2, 0))[0];
	    
	  my $p = $self->sparse_matrices->get_match_prob ($i, $j, $ii, $jj);

	  $col_acc_num += 2 * $p;
	  $col_acc_denom += 2;

	  $char_acc_num += 2 * $p;
	  $char_acc_denom += 2;

	  $col_ppv_num += 2 * $p;
	  $col_ppv_denom += 2;
	  
	  # both numerator and denominator for sn increase
	  $col_sn_num += $p;
	  $col_sn_denom += $self->sparse_matrices->get_sum_prob ($i, $j, $ii);

	  # both numerator and denominator for consistency increase
	  $col_cons_num += $p;
	  $col_cons_denom += $self->sparse_matrices->get_max_prob ($i, $j, $ii);

	  if ($use_log) {
	    $col_cert_num += 2 * $p;
	    $col_cert_denom += 2 * $self->sparse_matrices->get_altmax_prob ($i, $j, $ii, $jj);
	  } else {
	    $col_cert_num += 2 * (($p / $self->sparse_matrices->get_altmax_prob ($i, $j, $ii, $jj) >= $cert_threshold) ? 1 : 0);
	    $col_cert_denom += 2 * 1;
	  }

	}

      }

      # is it correctly aligned according to our threshold?
      if (!$stock->is_gap ($seq1, $col)) {
	++$total_chars;
	if ($char_acc_num / $char_acc_denom > $acc_threshold) { ++$exp_correct; }
	elsif ($char_acc_num / $char_acc_denom < 1.0 - $acc_threshold) { ++$exp_incorrect; }
      }

    }

    # increment alignment counters
    $acc_num += $col_acc_num;
    $acc_denom += $col_acc_denom;

    $ppv_num += $col_ppv_num;
    $ppv_denom += $col_ppv_denom;

    $sn_num += $col_sn_num;
    $sn_denom += $col_sn_denom;

    $cons_num += $col_cons_num;
    $cons_denom += $col_cons_denom;

    $cert_num += $col_cert_num;
    $cert_denom += $col_cert_denom;

    # sanity checks
    if (($col_acc_num / $col_acc_denom < 0) || ($col_acc_num / $col_acc_denom > 1)) {
      warn "Dangerous-looking Acc for column $col: $col_acc_num / $col_acc_denom.\n";
    }
    if (($col_ppv_denom > 0) and
	($col_ppv_num / $col_ppv_denom < 0) || ($col_ppv_num / $col_ppv_denom > 1)) {
      warn "Dangerous-looking PPV for column $col: $col_ppv_num / $col_ppv_denom.\n";
    }

    # modify per-column annotation as appropriate
    substr ($acc_annot, $col, 1) = get_markup_char ($col_acc_num / $col_acc_denom);
    substr ($ppv_annot, $col, 1) = ($col_ppv_denom > 0) ? get_markup_char ($col_ppv_num /$col_ppv_denom) : "-";
    substr ($sn_annot, $col, 1) = get_markup_char ($col_sn_num / $col_sn_denom);
    substr ($cons_annot, $col, 1) = get_markup_char ($col_cons_num / $col_cons_denom);
    if ($use_log) {
      substr ($cert_annot, $col, 1) = get_markup_char ($self->normalize_certainty ($col_cert_num / $col_cert_denom, $cert_threshold));
    } else {
      substr ($cert_annot, $col, 1) = get_markup_char ($col_cert_num / $col_cert_denom);
    }

  }

  # store total values
  $self->acc ($acc_num / $acc_denom);
  $self->sn ($sn_num / $sn_denom);
  $self->ppv ($ppv_num / $ppv_denom);
  $self->cons ($cons_num / $cons_denom);
  if ($use_log) {
    $self->cert ($self->normalize_certainty ($cert_num / $cert_denom, $cert_threshold));
  } else {
    $self->cert ($cert_num / $cert_denom);
  }
  $self->model_cert ($self->sparse_matrices->get_model_certainty());

  $self->exp_correct ($exp_correct);
  $self->exp_incorrect ($exp_incorrect);
  $self->total_chars ($total_chars);

  # store per-column annotation lines
  $self->acc_annot ($acc_annot);
  $self->sn_annot ($sn_annot);
  $self->ppv_annot ($ppv_annot);
  $self->cons_annot ($cons_annot);
  $self->cert_annot ($cert_annot);

}

# Normalize certainty to [0,1] with a logarithmic transform.
# Use the logarithm transform on values below $cert_threshold.
sub normalize_certainty {
  my ($class, $c, $cert_threshold) = @_;

  if (!defined $cert_threshold) { $cert_threshold = 5.0; }

  # if there was a better alternative, then return 0
  if ($c < 1.0) {
    return 0;
  }

  # if it was $cert_threshold times better than the best alternative, then return 1
  elsif ($c >= $cert_threshold) {
    return 1;
  }

  # else perform a logarithmic transform
  else {
    return (log ($c) / log ($cert_threshold));
  }

}

# Show alignment accuracy information.
sub show_accuracies {
  my ($self, $precision, $terse) = @_;

  if (!defined $precision) { $precision = 3; }
  if (!defined $terse) { $terse = 0; }

  # summary of accuracy metrics
  my $acc_summary = sprintf "%.${precision}f %.${precision}f %.${precision}f %.${precision}f %.${precision}f", $self->acc, $self->sn, $self->ppv, $self->cert, $self->cons;
  # model certainty
  my $model_cert = sprintf "%.${precision}f", $self->model_cert;
  # estimated number of definitely correct and incorrect characters`
  my $correctness = sprintf "%d %d %d", $self->exp_correct, $self->exp_incorrect, $self->total_chars;

  if (!$terse) {
    printf "Acc %.${precision}f\n", $self->acc;
    printf "Sn  %.${precision}f\n", $self->sn;
    printf "PPV %.${precision}f\n", $self->ppv;
    printf "Certainty %.${precision}f\n", $self->cert;
    printf "Consistency %.${precision}f\n", $self->cons;
    printf "Model certainty (slope) %.${precision}f\n", $self->model_cert;
    printf "Definitely correct characters %d / %d\n", $self->exp_correct, $self->total_chars;
    printf "Definitely incorrect characters %d / %d\n", $self->exp_incorrect, $self->total_chars;
    printf "%-10s %s %s %s\n", "Summary:", $acc_summary, $model_cert, $correctness;
  }
  else {
    print "$acc_summary $model_cert $correctness\n";
  }


}


# Annotate Stockholm alignment.
sub annotate_alignment {
  my ($self, $stock, $precision) = @_;
  
  # #=GF annotations
  $stock->add_gf ("Acc", sprintf "%.${precision}f", $self->acc);
  $stock->add_gf ("Sn", sprintf "%.${precision}f", $self->sn);
  $stock->add_gf ("PPV", sprintf "%.${precision}f", $self->ppv);
  $stock->add_gf ("Consistency", sprintf "%.${precision}f", $self->cons);

  # #=GC annotations
  $stock->gc->{"Acc"} = $self->acc_annot;
  $stock->gc->{"Sn"} = $self->sn_annot;
  $stock->gc->{"PPV"} = $self->ppv_annot;
  $stock->gc->{"Consistency"} = $self->cons_annot;

  if (!$stock->is_flush()) {
    croak "Alignment not flush.\n";
  }

}


# Convert a number in the interval [0-1] to an integer [0-9].
sub get_markup_char {
  my ($norm) = @_;

  my $safe = int (10 * $norm);
  $safe = ($safe > 9) ? 9 : $safe;
  my $char = substr ("0123456789", $safe, 1);

  return $char;
}


1
