#!/usr/bin/perl -w

=head1 NAME

GFF.pm

=head1 SYNOPSIS

Perl module encapsulating a single GFF entry.
Note that the GFF format, which specifies that coordinates
are 1-based and fully-closed, is enforced by this module.

The original GFF module was written by Ian Holmes.
This version has been extended by Robert Bradley.

=head1 METHODS

=cut

package GFF;

use strict;
use Carp;

=head2 new

    my $gff = GFF->new()

Creates an empty GFF object.

=cut
sub new {
  my $class = shift;

  my $self = {
	      'seqid' => '.',
	      'source' => '.',
	      'type' => '.',
	      'start' => '.',
	      'end' => '.',
	      'score' => '.',
	      'strand' => '.',
	      'phase' => '.',
	      'attributes' => '.',
              'attributes_hash' => undef,
	      @_
              };

  bless $self, ref($class) || $class;

  return $self;
}

=head2 seqid

    my $seqid = $gff->seqid

Returns the seqid field.

    $gff->seqid ($newseqid)

Sets the seqid field.

Similarly implemented for all GFF fields.

Catches all methods by default.

=cut
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD; # $AUTOLOAD contains the fully qualified name of the original subroutine
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq 'DESTROY';

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

=head2 from_string

    my $gff = GFF->from_string ($str)

Reads values from a string in GFF format.

=cut 
sub from_string {
  my $class = shift;

  my $self = $class;
  if (!ref ($class)) { $self = $class->new(); }

  my $str = shift if (@_);
  if (!defined $str) { return undef; }

  chomp ($str);

  # check for:
  # empty lines
  # comment and sequence header lines
  # sequence lines
  if ($str eq ""
      || $str =~ /^\s*(#|>)/
      || $str =~ /^\s*\S+\s*$/)
    { return undef; }

  my @gff = split /\t/, $str;
  $self->{"seqid"} = $gff[0];
  $self->{"source"} = $gff[1];
  $self->{"type"} = $gff[2];
  $self->set_start ($gff[3]);
  $self->{"end"} = $gff[4];
  $self->{"score"} = $gff[5];
  $self->{"strand"} = $gff[6];
  $self->{"phase"} = $gff[7];
  $self->{"attributes"} = $gff[8];

  return $self;
}

=head2 to_string

    print $gff->to_string()

Output to a string in GFF format.

=cut 
sub to_string {
  my $self = shift;
  return join ("\t", $self->seqid, $self->source, $self->type, 
	       $self->start, $self->end, $self->score,
	       $self->strand, $self->phase, (defined $self->attributes ? $self->attributes : '') . "\n");
}

=head2 to_file

    $gff->to_file ($myfile)

Output to a file in GFF format.

    $gff->to_file ($myfile, 1)

Append mode.

=cut
sub to_file {
  my ($self, $pFile, $pMode) = @_;

  my $modeOp = ">";
  if ($pMode) {$modeOp = ">>";}
  open(FILE, $modeOp . $pFile) or die "Could not open $pFile.\n";
  print FILE $self->to_string();
}

=head2 set_start

   $gff->set_start ($start)

Set the start coordinate.
Enforces start coordinate >= 1.

=cut
sub set_start {
  my ($self, $s) = @_;
  if ($s < 1) { warn "Setting start coordinate '$s' to 1.\n"; $s = 1; }
  $self->start ($s);
}

=head2 add_value

   $gff->add_value()

Add single entry (key and value) to attributes field.

=cut 
sub add_value {
  my ($self, $key, $value) = @_;

  # add value to string for attributes field
  if (($self->attributes eq ".") || !length ($self->attributes)) {
    $self->attributes ("$key=$value");
  } elsif ($self->attributes =~ /;$/) {
    $self->{attributes} .= "$key=$value";
  } else {
    $self->{attributes} .= ";$key=$value";
  }

  # store in attributes hash if appropriate
  if ($self->{attributes_hash}) { ${$self->{attributes_hash}}{$key} = $value; }
}

=head2 get_attributes_hash 

   $gff->get_attributes_hash()

Get attributes field as a hash from keys to values.
Creates it if it doesn't already exist.

=cut 
sub get_attributes_hash {
  my $self = shift;

  my $hashRef = $self->{attributes_hash};
  my %hash;
  if (!$hashRef) {
    my @pairs = split (/;/, $self->{attributes});
    foreach my $pair (@pairs) {
      my @keyVal = split (/=/, $pair);
      $hash{$keyVal[0]} = $keyVal[1];
    }
    $hashRef = \%hash;
    # set attributes_hash
    $self->{attributes_hash} = $hashRef;
  }

  return $hashRef;

}

=head2 get_value

   $gff->get_value ($key)

Get value in attributes field.

=cut 
sub get_value {
  my ($self, $key) = @_;

  return ${$self->get_attributes_hash()}{$key};

}

=head2 set_value

   $gff->set_value ($key)

Set value in attributes field.

=cut 
sub set_value {
  my ($self, $key, $value) = @_;

  my $hashref = $self->get_attributes_hash();
  # add the value if it's new
  if (!defined $hashref->{$key}) {
    $self->add_value ($key, $value);
  }
  # else parse out and replace the old value
  # (yes, this is very slow, but I can't be bothered to fix it now)
  else {
    my $oldstr = "$key=" . $hashref->{$key};
    my $newstr = "$key=" . $value;
    $self->{attributes} =~ s/$oldstr/$newstr/;
  }

}

1
