#!/usr/bin/perl -w

=head1 NAME

GFF::Database.pm

=head1 SYNOPSIS

Perl module encapsulating a GFF file ("database").

The original GFF::Database module was written by Ian Holmes.
This version has been extended by Robert Bradley.

=cut

package GFF::Database;

use strict;
use Carp;

use GFF;

=head2 new

    my $gffdb = GFF::Database->new()

Creates an empty GFF::Database object.

=cut
sub new {
  my ($class) = @_;

  my $self = {
	      'db' => []
	     };

  bless $self, ref ($class) || $class;

  return $self;
}

=head2 db

    my $records = $gffdb->db

Returns the list of GFF records contained.

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

=head2 from_file

    my $gffdb = GFF::Database->from_file ($filename)

Reads GFF file.

=cut 
sub from_file {
  my ($class, $filename) = @_;

  my $self = $class->new;

  $self->append_from_file ($filename);

  return $self;
}

=head2 append_from_file

    my $gffdb->append_from_file ($filename)

Appends from GFF file.

=cut 
sub append_from_file {
  my ($self, $filename) = @_;

  my $gff;
  local *FILE;
  local $_;
  open FILE, "<$filename" or croak "Couldn't open $filename: $!";
  while (<FILE>) {

    # skip empty lines
    chomp;
    if ($_ eq "") { next; }

    # skip comment and sequence header lines
    if (/^\s*(#|>)/) { next; }

    # skip sequence lines
    if (/^\s*\S+\s*$/) { next; }
      
    $gff = GFF->new unless defined $gff;
    if ($gff->from_string ($_)) {
      $self->add ($gff);
      $gff = undef;
    }

  }
  close FILE;

  # catch end case
  $self->add ($gff) if defined $gff;

}

=head2 add

    $gffdb->add ($gff)

Add a single GFF entry to database.

=cut
sub add {
  my ($self, $gff, $verbose) = @_;
  push @{$self->db}, $gff;
  warn "...loaded GFF record ", $gff->to_string() if $verbose;
  return $self;
}

=head2 to_string

    print $gff->to_string()

Output to a string in GFF format.

=cut 
sub to_string {
  my ($self) = @_;

  my $str = join ("", map ($_->to_string(), @{$self->db}));

  return $str;

}

=head2 find_by_group

    $gff->find_by_group ($key, $value)

Query for GFF records matching a given key/value in the group field.

=cut 
sub find_by_group {
  my ($self, $pKey, $pVal) = @_;

  my @gffList;
  foreach my $gff (@{$self->{db}})
  {
    my $group_hash = $gff->get_group_hash();
    if (${$group_hash}{$pKey} eq $pVal) { push(@gffList, $gff); }
  }

  return \@gffList;
}

=head2 find_by_type

    my $list = $gff->find_by_type ("type")

Get all entries with a particular type field.

=cut
sub find_by_type {
  my ($self, $type) = @_;

  my @list;
  foreach my $gff (@{$self->{db}}) {
    if ($gff->type eq $type) {
      push @list, $gff;
    }
  }

  return \@list;
}

1
