#!/usr/bin/perl -w

=head1 NAME

prettify.pl

=head1 SYNOPSIS

Prettify (make human-readable) an alignment.

Written by Robert Bradley.

=cut

use Stockholm;

my $usage = "\nUsage: $0 <alignment file>

Prints a human-readable alignment in Stockholm format.

\n";

# make everything look nice
# drop all-gap columns

my $file = shift;
if (!defined $file) { die $usage; }

my $stk = Stockholm->from_file ($file);
$stk->assert_flush();

# drop all-gap columns
$stk->drop_allgaps_columns();

# display
print $stk->to_string();
