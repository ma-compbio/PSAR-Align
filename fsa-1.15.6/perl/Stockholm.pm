#!/usr/bin/perl -w

=head1 NAME

Stockholm.pm

=head1 SYNOPSIS

Perl module encapsulating a Stockholm multiple alignment.
Can load Stockholm, multi-FASTA (MFA), CLUSTAL and MSF alignments.
Can print Stockholm and MFA aligments.

The original Stockholm module was written by Ian Holmes for his DART library.
This version has been modified and extended by Robert Bradley.

=head1 METHODS

=cut

package Stockholm;

use strict;
use vars '@ISA';

use Carp qw(carp cluck croak);

# define the gap alphabet
my $gapChars = '-._';

# define percent id annotation
my $percentid_annot = 'Percent_id';

# define gap fraction annotation
my $gapfraction_annot = 'Gap_fraction';


=head2 new

    my $stock = Stockholm->new();

Creates an empty Stockholm object.

=cut
sub new {
  my ($class) = @_;

  my $self = {
	      'seqname' => [],  # order of sequences
	      'seqdata' => {},  # sequence data itself
	      'gf' => {},       # Generic File annotation (freeform)
	      'gc' => {},       # Generic Consensus annotation (by-column)
	      'gs' => {},       # Generic Sequence annotation (by-sequence, freeform)
	      'gr' => {},       # Generic by-Row annotation (by-sequence, by-column)
	      'gforder' => []   # order of GF lines
	     };

  bless $self, ref ($class) || $class;

  return $self;
}

=head2 from_file

    my $stock = Stockholm->from_file ($filename);
    my $stock = Stockholm->from_file ($filename, 1);

Creates a new Stockholm object and reads it from a file.
Attempts to automatically detect the format of the input file
and parse accordingly.

If last argument is true, then call strip_leading_chr to remove
the 'chr' frequently prepended to chromosome names.

=cut
sub from_file {
  my ($class, $filename, $remove_chr) = @_;

  my $self;

  # by default don't remove the leading 'chr'
  if (!defined $remove_chr) { $remove_chr = 0; }

  # FASTA
  if ($class->detect_FASTA ($filename)) {
    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    $self = $class->from_filehandle_FASTA (\*FILE, $remove_chr);
    close FILE;
  }

  # CLUSTAL
  elsif ($class->detect_CLUSTAL ($filename)) {
    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    $self = $class->from_filehandle_CLUSTAL (\*FILE, $remove_chr);
    close FILE;
  }

  # MSF
  elsif ($class->detect_MSF ($filename)) {
    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    $self = $class->from_filehandle_MSF (\*FILE, $remove_chr);
    close FILE;
  }

  # Stockholm
  else {
    local *FILE;
    open FILE, "<$filename" or croak "Couldn't open $filename: $!";
    $self = $class->from_filehandle_Stockholm (\*FILE, $remove_chr);
    close FILE;
  }

  return $self;

}

=head2 detect_FASTA

    detect_FASTA ($filename);

Detects whether the alignment is FASTA format.
Looks for lines of the form
>seqname

=cut
sub detect_FASTA {
  my ($class, $filename) = @_;

  local *FILE;
  open FILE, "<$filename" or croak "Couldn't open '$filename'.\n";
  while (<FILE>) {
    if (/^\s*>\s*\S+/) { return 1; }
  }
  close FILE;

  return 0;
}

=head2 from_filehandle_FASTA

    $self->from_filehandle_FASTA ($filehandle);

Creates a new Stockholm object and reads it from a filehandle.

=cut
sub from_filehandle_FASTA {
  my ($class, $filehandle, $remove_chr) = @_;
  
  my $self = $class->new;

  # by default don't remove the leading 'chr'
  if (!defined $remove_chr) { $remove_chr = 0; }

  my $name = "";
  my $data = "";

  # perform a first pass through the data to put each alignment row
  # onto a single line
  # this *dramatically* speeds up reading in FASTA alignments
  local $_;
  my $str = "";
  while (<$filehandle>) {
    if (/^\s*>/) { $str .= "\n$_"; next; }
    chomp;
    $str .= $_;
  }

  # now read in this compressed format
  my $str_handle;
  open ($str_handle, "<",  \$str);
  while (<$str_handle>) {
    chomp;

    # if beginning of new sequence
    if (/^\s*>/) {

      # store previous sequence if available
      if ($name ne "") {
	push @{$self->seqname}, $name;
	$self->seqdata->{$name} = $data;
	$data = ""; # reset data
      }

      # store new name
      unless (/^\s*>\s*(\S+)/) { croak "Couldn't parse out sequence name from line: $_\n"; }
      $name = $1;
      if ($remove_chr) { $name = $self->strip_leading_chr ($name); }
    }

    # else store data
    else {
      $data .= $_;
      $data =~ s/\s//g; # remove whitespace
    }

  }
  close $str_handle;

  # catch last sequence
  push @{$self->seqname}, $name;
  $self->seqdata->{$name} = $data;
  
  return $self;

}

=head2 detect_CLUSTAL

    detect_CLUSTAL ($filename);

Detects whether the alignment is CLUSTAL format.

=cut
sub detect_CLUSTAL {
  my ($class, $filename) = @_;

  local *FILE;
  open FILE, "<$filename" or croak "Couldn't open '$filename'.\n";
  while (<FILE>) {
    if (/^\s*CLUSTAL/) { return 1; }
  }
  close FILE;

  return 0;
}

=head2 from_filehandle_CLUSTAL

    $self->from_filehandle_CLUSTAL ($filehandle);

Creates a new Stockholm object and reads it from a filehandle.

=cut
sub from_filehandle_CLUSTAL {
  my ($class, $filehandle, $remove_chr) = @_;
  
  my $self = $class->new;

  # by default don't remove the leading 'chr'
  if (!defined $remove_chr) { $remove_chr = 0; }

  local $_;
  while (<$filehandle>) {
    chomp;

    # skip conservation lines with * : .
    next if (/^(\s*[\*:\.]*\s*)+$/);

    # sequence data
    my @a = split;
    next unless @a == 2;
    my ($seq, $data) = @a;
    if ($remove_chr) { $seq = $self->strip_leading_chr ($seq); }

    if (defined $self->seqdata->{$seq}) {
      $self->seqdata->{$seq} .= $data;
    } else {
      push @{$self->seqname}, $seq;
      $self->seqdata->{$seq} = $data;
    }

  }

  return $self;
  
}

=head2 detect_MSF

    detect_MSF ($filename);

Detects whether the alignment is MSF format.
Looks for the special line containing
"MSF:", "Type:" , and "Check:" and ending with two dots.
Specification taken from:
http://biobug.life.nthu.edu.tw/predictprotein/Dexa/optin_msfDes.html

=cut
sub detect_MSF {
  my ($class, $filename) = @_;

  local *FILE;
  open FILE, "<$filename" or croak "Couldn't open '$filename'.\n";
  while (<FILE>) {
    # "MSF:", "Type:" , and "Check:" and ending with two dots.
    if (/(MSF:){1}.*(Type:){1}.*(Check:){1}.*(\.\.){1}/i) { return 1; }
  }
  close FILE;

  return 0;
}

=head2 from_filehandle_MSF

    $self->from_filehandle_MSF ($filehandle);

Creates a new Stockholm object and reads it from a filehandle.

=cut
sub from_filehandle_MSF {
  my ($class, $filehandle, $remove_chr) = @_;
  
  my $self = $class->new;

  # by default don't remove the leading 'chr'
  if (!defined $remove_chr) { $remove_chr = 0; }

  my $in_body = 0;

  local $_;
  while (<$filehandle>) {
    chomp;

    # alignment header
    if (!$in_body) {

      # if a sequence descriptor in the alignment description,
      # store the sequence name and initialize the sequence data
      if (/^\s*Name:\s*(\S+)\s*/) {
	my $name = $1;
	if ($remove_chr) { $name = $self->strip_leading_chr ($name); }
	push @{$self->seqname}, $name;
	if (defined $self->seqdata->{$name}) {
	  croak "Duplicate sequence name '$name' in MSF alignment descriptor; I'm quitting.\n";
	}
	$self->seqdata->{$name} = "";
      }

      # termination of header list
      elsif (/^\s*\/\/\s*$/) { $in_body = 1; }

    }

    # alignment body 
    else {

      # if a sequence line
      if (/^\s*(\S+)\s+(\S+.*)+\s*$/) {  # allows empty sequence lines
	my $name = $1;
	if ($remove_chr) { $name = $self->strip_leading_chr ($name); }
	if (!defined $self->seqdata->{$name}) {
	  warn "Skipping line because I don't recognize the sequence '$name':\n$_\n";
	}
	$self->seqdata->{$name} .= $2;
	$self->seqdata->{$name} =~ s/\s//g; # remove whitespace
      }

    }

  }

  return $self;
  
}

=head2 detect_Stockholm

    detect_Stockholm ($filename);

Detects whether the alignment is Stockholm format.
Doesn't allow empty sequence lines.

=cut
sub detect_Stockholm {
  my ($class, $filename) = @_;

  local *FILE;
  open FILE, "<$filename" or croak "Couldn't open '$filename'.\n";

  while (<FILE>) {

    # look for Stockholm header
    if (/\# STOCKHOLM/) { return 1; }

    # look for Stockholm comment lines
    elsif (/^\s*\#=[GF|GS|GC|GR]/) { return 1; }

    # skip unrecognized comment lines
    elsif (/^\s*\#.*$/) { next; }

    # skip alignment separators
    elsif (/^\s*\/\/\s*$/) { next; }

    # confirm that format of sequence lines is Stockholm-compatible
    elsif (/^\s*(\S+)\s+(\S+)\s*$/) { next; } # don't allow empty sequence lines

    # else declare it not Stockholm format
    else { return 0; }

  }
  close FILE;

  return 1;
}


=head2 from_filehandle_Stockholm

    my $stock = Stockholm->from_filehandle_Stockholm ($filehandle);

Creates a new Stockholm object and reads it from a filehandle.

=cut
sub from_filehandle_Stockholm {
  my ($class, $filehandle, $remove_chr) = @_;

  my $self = $class->new;

  # by default don't remove the leading 'chr'
  if (!defined $remove_chr) { $remove_chr = 0; }

  local $_;
  while (<$filehandle>) {
    last if $self->parse_input_line ($_, $remove_chr);
  }
  carp "Warning: alignment is not flush" unless $self->is_flush();

  return $self;
}

=head2 to_file

    $stock->to_file ($filename);

Writes a Stockholm object to a file.

=cut
sub to_file {
  my ($self, $filename, $maxcols) = @_;

  local *FILE;
  open FILE, ">$filename" or croak "Couldn't open '$filename' for writing: $!";
  print FILE $self->to_string ($maxcols);
  close FILE or croak "Couldn't close '$filename': $!";;

  return $filename;
}

=head2 to_string_FASTA

   print $fasta->to_string_FASTA();
   print $fasta->to_string_FASTA (0);

Returns the object as a multi-FASTA or ungapped FASTA formatted string.

=cut
sub to_string_FASTA {
  my ($self, $ungapped) = @_;
  
  my $fasta;

  # write MFA rather than ungapped FASTA by default
  if (!defined $ungapped) { $ungapped = 0; }
  
  foreach my $seq (@{$self->seqname}) {
    $fasta .= ">$seq\n";
    my $newseq = $self->seqdata->{$seq};
    if ($ungapped) { $newseq =~ s/[$gapChars]//g; } # remove gaps if requested
    $fasta .= $newseq . "\n";
  }

  return $fasta;

}

=head2 to_string

    print $stock->to_string ($maxcols)
    print $stock->to_string ($maxcols, ARG1=>VAL1, ARG2=>VAL2 ...)
    print $stock->to_string (MAXCOLS=>$maxcols, ARG1=>VAL1, ARG2=>VAL2 ...)

Returns the object as a Stockholm-formatted string.

ARGs can include...
        MAXCOLS    -- limit maximum number of columns (can also be specified as a single scalar arg)
        NOSEQDATA  -- don\'t print any sequence data (can be used this to compress output)

=cut
sub to_string {
  my ($self, @args) = @_;

  my (%args, $maxcols);
  if (@args % 2 == 1) {
    $maxcols = shift @args;
    %args = @args;
  } else {
    %args = @args;
    $maxcols = $args{'MAXCOLS'};
  }
  $maxcols = 80 unless defined $maxcols; # default 80-column output

  # init output array
  my @out;
  push @out, "# STOCKHOLM 1.0";

  # determine alignment columns, legend columns & effective columns per line
  my $acols = $self->columns();
  my $lcols = $self->lcols();
  my $colstep = $maxcols < 1 ? $acols : $maxcols - $lcols - 1;
  $colstep = $maxcols if $colstep < 1; # protect against negative and 0 colstep...

  # GF lines
  # check for gforder (insane, fragile Stockholm line ordering strikes again)
  if (@{$self->gforder} == map { (@$_) } values %{$self->gf}) { # gforder same number of lines as #=GF block?
    my %gfCursor = map (($_ => 0), keys %{$self->gf});
    foreach my $feature (@{$self->gforder}) {
      push @out, $self->prettify ($lcols, "#=GF $feature", $self->gf_($feature)->[$gfCursor{$feature}++]);
    }
  } else {
    @{$self->gforder} = ();	# gforder is useless, so flush it
    foreach my $feature (sort { $a cmp $b } keys %{$self->gf}) {
      push @out, $self->prettify ($lcols, "#=GF $feature", @{$self->gf_($feature)});
    }
  }

  # GS lines
  my @gs_seqname = @{$self->seqname};
  my %gs_seqname_hash = map (($_=>1), grep (!exists $self->seqdata->{$_}, map (keys(%{$self->gs_($_)}), keys %{$self->gs})));
  push @gs_seqname, keys %gs_seqname_hash;
  foreach my $feature (sort { $a cmp $b } keys %{$self->gs}) {
    my $hash = $self->gs_($feature);
    foreach my $seqname (grep (exists($hash->{$_}), @gs_seqname)) {
      push @out, $self->prettify ($lcols, "#=GS $seqname $feature", @{$hash->{$seqname}});
    }
  }

  my @gcfeat = sort { $a cmp $b } keys %{$self->gc};
  my @gr_seqname = @{$self->seqname};
  my %gr_seqname_hash = map (($_=>1), grep (!exists $self->seqdata->{$_}, map (keys(%{$self->gr_($_)}), keys %{$self->gr})));
  push @gr_seqname, keys %gr_seqname_hash;

  # Loop over columns
  for (my $col = 0; $col < $acols; $col += $colstep) {

    # GC lines
    foreach my $feature (@gcfeat) {
      push @out, $self->prettify ($lcols, "#=GC $feature",
				  substr ($self->gc_($feature), $col, $colstep));
    }
    for (my $i = 0; $i < @gr_seqname; ++$i) {
      my $seqname = $gr_seqname[$i];
      # Sequences
      #	    warn "Writing cols $col+$colstep for $seqname";
      push @out, $self->prettify ($lcols, $seqname,
				  substr ($self->seqdata->{$seqname}, $col, $colstep))
	if exists $self->seqdata->{$seqname}
	  && !$args{'NOSEQDATA'};
      # GR lines
      foreach my $feature (grep (exists ($self->gr->{$_}->{$seqname}), keys %{$self->gr})) {
	#		warn "Writing cols $col+$colstep for $seqname $feature";
	push @out, $self->prettify ($lcols, "#=GR $seqname $feature",
				    substr ($self->gr_($feature)->{$seqname}, $col, $colstep));
      }
    }
    push @out, "";
  }

  # alignment separator
  push @out, "//";

  # convert output array to string & return
  return join ("", map ("$_\n", @out));
}

=head2 copy

    my $newStock = $stock->copy();

Does a deep-copy, duplicating all information.

=cut
sub copy {
  my ($self) = @_;

  my $stock = Stockholm->new;

  # Sequence names & data
  @{$stock->seqname} = @{$self->seqname};
  %{$stock->seqdata} = %{$self->seqdata};

  #=GF
  while (my ($feature, $arrayRef) = each %{$self->gf}) {
    $stock->gf->{$feature} = [@$arrayRef];
  }
  @{$stock->gforder} = @{$self->gforder};

  #=GC
  while (my ($feature, $string) = each %{$self->gc}) {
    $stock->gc->{$feature} = $string;
  }

  #=GR
  while (my ($feature, $seqHash) = each %{$self->gr}) {
    $stock->gr->{$feature} = {%$seqHash};
  }

  #=GS
  while (my ($feature, $seqHash) = each %{$self->gs}) {
    $stock->gs->{$feature} = {%$seqHash};
  }

  # Return
  return $stock;
}

=head2 columns

    my $cols = $stock->columns()

Returns the number of columns in this alignment.
The value is recalculated each time because the
sequence or annotation data may have changed since
the method was last called.

=cut
sub columns {
  my ($self) = @_;

  return max (map (length ($_),
		   values (%{$self->seqdata}),
		   values (%{$self->gc}),
		   map (values (%$_), values (%{$self->gr}))));

}

=head2 sequences

    my $cols = $stock->sequences()

Returns the number of sequences (rows) in this alignment.

=cut
sub sequences {
  my ($self) = @_;
  return @{$self->seqname} + 0;
}

=head2 seqname

    my $rowName = $stock->seqname->[$rowIndex];

Returns a reference to an array of sequence names.

=head2 seqdata

    my $row = $stock->seqdata->{$rowName};

Returns a reference to a hash of alignment rows, keyed by sequence name.

=head2 gf

    my @gf = @{$stock->gf->{FEATURE}};
    my @gf = @{$stock->gf_FEATURE};

Returns a reference to an array of all the lines beginning '#=GF FEATURE ...'

=head2 gc

    my $gc = $stock->gc->{FEATURE};
    my $gc = $stock->gc_FEATURE;

Returns the line beginning '#=GC FEATURE ...'

=head2 gs

    my @gs = @{$stock->gs->{FEATURE}->{SEQNAME}};
    my @gs = @{$stock->gs_FEATURE->{SEQNAME}};

Returns a reference to an array of all the lines beginning '#=GS SEQNAME FEATURE ...'

=head2 gr

    my $gr = $stock->gr->{FEATURE}->{SEQNAME};
    my $gr = $stock->gr_FEATURE->{SEQNAME};

Returns the line beginning '#=GR SEQNAME FEATURE ...'

Catches all methods by default.

=cut
sub AUTOLOAD {
  my ($self, @args) = @_;
  my $sub = our $AUTOLOAD;
  $sub =~ s/.*:://;

  # check for DESTROY
  return if $sub eq "DESTROY";

  # check for GF, GC, GS, GR tag_feature accessors
  # e.g. $self->GF_ID
  # can also use $self->GF_('ID')
  if ($sub =~ /^(gf|gc|gs|gr)_(\S*)$/i) {
    my ($tag, $feature) = ($1, $2);
    $tag = lc $tag;
    $feature = shift(@args) unless length $feature;
    my $hash = $self->{$tag};
    if (@args) {
      $hash->{$feature} = shift(@args);	# TODO: check ref-type of arg
    }
    cluck "Warning: ignoring extra arguments to ${tag}_$feature"
      if @args > 0;
    if (!defined $hash->{$feature}) {
      $hash->{$feature} = 
	$tag eq "gf" ? [] :
	  $tag eq "gc" ? "" :
	    $tag eq "gs" ? {} :
	      $tag eq "gr" ? {} :
		croak "Unreachable";
    }
    return $hash->{$feature};
  }

  # check for ordinary accessors
  if (exists $self->{$sub}) {
    croak "Usage: $sub() or $sub(newValue)" if @args > 1;
    return
      @args
	? $self->{$sub} = $args[0]
	  : $self->{$sub};
  }

  # croak
  croak "Unsupported method: $sub";
}


=head2 prettify

    prettify ($lcols, $legend, @data);

Format a line.

=cut
sub prettify {
  my ($self, $lcols, $legend, @data) = @_;

  # This horribly inefficient/redundant series of transformations comes out with something I like (IH, 7/24/07)
  # Trim it down? pah! Like I have nothing better to do
  $legend = sprintf ("% ${lcols}s", $legend);
  $legend =~ s/^(\s+)(\#=\S\S)(\s+\S+)$/$2$1$3/;
  $legend =~ s/^(\s+)(\#=\S\S\s+\S+)(\s+\S+)$/$2$1$3/;
  $legend =~ s/^(\s\s\s\s\s)(\s+)([^\#]\S+)/$1$3$2/;

  return map ("$legend $_", @data);
}

=head2 lcols

    lcols();

Get legend width.
Subtract this from maxcols - 1 to get the number of columns
available for sequence display.

=cut
sub lcols {
  my ($self) = @_;
  my $lcols = max ($self->maxNameLen,
		   map(length("#=GF $_"), keys(%{$self->gf})),
		   map(length("#=GC $_"), keys(%{$self->gc})));
  while (my ($gr_key, $gr_hash) = each %{$self->gr}) {
    $lcols = max ($lcols, map(length("#=GR $gr_key $_"), keys(%$gr_hash)));
  }
  while (my ($gs_key, $gs_hash) = each %{$self->gs}) {
    $lcols = max ($lcols, map(length("#=GS $gs_key $_"), keys(%$gs_hash)));
  }
  return $lcols;
}

=head2 add_gf

    $stock->add_gf (FEATURE, $line)
    $stock->add_gf (FEATURE, $line1, $line2, ...)
    $stock->add_gf (FEATURE, @lines)
    $stock->add_gf (FEATURE, "$line1\n$line2\n$...")

Add '#=GF FEATURE' annotation, preserving the order of the '#=GF' lines.

Each list entry gets printed on its own line.

=cut
sub add_gf {
  my ($self, $feature, @data) = @_;
  unless (@data > 0) {
    carp "Missing parameters to 'add_gf'; nothing will be done";
    return;
  }

  foreach my $line (map {$_ eq '' ? '' : split ("\n", $_, -1)} @data) {
    push @{$self->gf_($feature)}, $line;
    push @{$self->gforder}, $feature;
  }
}

=head2 set_gf

    $stock->set_gf (FEATURE, $line)
    $stock->set_gf (FEATURE, $line1, $line2, ...)
    $stock->set_gf (FEATURE, @lines)
    $stock->set_gf (FEATURE, "$line1\n$line2\n$...")

Set a '#=GF FEATURE ...' annotation.
The difference between this and 'add_gf' is that here, if this sequence and
feature have an annotation already, it will be overwritten instead of appended to.

As with add_gf, each list entry gets printed on its own line.

=cut
sub set_gf {
  my ($self, $feature, @data) = @_;

  unless (defined $self->{gf}->{$feature}) {
    carp "No #=GF feature '$feature' to set, creating new";
    $self->add_gf ($feature, @data);
    return;
  }
  @data = map {$_ eq '' ? '' : split ("\n", $_, -1)} @data;

  ### go to insane lengths to preserve original line ordering

  my %gfCursor = map {$_ => 0} keys %{$self->{gf}};
  my (%new_gfFeature, @new_gforder);

  foreach my $curFeat (@{$self->{gforder}}) {
    $new_gfFeature{$curFeat} = [] unless defined $new_gfFeature{$curFeat};

    if ($curFeat eq $feature) {
      if (@data) {
	push (@{$new_gfFeature{$curFeat}}, shift @data);
	push (@new_gforder, $curFeat);
      }
    } else {
      push (@{$new_gfFeature{$curFeat}},
	    $self->{gf}->{$curFeat}->[$gfCursor{$curFeat}++]);
      push (@new_gforder, $curFeat);
    }
  }

  # still have lines left, insert them after last occurence of feature tag
  if (@data) {
    my $insertAt = $#new_gforder;
    $insertAt-- while $new_gforder[$insertAt] ne $feature;
    splice (@new_gforder, ++$insertAt, 0, map {$feature} @data);
    push (@{$new_gfFeature{$feature}}, @data);
  }

  $self->{gf} = {%new_gfFeature};
  $self->{gforder} = [@new_gforder];
}

=head2 get_gf

    my $gf = $stock->get_gf (FEATURE)

Get a #=GF annotation.
Concatenates and returns all lines beginning '#=GF FEATURE ...'

=cut
sub get_gf {
  my ($self, $feature) = @_;
  if (defined $self->{gf}->{$feature}) {
    return join ("\n", @{$self->{gf}->{$feature}});
  } else {
    #carp "No #=GF feature '$feature' in Stockholm object, returning undef";
    return undef;
  }
}

=head2 add_gs

    $stock->add_gs (SEQNAME, FEATURE, $line)
    $stock->add_gs (SEQNAME, FEATURE, $line1, $line2, ...)
    $stock->add_gs (SEQNAME, FEATURE, @lines)
    $stock->add_gs (SEQNAME, FEATURE, "$line1\n$line2\n$...")

Add a GS annotation: '#=GS SEQNAME FEATURE ...'

If such a line already exists for this sequence and feature, append to it.

Each list entry gets printed on its own line.

=cut
sub add_gs {
  my ($self, $seq, $feature, @text) = @_;
  unless (@text > 0) {
      carp "Missing parameters to 'add_gs' method; nothing will be done";
      return;
  }
  @text = map {$_ eq '' ? '' : split ("\n", $_, -1)} @text;
  push (@{$self->{gs}->{$feature}->{$seq}}, @text);
}

=head2 set_gs

    $stock->set_gs (SEQNAME, FEATURE, $line)
    $stock->set_gs (SEQNAME, FEATURE, $line1, $line2, ...)
    $stock->set_gs (SEQNAME, FEATURE, @lines)
    $stock->set_gs (SEQNAME, FEATURE, "$line1\n$line2\n$...")

Set a GS annotation: '#=GS SEQNAME FEATURE ...'

The difference between this and 'add_gs' is that here, if this sequence and
feature have an annotation already, it will be overwritten instead of appended to.

Each list entry gets printed on its own line.

=cut
sub set_gs {
  my ($self, $seq, $feature, @text) = @_;
  unless (@text > 0) {
      carp "Missing parameters to 'set_gs' method; nothing will be done";
      return;
  }
  @text = map {$_ eq '' ? '' : split ("\n", $_, -1)} @text;
  $self->{gs}->{$feature}->{$seq} = [@text];
}

=head2 get_gs

    my $gs = $stock->get_gs (SEQNAME, FEATURE)

Get a #=GS annotation.
Concatenates and returns all lines beginning '#=GS SEQNAME FEATURE ...'

=cut
sub get_gs {
  my ($self, $seq, $feature) = @_;
  unless (defined $feature) {
      carp "Missing parameters to 'get_gs' method; nothing will be done";
      return undef;
  }
  if (defined $self->{gs}->{$feature}->{$seq}) {
      return join ("\n", @{$self->{gs}->{$feature}->{$seq}});
  }
  else {
      #carp "Annotation for sequence $seq, feature $feature not found; returning empty string";
      return undef;
  }
}

=head2 drop_columns

    my $success = $stock->drop_columns (@columns)

Drops a set of columns from the alignment and the GC and GR annotations.  Note
that this may break annotations where the column annotations are not
independent (e.g. you might drop one base in an RNA base pair, but not the
other).  The caller is responsible for making sure the set of columns passed
in does not mess up the annoation.

Arguments: a list of scalar, zero-based column indices to drop.

Returns success of operation (1 for success, 0 for failure)

=cut
sub drop_columns {
  my ($self, @columns) = @_;

  unless ($self->is_flush()) {
    carp "File is not flush; cannot (should not) drop columns";
    return 0;
  }

  # drop from sequences and per-column annotations
  foreach my $seq_or_annot ($self->{seqdata}, $self->{gc})
    {
      # for each sequence or annotation string
      foreach my $key (keys %$seq_or_annot)
	{
	  # convert string to array form
	  my @data = split (//, $seq_or_annot->{$key});
	  foreach my $col (@columns)
	    {
	      if ( ($col < 0) or ($col >= @data) ) {
		carp "Column $col is out of range [1, ", scalar(@data),
		  "] in $key; can't drop, skipping";
	      }
	      else {
		# blank out the column that should get dropped
		$data[$col] = '';
	      }
	    }
	  # re-assemble back to string
	  $seq_or_annot->{$key} = join ('', @data);
	}
    }

  # drop from per-column/per-sequence annotations
  foreach my $feat_key (keys %{$self->{gr}})
    {
      foreach my $seq_key (keys %{$self->{gr}->{$feat_key}})
	{
	  my @data = split (//, $self->{gr}->{$feat_key}->{$seq_key});
	  foreach my $col (@columns)
	    {
	      if ( ($col < 0) or ($col >= @data) ) {
		carp "Column $col is out of range [1, ", scalar(@data),
		  "] in $feat_key, $seq_key; can't drop, skipping";
	      }
	      else {
		# blank out the column that should get dropped
		$data[$col] = '';
	      }
	    }
	  # re-assemble back to string
	  $self->{gr}->{$feat_key}->{$seq_key} = join ('', @data);
	}
    }

  return 1;
}

# parse line of Stockholm format file
# returns true if it finds a separator
sub parse_input_line {
  my ($self, $line, $remove_chr) = @_;

  # by default don't remove the leading 'chr'
  if (!defined $remove_chr) { $remove_chr = 0; }

  $line = $_ unless defined $line;

  # "#=GF [feature] [data]" line
  if ($line =~ /^\s*\#=GF\s+(\S+)\s+(\S.*)\s*$/) { 
    my ($feature, $data) = ($1, $2);
    $self->add_gf ($feature, $data); # preserve order of #=GF lines, for crazy context-sensitive Stockholm semantics
  }

  # "#=GC [feature] [data]" line
  elsif ($line =~ /^\s*\#=GC\s+(\S+)\s+(\S+)\s*$/) { 
    my ($feature, $data) = ($1, $2, $3);
    $self->gc->{$feature} .= $data; 
  }

  # "#=GS [seqname] [feature] [data]" line
  elsif ($line =~ /^\s*\#=GS\s+(\S+)\s+(\S+)\s+(\S.*)\s*$/) { 
    my ($seqname, $feature, $data) = ($1, $2, $3);
    if ($remove_chr) { $seqname = $self->strip_leading_chr ($seqname); }
    my $gs = $self->gs_($feature);
    $gs->{$seqname} = [] unless exists $gs->{$seqname};
    push @{$gs->{$seqname}}, $data;
  }

  # "#=GR [seqname] [feature] [data]" line
  elsif ($line =~ /^\s*\#=GR\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) { 
    my ($seqname, $feature, $data) = ($1, $2, $3);
    if ($remove_chr) { $seqname = $self->strip_leading_chr ($seqname); }
    my $gr = $self->gr_($feature);
    $gr->{$seqname} = "" unless exists $gr->{$seqname};
    $gr->{$seqname} .= $data; 
  }

  # Unrecognised "#" line
  elsif ($line =~ /^\s*\#.*$/) {
    #	warn "Comment line $line";
  }

  # Alignment separator: return true to indicate loop exit
  elsif ($line =~ /^\s*\/\/\s*$/) { 
    return 1;
  }

  # Sequence line
  elsif ($line =~ /^\s*(\S+)\s*(\S*)\s*$/) { # allows empty sequence lines
    my ($seqname, $data) = ($1, $2);
    if ($remove_chr) { $seqname = $self->strip_leading_chr ($seqname); }
    unless (defined $self->seqdata->{$seqname}) {
      $self->seqdata->{$seqname} = "";
      push @{$self->seqname}, $seqname;	# preserve sequence order, for tidiness
    }
    $self->seqdata->{$seqname} .= $data;
  } elsif ($line =~ /\S/) { 
    carp "Ignoring line: $_";
  }

  # This line wasn't a alignment separator: return false
  return 0;
}

=head2 empty

    my $isEmpty = $stock->empty

Returns true if the alignment is empty (i.e. there are no sequences).

=cut
sub empty {
  my ($self) = @_;
  return @{$self->seqname} == 0;
}

=head2 ungapped_length

   my $seqlength = $stock->ungapped_length (SEQNAME)

Calculates the length of the ungapped sequence.

=cut
sub ungapped_length {
  my ($self, $name) = @_;

  if (!defined $self->seqdata->{$name}) { croak "No sequence exists for '$name':\n", $self->to_string; }

  my $numgaps = $self->seqdata->{$name} =~ tr/-._//;
  return (length ($self->seqdata->{$name}) - $numgaps);

}

=head2 subseq

    my $subseq = $stock->subseq (SEQNAME, $startPos)
    my $subseq = $stock->subseq (SEQNAME, $startPos, $length)

Returns a subsequence of the named alignment row.

Note that the coordinates are with respect to the alignment columns (starting at zero), NOT sequence positions.
That is, gaps are counted.

If $length is omitted, this call returns the subsequence from $startPos to the end of the row.

=cut
sub subseq {
  my ($self, $name, $start, $len) = @_;

  if (!defined $self->seqdata->{$name}) { croak "No sequence exists for '$name'."; }

  if (!defined $len) { carp "Setting undefined length to 0."; $len = 0; }

  my $seqlength = length $self->seqdata->{$name};

  # check start position sane
  if (!defined $start) { croak "Start position undefined."; }
  elsif (($start > $seqlength) or ($start < 0)) {
    carp
      "Subsequence start position $start is outside of sequence bounds ",
	"[0, " . $seqlength-1 . "]; returning empty string for subsequence of $name\n";
    return '';
  }

  return substr ($self->seqdata->{$name}, $start, $len);
}

=head2 get_column

    my $column = $stock->get_column ($colNum)

Extracts the specified column and returns it as a string.

The column number $colNum uses zero-based indexing.

=cut
sub get_column {
  my ($self, $col_num) = @_;
  my $col_as_string;

  unless ($self->is_flush()) {
    carp "File is not flush; cannot (should not) get column";
    return 0;
  }

  my $align_length = $self->columns();

  if ( ($col_num >= $align_length) or ($col_num < 0) ) {
    carp
      "Column index $col_num is outside of alignment bounds ",
	"[1, $align_length]; returning empty string for column $col_num\n";
    return '';
  }
  else {
    foreach my $seq (@{$self->{seqname}}) {
      $col_as_string .= substr ($self->{seqdata}->{$seq}, $col_num, 1);
    }
    return $col_as_string;
  }
}

=head2 map_align_to_seq_coords

    my ($start, $end) = $stock->map_align_to_seq_coords ($colStart, $colEnd, $seqName);
    croak "Failed" unless defined ($start);

    my ($start, $end) = $stock->map_align_to_seq_coords ($colStart, $colEnd, $seqName, $seqStart);
    croak "Failed" unless defined ($start);

Given a range of column indices ($colStart, $colEnd), returns the coordinates
of sequence $seqName contained within that range, accounting for gaps.

Returns start and end coordinates as a 2-ple list, or returns the empty list
if there is no sequence in the given range (i.e. range contains all gaps in
$seqName).

You can provide $seqStart, which is the coordinate of the first non-gap
character in $seqName.  If you don\'t, the default is that the sequence starts
at 0 (0-based indexing).

=cut
sub map_align_to_seq_coords {
  my ($self, $colStart, $colEnd, $seqName, $seqStart) = @_;
  $seqStart = 0 unless defined $seqStart;

  my $seq = $self->subseq ($seqName, $colStart, ($colEnd - $colStart + 1));

  if (my @seq = $seq =~ /[^$gapChars]/g)
  {
    my $rangeStart;

    if ($colStart == 0) {
      # special case: range starts in the first column
      $rangeStart = $seqStart;
    }
    else {
      my $prefix = $self->subseq ($seqName, 0, $colStart);
      $rangeStart = $seqStart + length (join ('', ($prefix =~ /[^$gapChars]/g)));
    }

    return ($rangeStart, ($rangeStart + @seq - 1));
  }
  else {
    # no sequence in desired range (i.e. all gaps)
    return ();
  }
}

=head2 map_seq_to_align_coords

    my ($start, $end) = $stock->map_seq_to_align_coords ($seq, $start, $end);

Given a range of sequence position indices ($start, $end),
returns the alignment coordinates corresponding to that range.
Returns undef if can't map coordinates.

Returns start and end coordinates as a duple.

=cut
sub map_seq_to_align_coords {
  my ($self, $seq, $start, $end) = @_;

  # sanity check
  if (!defined $self->seqdata->{$seq}) { croak "Can't find sequence '$seq'.\n"; }
  if ($end < $start) { croak "Invalid coordinates.\n"; }

  # do the mapping with a linear search
  my ($s, $e);    # alignment coordinates
  my $seqpos = 0; # position within the sequence

  # NB: This operation can be quite slow,
  # although the time hit is worth it if the desired
  # coordinates are deep into the alignment
  # (substr to pull out each character becomes expensive)
  # unpack seems to be slightly faster than split
  my @data = unpack 'a' x  length $self->seqdata->{$seq}, $self->seqdata->{$seq};

  my $cols = $self->columns();
  for (my $alignpos = 0; $alignpos < $cols; ++$alignpos) {

    # keep going if it's a gap
    if ($data[$alignpos] =~ /[$gapChars]/) { next; }

    # else see if we've reached a boundary
    if ($seqpos == $start) { $s = $alignpos; }
    if ($seqpos == $end) { $e = $alignpos; last; }

    # increment coords
    ++$seqpos;

  }

  return ($s, $e);
}

=head2 build_seq_to_align_coords_map

    my $map = $stock->build_seq_to_align_coords_map ($seq);

Returns a map from sequence to alignment coordinates
for a particular sequence as an array reference.

=cut
sub build_seq_to_align_coords_map {
  my ($self, $seq, $verbose) = @_;

  if ($verbose) {
    print STDERR "Building coordinate index for '$seq'...";
  }

  # sanity check
  if (!defined $self->seqdata->{$seq}) { croak "Can't find sequence '$seq'.\n"; }

  my @map;

  # do the mapping with a linear search
  my ($s, $e);    # alignment coordinates
  my $seqpos = 0; # position within the sequence

  # assemble sequence data as character arrays for fast access
  # (avoid having to use substr, which is slow)
  # unpack seems to be slightly faster than split
  my @data = unpack 'a' x  length $self->seqdata->{$seq}, $self->seqdata->{$seq};

  my $cols = $self->columns();
  for (my $alignpos = 0; $alignpos < $cols; ++$alignpos) {

    # keep going if it's a gap
    if ($data[$alignpos] =~ /[$gapChars]/) { next; }

    # else store the record
    $map[$seqpos++] = $alignpos;

  }

  if ($verbose) {
    print STDERR "done.\n";
  }

  return \@map;

}


=head2 is_gap

    $stock->is_gap (SEQNAME, $colNum)

Return true if a given (SEQNAME,column) coordinate is a gap.

=cut
sub is_gap {
  my ($self, $row_name, $col) = @_;
  my $c = substr ($self->seqdata->{$row_name}, $col, 1);
  return $c =~ /[$gapChars]/;
}

=head2 is_allgaps_col

    $stock->is_allgaps_col ($col)

Returns true if a given column is all gaps.

=cut
sub is_allgaps_col {
  my ($self, $c) = @_;
  my $is = 1;
  foreach my $seqname (@{$self->seqname}) {
    if (!$self->is_gap ($seqname, $c)) {
      $is = 0;
    }
  }
  return $is;
}

=head2 drop_allgaps_columns

    $stock->drop_allgaps_columns()

Drops columns which are all gaps from an alignment.

=cut
sub drop_allgaps_columns {
  my ($self) = @_;

  my $cols = $self->columns();
  my @drop; # columns to drop
  for (my $c = 0; $c < $cols; ++$c) {
    if ($self->is_allgaps_col ($c)) {
      push @drop, $c;
    }
  }

  $self->drop_columns (@drop);

}

=head2 subalign

    my $substock = $stock->subalign ($start, $len)
    my $substock = $stock->subalign ($start, $len, $strand)

Subalignment accessor.
Returns a Stockholm object representing a sub-alignment of the current alignment,
starting at column $start (zero-based) and containing $len columns.

If $strand is supplied and is equal to '-', then the sub-alignment will be reverse-complemented.

=cut
sub subalign {
  my ($self, $start, $len, $strand) = @_;

  # check strand sane
  if (!defined $strand) {
    $strand = '+';
  } elsif (($strand ne '-') && ($strand ne '+')) {
    warn "Unrecognized strand '$strand'; using positive strand.\n";
    $strand = '+';
  }

  # init subalign
  my $subalign = Stockholm->new;

  # sequence data
  foreach my $seqname (@{$self->seqname}) {
    $subalign->seqdata->{$seqname} = substr ($self->seqdata->{$seqname}, $start, $len);
    if ($strand eq '-') {
      $subalign->seqdata->{$seqname} = revcomp ($subalign->seqdata->{$seqname});
    }
  }
  @{$subalign->seqname} = @{$self->seqname};

  # GF lines
  while (my ($feature, $arrayRef) = each %{$self->gf}) {
    $subalign->gf->{$feature} = [@$arrayRef];
  }
  @{$subalign->gforder} = @{$self->gforder};

  # GC lines
  foreach my $tag (keys %{$self->gc}) {
    $subalign->gc->{$tag} = substr ($self->gc->{$tag}, $start, $len);
    if ($strand eq '-') {
      $subalign->gc->{$tag} = reverse $subalign->gc->{$tag};
    }
  }

  # GR lines
  foreach my $tag (keys %{$self->gr}) {
    foreach my $seqname (keys %{$self->gr->{$tag}}) {
      $subalign->gr->{$tag}->{$seqname} = substr ($self->gr->{$tag}->{$seqname}, $start, $len);
      if ($strand eq '-') {
	$subalign->gr->{$tag}->{$seqname} = reverse ($subalign->gr->{$tag}->{$seqname});
      }
    }
  }


  #=GS
  while (my ($feature, $seqHash) = each %{$self->gs}) {
    $subalign->gs->{$feature} = {%$seqHash};
  }

  return $subalign;
}

=head2 concatenate

    $stock->concatenate ($stock2)

Concatenates another alignment onto the end of this one.

=cut
sub concatenate {
  my ($self, $stock) = @_;

  # get widths
  my $cols = $self->columns;
  my $catcols = $stock->columns;

  # Sequence data
  foreach my $seqname (@{$stock->seqname}) {
    if (exists $self->seqdata->{$seqname}) {
      $self->seqdata->{$seqname} .= $stock->seqdata->{$seqname};
    } else {
      push @{$self->seqname}, $seqname;
      $self->seqdata->{$seqname} = "." x $cols . $stock->seqdata->{$seqname};
    }
  }
  foreach my $seqname (@{$self->seqname}) {
    unless (exists $stock->seqdata->{$seqname}) {
      $self->seqdata->{$seqname} .= "." x $catcols;
    }
  }

  # GC lines
  foreach my $tag (keys %{$stock->gc}) {
    if (exists $self->gc->{$tag}) {
      $self->gc->{$tag} .= $stock->gc->{$tag};
    } else {
      $self->gc->{$tag} = "." x $cols . $stock->gc->{$tag};
    }
  }
  foreach my $tag (keys %{$self->gc}) {
    unless (exists $stock->gc->{$tag}) {
      $self->gc->{$tag} .= "." x $catcols;
    }
  }

  # GR lines
  foreach my $tag (keys %{$stock->gr}) {
    if (exists $self->gr->{$tag}) {
      foreach my $seqname (keys %{$stock->gr->{$tag}}) {
	if (exists $self->gr->{$tag}->{$seqname}) {
	  $self->gr->{$tag}->{$seqname} .= $stock->gr->{$tag}->{$seqname};
	} else {
	  $self->gr->{$tag}->{$seqname} = "." x $cols . $stock->gr->{$tag}->{$seqname};
	}
      }
    } else {
      $self->gr->{$tag} = {};
      foreach my $seqname (keys %{$stock->gr->{$tag}}) {
	$self->gr->{$tag}->{$seqname} = "." x $cols . $stock->gr->{$tag}->{$seqname};
      }
    }
  }
  foreach my $tag (keys %{$self->gr}) {
    foreach my $seqname (keys %{$self->gr->{$tag}}) {
      unless (exists $stock->gr->{$tag} && exists $stock->gr->{$tag}->{$seqname}) {
	$self->gr->{$tag}->{$seqname} .= "." x $catcols;
      }
    }
  }

  # GF and GS lines
  foreach my $tag (keys %{$stock->gf}) {
    $self->gf->{$tag} = [] unless exists $self->gf->{$tag};
    push @{$self->gf->{$tag}}, @{$stock->gf->{$tag}};
  }

  foreach my $tag (keys %{$stock->gs}) {
    $self->gs->{$tag} = {} unless exists $self->gs->{$tag};
    foreach my $seqname (keys %{$stock->gs->{$tag}}) {
      $self->gs->{$tag}->{$seqname} = [] unless exists $self->gs->{$tag}->{$seqname};
      push @{$self->gs->{$tag}->{$seqname}}, @{$stock->gs->{$tag}->{$seqname}};
    }
  }
}

# add_row method
sub add_row {
    my ($self, $seqname, $rowdata) = @_;
    die "Attempt to add duplicate row" if exists $self->seqdata->{$seqname};
    push @{$self->seqname}, $seqname;
    $self->seqdata->{$seqname} = $rowdata;
    return @{$self->seqname} - 1;  # new row index
}


# Reverse complement a sequence
sub revcomp {
  my ($arg0, $arg1) = @_;
  my $seq;

  if (defined ($arg1)) {  # must have been called externally (as method or package sub)
    $seq = $arg1;
  } else {
    $seq = $arg0;
  }

  $seq =~ tr/acgtuACGTU/tgcaaTGCAA/;
  $seq = reverse $seq;
  return $seq;
}

# Determine maximum sequence name length
sub maxNameLen {
  my ($self) = @_;
  return max (map ( length ($_), @{$self->seqname}));
}

=head2 is_flush

    my $flush = $stock->is_flush()

Returns true if the alignment is "flush", i.e. all sequence and annotation lines have the same length.

=cut
sub is_flush {
  my ($self) = @_;
  my $columns = $self->columns();

  # do stuff that has one key

  map {
    if (length ($self->{seqdata}->{$_}) != $columns) {
      carp "Sequence $_ is not flush!";
      return 0;
    }
  } keys %{$self->{seqdata}};

  map {
    if (length($self->{gc}->{$_}) != $columns) {
      carp "#=GC annotation for feature $_ is not flush!";
      return 0;
    }
  } keys %{$self->{gc}};

  # do stuff that has two keys

  map {
    if (length($_) != $columns) {
      carp "#=GR annotation is not flush! ($_)";
      return 0;
    }
  } map { values %{$self->{gr}->{$_}} } keys %{$self->{gr}};

  return 1;  # if we got this far, everything must be flush
}

=head2 assert_flush

    $stock->assert_flush()

Dies if the alignment isn't "flush."

=cut
sub assert_flush {
  my ($self) = @_;

  unless ($self->is_flush()) {
    die "Alignment not flush (rows are not the same length).\n";
  }
}


# Determine the list maximum value.
sub max  {
  my ($x, @y) = @_;
  foreach my $y (@y) 
  { 
    $x = $y if !defined($x) || (defined($y) && $y > $x)
  }
  return $x;
}

=head2 NH

    $self->NH();
    $self->NH ($newtree);

Get/set the New Hampshire tree.

=cut
sub NH {
  my ($self, $newval) = @_;
  if (defined $newval) {
    @{$self->gf_NH} = ($newval);
  }
  return join ("", @{$self->gf_NH});
}

=head2 gapfraction

   $self->gapfraction()

Calculate fraction of gaps (# gaps / total # of characters)
in the alignment.

=cut
sub gapfraction {
  my ($self) = @_;

  my ($num, $denom) = (0, 0);
  map { $num += $self->ungapped_length ($_); $denom += length ($self->seqdata->{$_}); } keys %{$self->seqdata};

  return ($num / $denom);
}

=head2 percentid_annot

    $self->percentid_annot()

Annotation for #=GF Percent_id line.

=cut
sub percentid_annot {
  my $self = shift;
  return $percentid_annot;
}

=head2 gapfraction_annot

    $self->gapfraction_annot()

Annotation for #=GF Gap_fraction line.

=cut
sub gapfraction_annot {
  my $self = shift;
  return $gapfraction_annot;
}

=head2 percentid

    $self->percentid();
    $self->percentid ($seqs);

Calculate percent id for two or more sequences of the alignment.
Calculates per-column percent id as:
(# of times most common non-gap character in column appears) / (# non-gap characters in column)
Reports the average per-column percent id.

Considers gaps during calculation if requested.

=cut
sub percentid {
  my ($self, $seqs, $countgaps) = @_;

  # must be at least 2 seqs
  # ($seqs = ref to array of sequence names)
  if (!defined $seqs) { $seqs = $self->seqname; }
  unless (@$seqs >= 2) { croak "Must be at least 2 sequences.\n"; }

  # don't count gaps by default
  if (!defined $countgaps) { $countgaps = 0; }

  # assemble sequence data as character arrays for fast access
  my $data;
  foreach my $seq (@$seqs) {
    croak "Sequence '$seq' not defined.\n" unless exists $self->seqdata->{$seq};
    @{$data->{$seq}} = unpack 'a' x  length $self->seqdata->{$seq}, $self->seqdata->{$seq};
  }

  # assert flush
  croak "Alignment is not flush.\n" unless $self->is_flush();

  my $sum = 0;
  my $n = 0;

  for (my $c = 0; $c < $self->columns; ++$c) {
    my %f;
    my $count = 0; # ungapped
    my $countall = 0; # gaps as well
    for (my $i = 0; $i < @$seqs; ++$i) {
      my $char = lc $data->{$seqs->[$i]}->[$c];
      if ($char ne '-' && $char ne '.') {
	++$f{$char};
	++$count;
      }
      ++$countall;
    }
    # unless $countgaps, only count columns with > 1 non-gap character
    if (!$countgaps && ($count > 1)) {
      my @sym = sort { $f{$b} <=> $f{$a} } keys %f;

      my $id;
      # if no character appears more than once, then % id = 0
      if ($f{$sym[0]} == 1) {
	$id = 0;
      }
      # else calculate % id as (# most common character) / (# non-gap characters in column)
      else {
	$id = $f{$sym[0]} / $count;
      }
      $sum += $id;
      ++$n;
    }
    # if $countsgaps, count columns with >= 1 non-gap character
    elsif ($countgaps && ($count >= 1)) {
      my @sym = sort { $f{$b} <=> $f{$a} } keys %f;

      my $id;
      # if no character appears more than once, then % id = 0
      if ($f{$sym[0]} == 1) {
	$id = 0;
      }
      # else calculate % id as (# most common character) / (# total characters in column)
      else {
	$id = $f{$sym[0]} / $countall;
      }
      $sum += $id;
      ++$n;
    }

  }

  if ($n == 0) {
    return 0;
  }
  
  return $sum / $n;
}

=head2 strip_leading_chr

    strip_leading_chr ("chr2R");

Removed the leading 'chr' frequently prepended to chromosome names.

=cut
sub strip_leading_chr {

  my ($class, $str) = @_;

  $str =~ s/^chr//;

  return $str;

}

1
