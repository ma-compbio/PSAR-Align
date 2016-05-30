#!/usr/bin/perl

use strict;
use warnings;

my $sample_dir = shift;

# parse alignment samples
my @files = <$sample_dir/*.fa>;

my @spcs = ();
my %hs_seqs = ();
my $sampleid = 0;
foreach my $file (@files) {
	$sampleid++;
	my ($spc, $seq) = ("", "");
	open(F,$file);
	while(<F>) {
		chomp;
		if (length($_) == 0 || $_ =~ /^#/) { next; }

		if ($_ =~ /^>(\S+)/) {
			if (length($seq) > 0) {
				$hs_seqs{$sampleid}{$spc} = $seq;
				$seq = "";
			}
			$spc = $1;
			if ($sampleid == 1) { push(@spcs, $spc); }
		} else {
			$seq .= $_;	
		}
	}
	close(F);
	if (length($seq) > 0) { 
		$hs_seqs{$sampleid}{$spc} = $seq; 
	}
}

my $total = $sampleid;
my %hs_count = ();

foreach my $sid (sort {$a<=>$b} keys %hs_seqs) {
	for (my $i = 0; $i < $#spcs; $i++) {
		my $spc1 = $spcs[$i];
		my $seq1 = $hs_seqs{$sid}{$spc1};
		my @chars1 = split '', $seq1;	

		for (my $j = $i+1; $j <= $#spcs; $j++) {
			my $spc2 = $spcs[$j];
			my $seq2 = $hs_seqs{$sid}{$spc2};

			# compare seq1 and seq2
			my ($ii, $jj) = (0,0);
			my @chars2 = split '', $seq2;
			for (my $ai = 0; $ai <= $#chars1; $ai++) {
				my ($char1, $char2) = ($chars1[$ai], $chars2[$ai]);
			
				if ($char1 ne "-" && $char2 ne "-") {
					my $cnt = $hs_count{$i}{$j}{$ii}{$jj};
					if (defined($cnt)) {
						$hs_count{$i}{$j}{$ii}{$jj}++;
					} else {
						$hs_count{$i}{$j}{$ii}{$jj} = 1;
					}
				}	

				if ($char1 ne "-") { $ii++; }	
				if ($char2 ne "-") { $jj++; }	
			}	

		}
	}
}

my %hs_gapprobi = ();
my %hs_gapprobj = ();

# match posteriors
print "; match posteriors\n";
foreach my $i (sort {$a<=>$b} keys %hs_count) {
	my $rhsj = $hs_count{$i};
	foreach my $j (sort {$a<=>$b} keys %$rhsj) {
		my $rhsii = $$rhsj{$j};
		foreach my $ii (sort {$a<=>$b} keys %$rhsii) {
			my $rhsjj = $$rhsii{$ii};
			foreach my $jj (sort {$a<=>$b} keys %$rhsjj) {
				my $cnt = $$rhsjj{$jj};
				my $prob = $cnt/$total;
				if ($prob > 0.0) {
					printf "(%d, %d) ~ (%d, %d) => %.6f\n", $i, $ii, $j, $jj, $prob;
				}

				if (defined($hs_gapprobi{$i}{$ii}{$j})) {
					$hs_gapprobi{$i}{$ii}{$j} += $prob;	
				} else {
					$hs_gapprobi{$i}{$ii}{$j} = $prob;	
				}
				
				if (defined($hs_gapprobj{$i}{$j}{$jj})) {
					$hs_gapprobj{$i}{$j}{$jj} += $prob;	
				} else {
					$hs_gapprobj{$i}{$j}{$jj} = $prob;	
				}
			}
		}
	}
}

# gap posteriors
print "\n; gap posteriors\n";
foreach my $i (sort {$a<=>$b} keys %hs_gapprobi) {
	my $rhsii = $hs_gapprobi{$i};
	foreach my $ii (sort {$a<=>$b} keys %$rhsii) {
		my $rhsj = $$rhsii{$ii};
		foreach my $j (sort {$a<=>$b} keys %$rhsj) {
			my $prob = $$rhsj{$j};
			if (1-$prob > 0.0) {
				printf "(%d, %d) ~ (%d, -1) => %.6f\n", $i, $ii, $j, 1-$prob;
			}
		}
	}
}
foreach my $i (sort {$a<=>$b} keys %hs_gapprobj) {
	my $rhsj = $hs_gapprobj{$i};
	foreach my $j (sort {$a<=>$b} keys %$rhsj) {
		my $rhsjj = $$rhsj{$j};
		foreach my $jj (sort {$a<=>$b} keys %$rhsjj) {
			my $prob = $$rhsjj{$jj};
			if (1-$prob > 0.0) {
				printf "(%d, -1) ~ (%d, %d) => %.6f\n", $i, $j, $jj, 1-$prob;
			}
		}
	}
}
